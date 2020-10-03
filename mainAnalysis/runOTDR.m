% FUNCTION for running oTDR for given set of datasets and predictors.
% Called by oTDR.m
% TODO: add docstring
function [Summary,errPerIter,optimGradFinal] = runOTDR(R,codedParams,dataset_times,codedParams4Dataset,...
    orthFlg,dataTensor,dataTensorExcluded,commCond,commCond4Proj,trCountmtx,...
    bitSerialOTDR,bit4Bootstrap,RARand,vect4Angle,spcConstraint,...
    varargin) 
%
% vect4Angle -- N x P matrix specifying the N dimensions of P vectors. All
%       angles between the provided vectors and the discovered
%       non-orthogonalised sRAs (i.e., sRAStar) are measured with
%       associated p-values using the RARand (biased vector) null model.
%       (N dimenions must match the N dimensions provided in R.)

%% set defaults for optional vars

% scalar specifying threshold for goodness-of-fit in fitting parametric
% distribution to empirical null model for purposes of computing p-values.
% When threshold is exceeded, error is thrown.
gofThresh = 0.17;

%% collect optionally provided parameters
warnopts(assignopts(who, varargin));

%% set default values

if bitSerialOTDR
    % It used to be when fitting to serial data (in TDRvarAnalysis.m), the null distributions were not well fit
    % by gamma (the default) and we used empirical instead. However, now we
    % pass random vectors based on the full dataspace and thus the
    % null distributions are better behaved and we use the default (gamma)
    % distribution.
%     nullStatModel = 'empirical';
    nullStatModel = [];    
else
    % leave empty to use default model (gamma)
    nullStatModel = [];
end
    
% logical on whether to plot data and fits of distribution of null angles
% between sRAs
bitPlotRandAngleFit = false;

% logical on whether to compute statistics on the angles between subspaces,
% which takes a great deal of time to generate the null distributions
bitSubspaceStats = false;

% RARand is a matrix of random vectors that is necessary when performing
% serial oTDR, since the input data R does not contain all the data to
% define the dataspace. However, when performing regular oTDR, do not use
% the pre-existing random vectors, as it takes more time to compute the
% variance they explain than to generate them on the fly and compute 1 by 1
% (this is a quirk of varianceIXs.m)
if ~bitSerialOTDR
    RARand = [];
end

% initialize final gradient of optimization error function. When zero,
% optimization converged within max interations. When >0, it did not and
% the final grad of error function is returned.
% optimGradFinal = 0;

%% find PCs based on dataTensor
[T, N, C] = size(dataTensor);
bigA = reshape(permute(dataTensor, [1 3 2]), T*C, N);
[PCs,~,eigenVal] = pca(bigA);

%% generate custom datasets and predictors

RBar = [];
t4SRAAll = {};
nDatasetMax = length(R); % inferred from number of datasets provided.
% loop through datasets, generating data only when predictors
% were specified
for i = 1:nDatasetMax
    if i <= length(codedParams4Dataset) && ~isempty(codedParams4Dataset{i})
        RBar{i} = squeeze(mean(R{i}, 1)).';
        
        % build cell array of times used for each predictor based on times
        % of dataset
%         eval(sprintf('t4SRAAll(end+1:end+length(codedParams4Dataset{i})) = {t%d};',i));
        t4SRAAll(end+1:end+length(codedParams4Dataset{i})) = dataset_times(i);
    else
        RBar{i} = [];
    end
end

% define the datasets. May use 2 dataset (generally with Benefit in dataset
% 1 and shared Choice and ER +/- Benefit in dataset 2), or 3 datasets
% (generally with Benefit in dataset 1, Choice in dataset 2, and ER in
% dataset 3).
% RBar{1} = squeeze(mean(R{1}, 1)).';
% RBar{2} = squeeze(mean(R{2}, 1)).';
% RBar{3} = squeeze(mean(R{3}, 1)).';

% dataTensorC = dataTensor(all_times>=1, :, :);
% for i = 1:size(dataTensorC,1)
%     RRBar{i} = squeeze(dataTensorC(i, :, :))';
% end

% loop through sets of predictors specified for each dataset, building
% matrix of coded parameters for each dataset
CC = [];
for i = 1:nDatasetMax
    if i <= length(codedParams4Dataset)
        CC{i} = codedParams(:,codedParams4Dataset{i});
    else
        CC{i} = [];
    end
end

% pass the benefit predictor only for the first dataset
% CC{1} = codedParams(:,1:3);

% pass the benefit and ER predictor for the first dataset so as to increase
% specificity of the Benefit sRA for the benefit signal. This is not
% helpful in reducing variance explained by ER.
% CC{1} = codedParams(:,[1,3]);

% by passing columns 2 and 3, the objective will solve for choice and
% expected reward simultaneously in a shared choice and ER dataset
% CC{2} = codedParams(:,2:3); 

% by passing columns 1, 2 and 3, the objective will solve for benefit,
% choice and expected reward simultaneously in a shared choice and ER
% dataset
% CC{2} = codedParams(:,1:3);

% pass the ER predictor as a third set of predictors when including a 3rd
% dataset. In this case, the objective will solve for ony ER only in the
% 3rd dataset
% CC{3} = codedParams(:,3); 

% some datasets maybe empty (such as when no predictors were specified for
% the dataset). Loop backwards through datasets and remove empty ones
for i = length(RBar):-1:1
    if isempty(RBar{i})
        % check that coded params are also empty
        if ~isempty(CC{i})
            error('When a dataset is skipped, no coded parameters can be specified for that dataset')
        end
        % eliminate dataset and corresponding coded params
        RBar(i) = [];
        CC(i) = [];
    end
end


%% Variables needed for all oTDR Runs

% DK -- this appears to be an effort to linearize orthFlg, which is a
% cell array of potentially different sized vectors. Since vectors
% of different sizes cannot be combined into a matrix, we perform this
% linearization below.
% orthSignals = reshape(vertcat(orthFlg{:})', [],1)>0;
orthSignals = horzcat(orthFlg{:}) > 0;

% store coded param index for each on-target and off-target sRA
foo = horzcat(codedParams4Dataset{:});
param4SRA = foo(orthSignals);
param4SRA_offTarget = foo(~orthSignals);

% store times used for each on-target and off-target sRA
t4SRA = t4SRAAll(orthSignals);
t4SRA_offTarget = t4SRAAll(~orthSignals);

clear foo t4SRAAll 

%% compute sRA (where select sRAs are orthogonalized)
if ~bit4Bootstrap
%     if bitSerialOTDR
%         
%         loocvFlg = true; % no regularization
%         
%         % use TDR with single predictor
%         sRA = NaN(N,size(codedParams,2));
%         for k = 1:size(codedParams,2)
%             [foo] = runTDR(permute(RBar{:},[3 2 1]), N, [codedParams(:,k), ones(C,1)], trCountmtx, loocvFlg);
%             % collect just the dRA corresponding to the coded param (not
%             % the constant term)
%             sRA(:,k) = squeeze(foo(:,:,1));
%             clear foo
%         end
%         
%         sRA_offTarget = [];
%         errPerIter = NaN;
%         optimGradFinal = 0;
%     else
        [sRAAll,~,errPerIter,optimGradFinal] = optimizeOTDR(RBar, CC, trCountmtx, spcConstraint, orthFlg);

        % extract on-target (orthogonalized) sRAs
        sRA = sRAAll(:, orthSignals);
        % [~, sRA] = normVects(sRA); already done inside of optimizeOTDR

        % extract off-target (non-orthogonalized) sRAs
        sRA_offTarget = sRAAll(:, ~orthSignals);

        clear sRAAll 
%     end
end

%% compute nonOrthogonal sRA (no sRAs are orthogonalized)
% do not perform for serial oTDR
if ~bitSerialOTDR
    orthFlgStar = [];
    for i = 1:length(CC);
        orthFlgStar{i} = false(1, size(CC{i},2));
    end
    % note that errPerIter and optimGradFinal will be overwritten by this
    % call. However, in no case should this call to optimizeOTDR
    % (bit4Bootstrap && ~bitSerialOTDR) and above call to optimizeOTDR
    % (~bit4Bootstrap && bitSerialOTDR) be made in the same call to
    % runOTDR.m
    [sRAStarAll,~,errPerIter,optimGradFinal] = optimizeOTDR(RBar, CC, trCountmtx, [], orthFlgStar);
    % normalize sRAs to unit vectors
    [~, sRAStarAll] = normVects(sRAStarAll); 
        
    % extract on-target sRAs
    sRAStar = sRAStarAll(:, orthSignals);
    
    % extract off-target sRAs
    sRAStar_offTarget = sRAStarAll(:, ~orthSignals);
    
    clear sRAStarAll
end
%% eliminate off-target sRAs -- DEFUNCT

% % because we may include off-target predictors to soak-up variance (e.g.,
% % benefit predictor in dataset 2 or 3), we have to remove the associated
% % sRAs, as these off-target predictors are ignored. Since as of now, we
% % hardcode which predictors are included in which dataset (see above), we
% % have to guess which columns of sRAs correspond to which predictors.
% if size(sRA,2) == 4
%     if size(CC{2},2) == 3 && size(CC{1},2) == 1
%         % second sRA corresponds to off-target benefit predictor and will
%         % be eliminated.
%         sRA(:,2) = [];
%         sRAStar(:,2) = [];
%     else
%         error('unknown sRA to predictor mapping')
%     end
% elseif size(sRA,2) == 5
%     if size(CC{2},2) == 3 && size(CC{1},2) == 2
%         % second sRA corresponds to off-target ER predictor and third sRA
%         % corresponds to off-target benefit predictor and will be
%         % eliminated.
%         sRA(:,2:3) = [];
%         sRAStar(:,2:3) = [];
%     else
%         error('unknown sRA to predictor mapping')
%     end
%     
% end

%%
% bigA_projsRA1 = bsxfun(@minus, bigA, mean(bigA))*sRA(:,1);
% bigA_projsRA2 = bsxfun(@minus, bigA, mean(bigA))*sRA(:,2);
% bigA_projsRA3 = bsxfun(@minus, bigA, mean(bigA))*sRA(:,3);
% bigA_projPC1 = bsxfun(@minus, bigA, mean(bigA))*PCs(:,1);
% bigA_projPC2 = bsxfun(@minus, bigA, mean(bigA))*PCs(:,2);
% bigA_projPC3 = bsxfun(@minus, bigA, mean(bigA))*PCs(:,3);
% 
% A_projsRA1 = reshape(bigA_projsRA1, T, C);
% A_projsRA2 = reshape(bigA_projsRA2, T, C);
% A_projsRA3 = reshape(bigA_projsRA3, T, C);
% A_projsPC1 = reshape(bigA_projPC1, T, C);
% A_projsPC2 = reshape(bigA_projPC2, T, C);
% A_projsPC3 = reshape(bigA_projPC3, T, C);


% if ~isempty(DataSingle)
%    bigASingle = preprocess4TDR(data2Tensor(DataSingle), preprocessingSummary);
%    bigASingle_projsRA1 = bsxfun(@minus, bigASingle, mean(bigA))*sRA(:,1);
%    Asingle_projsRA1 = reshape(bigASingle_projsRA1, T, 1);
%    % DK -- extend to produce projections for the singleton on all sRAs:
%    Asingle_projsRA2 = reshape(bsxfun(@minus, bigASingle, mean(bigA))*sRA(:,2), T, 1);
%    Asingle_projsRA3 = reshape(bsxfun(@minus, bigASingle, mean(bigA))*sRA(:,3), T, 1);  
% end
% 

%% variance analysis of sRAs

if ~bit4Bootstrap
    
    % on-target sRA variance -- separately including common condition
    % response (times x neurons). Also store matrix of random vectors.
    % These are not specific to any one variance analysis, and thus are
    % extracted once and serve for the entire dataset.
    [Summary.sRA.varAnalysis, RARand] = TDRvarAnalysis(sRA, dataTensor, codedParams, commCond, nullStatModel, RARand, gofThresh);

%     % do not store Random vectors when doing serial oTDR
    if ~bitSerialOTDR
        Summary.RARand = RARand;
        % reset RARand so that subsequent calls to TDRvarAnalysis do not
        % use it
        RARand = [];
    end
    
    % off-target sRA variance
    if ~isempty(sRA_offTarget)
        [Summary.sRA_offTarget.varAnalysis] = TDRvarAnalysis(sRA_offTarget, dataTensor, codedParams, [], nullStatModel, RARand, gofThresh);
    end
end

%% variance of sRA Star
% do not perform for serial oTDR or when bootstrapping
if ~bit4Bootstrap && ~bitSerialOTDR

    % ontarget sRAStar variance (variance of non-orthogonalized, on-target
    % sRAs)
    [Summary.sRAStar.varAnalysis] = TDRvarAnalysis(sRAStar, dataTensor, codedParams,[], nullStatModel, RARand, gofThresh);
end

%% variance of PCs
% do not perform for serial oTDR or when bootstrapping
if ~bit4Bootstrap && ~bitSerialOTDR
    % PC variance
    [Summary.PCs.varAnalysis] = TDRvarAnalysis(PCs(:,1:size(sRA,2)), dataTensor, codedParams, commCond, nullStatModel, RARand, gofThresh);
end

%% PC-space of the common condition response

% do not perform for serial oTDR or when bootstrapping
if ~bit4Bootstrap && ~bitSerialOTDR

    % find PC space asssociated with common condition response
    [Summary.PCs_CC.RA] = pca(commCond);
    
    % signal-of-interest variance explained by CC PCs
    [Summary.PCs_CC.varAnalysis] = TDRvarAnalysis(Summary.PCs_CC.RA(:,1:size(sRA,2)), dataTensor, codedParams, commCond, nullStatModel, RARand, gofThresh);

    % compute variance explained by PCs (for simplicity, only save the number
    % of PCs as sRAs)
    [Summary.PCs_CC.varAnalysis.CC] = commonCondVar(Summary.PCs_CC.RA(:,1:size(sRA,2)),commCond);    

    % compute the associated p-value (using null distribution from .sRA
    % substruct) For now, compute EMPIRICAL p-value, since the distribution of
    % variance does not match a cannonical distribition
    Summary.PCs_CC.varAnalysis.CC.pV_RA = NaN(size(Summary.PCs_CC.varAnalysis.CC.V_RA));
    for k = 1:length(Summary.PCs_CC.varAnalysis.CC.V_RA)
        Summary.PCs_CC.varAnalysis.CC.pV_RA(k) = sigTest(Summary.PCs_CC.varAnalysis.CC.V_RA(k), Summary.sRA.varAnalysis.CC.V_rand, 'upper', 'empirical');
    end

    % store projection of full data response onto PCs_CC (with mean
    % subtraction)
    Summary.PCs_CC.proj = projData(dataTensor,Summary.PCs_CC.RA(:,1:size(sRA,2)));

    % STORE projection of full data response onto PCs_CC (NO mean
    % subtraction)
    Summary.PCs_CC.projNoMeanSubtraction = projData(bsxfun(@plus,dataTensor,commCond4Proj),...
        Summary.PCs_CC.RA(:,1:size(sRA,2)),false);
    
    % store projection of CC response onto PCs_CC (no mean subtraction)
    Summary.PCs_CC.projCC = projData(commCond,Summary.PCs_CC.RA(:,1:size(sRA,2)),false);
end


%% STORE sRAs on-target
if ~bit4Bootstrap
    Summary.sRA.RA = sRA;
    Summary.sRA.param4RA = param4SRA;
    Summary.sRA.t4RA = t4SRA;
    % compute projections on the fly for on-target sRAS
    Summary.sRA.proj = projData(dataTensor,sRA);
    % projections of excluded conditions
    if ~isempty(dataTensorExcluded)
        Summary.sRA.projExcluded = projData(dataTensorExcluded,sRA);
    end

    %%%% STORE sRA on-target projections WITHOUT mean subtraction. Just add
    %%%% back common condition response
    Summary.sRA.projNoMeanSubtraction = projData(bsxfun(@plus,dataTensor,commCond4Proj),sRA,false);

    %%%% STORE projection of common condition response onto sRAs (no mean
    %%%% subtraction): 
    % do not perform for serial oTDR
    if ~bitSerialOTDR
        Summary.sRA.projCC = projData(commCond,sRA,false);
    end
    
    %%%% STORE sRA off-target
    if ~isempty(sRA_offTarget)
        Summary.sRA_offTarget.RA = sRA_offTarget;
        Summary.sRA_offTarget.param4RA = param4SRA_offTarget;
        % compute projections on the fly for off target sRAs
        Summary.sRA_offTarget.proj = projData(dataTensor,sRA_offTarget);
        Summary.sRA_offTarget.t4RA = t4SRA_offTarget;
        % projections of excluded conditions 
        if ~isempty(dataTensorExcluded)
            Summary.sRA_offTarget.projExcluded = projData(dataTensorExcluded,sRA_offTarget);
        end
    end
end

%% STORE sRA Star (non-orthogonalized, on-target)
% do not perform for serial oTDR
if ~bitSerialOTDR
    Summary.sRAStar.RA = sRAStar;
    % compute projections on the fly for on-target sRA star (non-orthogonalized)
    Summary.sRAStar.proj = projData(dataTensor,sRAStar);
    % redundant with .sRA field, but useful to have in same relative location:
    Summary.sRAStar.param4RA = param4SRA;
    Summary.sRAStar.t4RA = t4SRA;
    % projections of excluded conditions 
    if ~isempty(dataTensorExcluded)
        Summary.sRAStar.projExcluded = projData(dataTensorExcluded,sRAStar);
    end

    %%% COMPUTE the angle between sRAs for sRASTAR
    if ~bit4Bootstrap
        
        %%%%%%%%
        % Function for computing angle and p-values between sRAs.
        
        % Done first between the sRAs for the current epoch.
        Summary.sRAStar.angleAnalysis = ...
             angleFn(Summary.sRAStar.RA,[],Summary.RARand,bitPlotRandAngleFit);
         
        % Then between current epoch sRAs and alternative vectors provided
        % in vect4Angle:
        if ~isempty(vect4Angle)
            Summary.sRAStar.angleAnalysis_withVect4Angle = ...
                angleFn(Summary.sRAStar.RA,vect4Angle,Summary.RARand,bitPlotRandAngleFit);
        end

    end
    
    %%%% STORE sRA Star off-target
    if ~isempty(sRAStar_offTarget)
        Summary.sRAStar_offTarget.RA = sRAStar_offTarget;
        % compute projections on the fly
        Summary.sRAStar_offTarget.proj = projData(dataTensor,sRAStar_offTarget);
    end
end

%% STORE PCs
% do not perform for serial oTDR or when bootstrapping
if ~bit4Bootstrap && ~bitSerialOTDR

    % store PC vectors in high D space
    Summary.PCs.RA = PCs;
    % store PC eigen values
    Summary.PCs.eigenVal = eigenVal;
    % compute projectons on the fly (with mean subtraction)
    Summary.PCs.proj = projData(dataTensor,PCs(:,1:size(sRA,2)));
    % projections of excluded conditions 
    if ~isempty(dataTensorExcluded)
        Summary.PCs.projExcluded = projData(dataTensorExcluded,PCs(:,1:size(sRA,2)));
    end
   
    % STORE projection of full data response onto PCs (NO mean subtraction)
    Summary.PCs.projNoMeanSubtraction = projData(bsxfun(@plus,dataTensor,commCond4Proj),...
        Summary.PCs.RA(:,1:size(sRA,2)),false);
    
    % store projection of common condition response onto PCs (no mean subtraction)
    Summary.PCs.projCC = projData(commCond,Summary.PCs.RA(:,1:size(sRA,2)),false);
    
    % compute common condition variance explained by PCs (for simplicity, only save the number
    % of PCs as sRAs)
    [Summary.PCs.varAnalysis.CC] = commonCondVar(Summary.PCs.RA(:,1:size(sRA,2)),commCond);    

    % compute the associated p-value (using null distribution from .sRA
    % substruct) For now, compute EMPIRICAL p-value, since the distribution of
    % variance does not match a cannonical distribition
    Summary.PCs.varAnalysis.CC.pV_RA = NaN(size(Summary.PCs.varAnalysis.CC.V_RA));
    for k = 1:length(Summary.PCs.varAnalysis.CC.V_RA)
        Summary.PCs.varAnalysis.CC.pV_RA(k) = sigTest(Summary.PCs.varAnalysis.CC.V_RA(k), Summary.sRA.varAnalysis.CC.V_rand, 'upper', 'empirical');
    end
    
end

%% generate random subspaces and compute probability of alignment index

% do not perform for serial oTDR or when bootstrapping
if ~bit4Bootstrap && ~bitSerialOTDR
    nSRA = size(sRA,2);
    
    if bitSubspaceStats
        nSample = 1000;
        [~,foo] = sampleRandSubspaces(nSRA, cov(bigA), [], nSample);
        
        SS_rand = NaN([size(foo{1}{1}) nSample]);
        % extract random subspaces
        for i = 1:length(foo)
            SS_rand(:,:,i) = foo{i}{:};
        end
        
        % instantiate
        nSSPair = (nSample^2-nSample)/2; % num of subspace pairs
        SS_angle = NaN(nSSPair,1);
        SS_alignIx = NaN(nSSPair,1);
        
        % loop through all random subspaces, computing subspace angle and
        % alignment index for each pair
        pairN = 0;
        for i = 1:nSample-1
            for j = i+1:nSample
                pairN = pairN + 1;
                SS_angle(pairN) = subspace(SS_rand(:,:,i),SS_rand(:,:,j));
                SS_alignIx(pairN) = align_ix(SS_rand(:,:,i),SS_rand(:,:,j));
            end
        end
    end
    
    
    % compute angles and stats for pairs of subspaces
    for i = 1:3
        switch i
            case 1
                pairName = 'sRA vs. PCs';
                fnA = 'sRA';
                fnB = 'PCs';
            case 2
                pairName = 'sRA vs. CC PCs';
                fnA = 'sRA';
                fnB = 'PCs_CC';
            case 3
                pairName = 'PCs vs CC PCs';
                fnA = 'PCs';
                fnB = 'PCs_CC';
        end
        
        subspaceComp(i).pairName = pairName;
        subspaceComp(i).angle = subspace(Summary.(fnA).RA(:,1:nSRA),...
            Summary.(fnB).RA(:,1:nSRA));
        subspaceComp(i).alignIx = align_ix(Summary.(fnA).RA(:,1:nSRA),...
            Summary.(fnB).RA(:,1:nSRA));
        
        if bitSubspaceStats
            subspaceComp(i).angleP = sigTest(subspaceComp(i).angle,SS_angle,...
                'upper','empirical');
            subspaceComp(i).alignIxP = sigTest(subspaceComp(i).alignIx,SS_alignIx,...
                'lower','empirical');
        end
    end
    
    % save to struct
    Summary.subspaceComp = subspaceComp;
end

end % main function
%% SUBFUNCTION for computing angles and p-values between sRAs

function aa = angleFn(RA1,RA2,RARand,bitPlotRandAngleFit)

        % determine if computing angles between single set of vectors or
        % between two sets of vectors
        if isempty(RA2)
            % computing angles between single set of vectors
            bitCorrXX = true;
            % set RA2 to RA1
            RA2 = RA1;
        else
            % computing angles between sets of vectors
            bitCorrXX = false;
        end
        
        % correlations
        [aa.corr, aa.pCorr] = corr(RA1,RA2);
        % angles
        aa.angle = real(acos((RA1'*RA2))*180/pi); %% angle 0-180
        
        % if computing between angles of single set, then covert to upper
        % triangle (returning a vector of values). Otherwise return full
        % matrix of all pairwise angles
        if bitCorrXX
            goo = logical(triu(ones(size(aa.corr)),1));
            aa.corr = aa.corr(goo);
            aa.pCorr = aa.pCorr(goo);
            aa.angle = aa.angle(goo);
        end
        
        % also compute folded angle (0 < theta < 90).
        % when reflecting theta, the distribution tends to pile up at <90
        % deg. This is difficult to fit. And so instead we invert the
        % distribution as theta' = 90 - theta' with "offset" = 90 and
        % "scale" = -1. The resulting distribution is >0 with most values
        % around 0 and appears as an exponential distribution that is well
        % captured by a gamma fit.
        aa.angleFolded = 90 - (aa.angle - 2*max(aa.angle-90,0));
        clear goo
        
        %%% compute null distribution of angles
        % pairwise angles between random vectors
        foo = real(acos((RARand'*RARand))*180/pi);
        % take upper triangle
        foo = foo(logical(triu(ones(size(foo)),1)));
        % tranform to folded values. See above note about additional offset and
        % scale transformation.
        fooF = 90 - (foo - 2*max(foo-90,0));
        
        % compute p-value for each signal against null hypothesis that angles are
        % SMALLER (i.e., further from 90) or (tested separately) LARGER (i.e.,
        % closer to 90) than expected by chance. Do this twice. Once assuming an
        % unfolded distribution (0 < theta < 180) and once assuming a folded
        % distribution (0 < theta < 90). Note that testing if the angle is larger
        % than expected only applies to the folded distribution, since we cannot
        % test if the angle is closer to 90deg if the null distribution is
        % symmetric around (or nearly around) 90deg.
        aa.pAngle = NaN(size(aa.angle));
        aa.pAngleFolded = NaN(size(aa.angle));
        aa.pAngleFolded_large = NaN(size(aa.angle));
        aa.nullH_model = 'gauss';
        aa.nullH_param = [];
        aa.nullHFolded_model = 'halfgaus';
        aa.nullHFolded_param = [];
        clear temp
        for i = 1:size(aa.pAngle,1)
            for j = 1:size(aa.pAngle,2)
                if aa.angle(i,j) >= 90
                    goo = 'upper';
                else
                    goo = 'lower';
                end
            
                % unfolded p-value
                [aa.pAngle(i,j),aa.nullH_param] = sigTest(aa.angle(i,j),foo,goo,'gauss','prms',aa.nullH_param);

                % folded p-value
                [aa.pAngleFolded(i,j),aa.nullHFolded_param] = sigTest(aa.angleFolded(i,j),fooF,'upper','halfgauss','prms',aa.nullHFolded_param);

                % folded p-value
                [aa.pAngleFolded_large(i,j)] = sigTest(aa.angleFolded(i,j),fooF,'lower','halfgauss','prms',aa.nullHFolded_param);
            end            
        end
        
        % Need to convert pAngle to 2-tailed probability given we did not have an
        % apriori assumption whether theta >> 90 or theta << 90. Assume that null
        % distribution is symmetric about 90, and thus the 2-tailed probability = 2
        % * 1-tailed probability.
        aa.pAngle = 2 * aa.pAngle;
        % un-transform folded angle back to 0 deg offset
        aa.angleFolded = 90 - aa.angleFolded;
        
        % plot histogram and model fit of random vector angles if requested for
        % debugging purposes
        if bitPlotRandAngleFit
            % compute histogram -- only necessary for 1 angle since they all share the
            % same underlying distribution
            [temp.freq,temp.bin] = hist(foo,30);
            [temp.freqF,temp.binF] = hist(fooF,30);
            
            figure('Name','Random vector angles and fits');
            
            % unfolded
            subplot(1,2,1)
            plot(temp.bin,temp.freq / sum(temp.freq) / diff(temp.bin(1:2)),'Color',[0.5 0.5 0.5])
            hold on
            plot(temp.bin,pdf('normal',temp.bin,aa.nullH_param(1),aa.nullH_param(2)),'Color','k');
            
            % folded
            subplot(1,2,2)
            plot(temp.binF,temp.freqF / sum(temp.freqF) / diff(temp.binF(1:2)),'Color',[0.5 0.5 0.5])
            hold on
            plot(temp.binF,2*pdf('normal',temp.binF,aa.nullHFolded_param(1),aa.nullHFolded_param(2)),'Color','k');
        end


end % end subfunction


