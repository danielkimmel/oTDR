function [Summary,RARand] = TDRvarAnalysis(RAs, dataTensor, codedParams, varargin)
% TODO: update docstring
%
% [Summary] = TDRvarAnalysis(RAs, dataTensor, codedParams, [CC],[nullStatModel],[RARand],[gofThresh])
%
% [CC] -- (optional) T x M matrix of common condition response
%       across M neurons and T time points.
%
% , varRand_t, RSVRand1_t, RSVRand2_t, RSVRand3_t, ISVRand1_t, ISVRand2_t, ISVRand3_t
%
% optional inputs
% ...
%
% OUTPUT:
% Summary -- struct containing numerous fields related to variance
% explained by the RAs
%
% RARand -- M x S matrix of S random vectors biased to the subspace of the
% M-dimensional data 

if length(varargin) > 0
    CC = varargin{1}; % common condition response 
else
    CC = [];
end
if length(varargin) > 1
    nullStatModel = varargin{2}; % string for model to fit to null distribution. See sigTest.m for options.
end
if length(varargin) <= 1 || isempty(nullStatModel)
    nullStatModel = 'gamma';
end

% matrix of random vectors (1 per column) to be passed to varianceIXs.m
% The main advantage of passing the vectors rather than computing de novo
% is both a) time, and b) when the dataspace is not represented by
% dataTensor. In this latter case, we MUST provide the random vectors. This
% occurs, for example, when computing serial sRAs in which the call to
% runOTDR.m does not have the full dataTensor. Strangely, the code is
% faster when generating random vectors on the fly. This has to do with the
% call to varianceIXs.m, which appears to run faster when computing
% variance for a single vector one at a time rather than a matrix of
% vectors.
if length(varargin) > 2
    RARand = varargin{3}; 
else
    RARand = [];
end

% scalar specifying threshold for goodness-of-fit in fitting parametric
% distribution to empirical null model for purposes of computing p-values.
% When threshold is exceeded, error is thrown.
if length(varargin) > 3
    gofThresh = varargin{4}; 
else
    gofThresh = 0.17;
end

%%
nRA = size(RAs, 2); % numher of regression axes
nPred = size(codedParams,2); % number of predictors
[T, N, C] = size(dataTensor);
varIxs = varianceIXs(eye(N), dataTensor, codedParams, false, CC);
varNeu_t = varIxs.var_t; % variance across conditions, per time
totalVar_t = sum(varNeu_t,2); 
Summary.varNeu_t = varNeu_t;
Summary.totalVar_t = totalVar_t;
if ~isempty(CC)
    % variance across time is only computed when common condition response
    % CC is provided
    varNeu_cond = varIxs.var_cond; % variance across time, per condition
    totalVar_cond = sum(varNeu_cond,2);
    Summary.varNeu_cond = varNeu_cond;
    Summary.totalVar_cond = totalVar_cond;
end

bigA = reshape(permute(dataTensor,[1 3 2]), [], N);
numSamples = 10000;
% do not need to generate random vectors if already provided
if isempty(RARand)
    [varIxs_rand,RA_rand] = sampleRandSubspaces(1, cov(bigA), 'varianceIXs', numSamples, dataTensor, codedParams, true, CC);
    
    % initialize
    varRand_t = NaN(T,numSamples);
    RSVRand_t = NaN(T,numSamples,nPred);
    ISVRand_t = NaN(size(RSVRand_t));
    RSVRand_partial_t = NaN(T,numSamples,nPred,nPred-1);
    ISVRand_partial_t = NaN(size(RSVRand_partial_t));
    varRand_cond = NaN(C,numSamples); % for variance over time
    RARand = NaN(N,numSamples); % vectors
    % difference in sqaured correlation coefficients r2_bp - r2_ab,
    % where r2_bp is correlation between off-target signal b (dimension
    % 4) and r2_ab is correlation between off-target signal b and
    % on-target signal a (dimension 3).
    diffCorr_proj_signal_Rand = NaN(T,numSamples,nPred,nPred-1);
    
    for s = 1:numSamples
        varRand_t(:, s) = varIxs_rand{s}.var_t;
        % NOTE: The disjoint relationship between RSV and ARSV (r > 0 vs. r
        % < 0) inflates the frequency of RSV = 0 (when r < 0) and ARSV = 0
        % (when r > 0), which unnaturally distorts the null distributions
        % of these values. We do not believe there is a fundamental
        % neurobiological difference between r > 0 and r < 0. Thus, small
        % fluctuations in r near r = 0 are responsible for an absolute
        % classification of RSV vs. ARSV. The null distrubition should not
        % be sensitive to these small fluctuations and should instead
        % reflect positive and negative relationships equally. We therefore
        % define the null distribution of RSV as RSV_hat = VE_hat * r^2,
        % regardless of the sign of r. This can be achieved post-hoc by
        % adding RSV and ARSV, since they are disjoint. We'll check the
        % disjoint nature just in case.
%         if ~all(all(varIxs_rand{s}.RSV_t==0 | varIxs_rand{s}.aRSV_t==0))
%             error('Values of RSV and aRSV should be disjoint. At no point should both values be non-zero')
%         end
        RSVRand_t(:, s, :) = varIxs_rand{s}.RSV_t;
%         RSVRand_t(:, s, :) = varIxs_rand{s}.RSV_t + varIxs_rand{s}.aRSV_t;
        ISVRand_t(:, s, :) = varIxs_rand{s}.ISV_t;
        RSVRand_partial_t(:, s, :, :) = varIxs_rand{s}.RSV_t_partial;
%         RSVRand_partial_t(:, s, :, :) = varIxs_rand{s}.RSV_t_partial + varIxs_rand{s}.aRSV_t_partial;
        ISVRand_partial_t(:, s, :, :) = varIxs_rand{s}.ISV_t_partial;
        RARand(:,s) = RA_rand{s}{:};
        diffCorr_proj_signal_Rand(:, s, :, :) = varIxs_rand{s}.diffCorr_proj_signal; 
        if ~isempty(CC)
            % variance across time is only computed when common condition response
            % CC is provided
            varRand_cond(:, s) = varIxs_rand{s}.var_cond;
        end
    end
    
else
    % check that sufficient samples are provided and in right
    % dimensionality
    if size(RARand,2) < numSamples
        error('Provided random vectors have fewer samples than requested')
    elseif size(RARand,2) > numSamples
        % if more samples are provided than needed, pair down
        RARand(:,numSamples+1:end) = [];
    end
    if size(RARand,1) ~= N
        error('Dimensionality of random vectors does not match dimensionality of data')
    end
    
    % compute variance explained by random samples
    varIxs_rand = varianceIXs(RARand, dataTensor, codedParams, true, CC);
    
    % collect
    varRand_t = varIxs_rand.var_t;
    RSVRand_t = varIxs_rand.RSV_t;
    ISVRand_t = varIxs_rand.ISV_t;
    RSVRand_partial_t = varIxs_rand.RSV_t_partial;
    ISVRand_partial_t = varIxs_rand.ISV_t_partial;
    diffCorr_proj_signal_Rand = varIxs_rand.diffCorr_proj_signal;
    if ~isempty(CC)
        % variance across time is only computed when common condition response
        % CC is provided
        varRand_cond = varIxs_rand.var_cond;
    end
end

Summary.varRand_t = varRand_t;
Summary.RSVRand_t = RSVRand_t;
Summary.ISVRand_t = ISVRand_t;
Summary.RSVRand_t = RSVRand_partial_t;
Summary.ISVRand_t = ISVRand_partial_t;
Summary.diffCorr_proj_signal_Rand = diffCorr_proj_signal_Rand;
if ~isempty(CC)
    Summary.varRand_cond = varRand_cond;
end

%% variance indexes for RA for RSV and ISV dimensionality is (times x axes x signals)

% initialize
V_RA = NaN(T,nRA);
VE_RA = NaN(size(V_RA));
pV_RA = NaN(size(V_RA));
RSV_RA = NaN(T,nRA,nPred);
RSVE_RA = NaN(size(RSV_RA));
pRSV_RA = NaN(size(RSV_RA));
% ARSV_RA = NaN(size(RSV_RA));
% ARSVE_RA = NaN(size(RSV_RA));
% pARSV_RA = NaN(size(RSV_RA)); 
ISV_RA = NaN(size(RSV_RA));
ISVE_RA = NaN(size(RSV_RA));
pISV_RA = NaN(size(RSV_RA));
corr_proj_signal = NaN(size(RSV_RA));

% initialize partial correlations
RSV_RA_partial = NaN(T,nRA,nPred,nPred-1);
RSVE_RA_partial = NaN(size(RSV_RA_partial));
pRSV_RA_partial = NaN(size(RSV_RA_partial));
% ARSV_RA_partial = NaN(size(RSV_RA_partial));
% ARSVE_RA_partial = NaN(size(RSV_RA_partial));
% pARSV_RA_partial = NaN(size(RSV_RA_partial));
ISV_RA_partial = NaN(size(RSV_RA_partial));
ISVE_RA_partial = NaN(size(RSV_RA_partial));
pISV_RA_partial = NaN(size(RSV_RA_partial));

% initialize difference in sqaured correlation coefficients r2_bp - r2_ab,
% where r2_bp is correlation between off-target signal b (dimension 4) and
% r2_ab is correlation between off-target signal b and on-target signal a
% (dimension 3).
diffCorr_proj_signal = NaN(size(RSV_RA_partial));
pDiffCorr_proj_signal = NaN(size(RSV_RA_partial));

% initialize for variance across time
V_RA_cond = NaN(C,nRA);
VE_RA_cond = NaN(C,nRA);
pV_RA_cond = NaN(size(V_RA_cond));

for k = 1:nRA
    varIxs = varianceIXs(RAs(:, k), dataTensor, codedParams,1,CC); % argument "1" (logical) to perform partial correlations in addition to full correlations
    V_RA(:,k) = varIxs.var_t;
    VE_RA(:,k) = bsxfun(@times, V_RA(:,k), 100./totalVar_t);
    RSV_RA(:,k,:) = varIxs.RSV_t(:, :, :);
    RSVE_RA(:,k,:) = bsxfun(@times, RSV_RA(:,k,:), 100./totalVar_t);
%     ARSV_RA(:,k,:) = varIxs.aRSV_t(:, :, :);
%     ARSVE_RA(:,k,:) = bsxfun(@times, ARSV_RA(:,k,:), 100./totalVar_t);
    ISV_RA(:,k,:) = varIxs.ISV_t(:, :, :);
    ISVE_RA(:,k,:) = bsxfun(@times, ISV_RA(:,k,:), 100./totalVar_t);
    pV_RA(:,k) = sigTest(V_RA(:,k), varRand_t, 'upper', nullStatModel,'errThresh',gofThresh);
    for kk = 1:size(RSV_RA,3)
        pRSV_RA(:,k,kk) = sigTest(RSV_RA(:, k, kk), RSVRand_t(:,:,kk), 'upper', nullStatModel,'errThresh',gofThresh);
%         pARSV_RA(:,k,kk) = sigTest(ARSV_RA(:, k, kk), RSVRand_t(:,:,kk), 'upper', nullStatModel);
        pISV_RA(:,k,kk) = sigTest(ISV_RA(:, k, kk), ISVRand_t(:,:,kk), 'upper', nullStatModel,'errThresh',gofThresh);
    end
    corr_proj_signal(:,k,:) = varIxs.Cr_t(:, :, :);

    if ~isempty(CC)
        % variance across time is only computed when common condition response
        % CC is provided
        V_RA_cond(:,k) = varIxs.var_cond;
        VE_RA_cond(:,k) = bsxfun(@times, V_RA_cond(:,k), 100./totalVar_cond);
        % gamma distribution is less of a good fit for cross-time variance,
        % thus requiring higher error threshold
        pV_RA_cond(:,k) = sigTest(V_RA_cond(:,k), varRand_cond, 'upper', nullStatModel,...
            'errThresh',gofThresh);
    end
    
    % collect partial RSV, ISV
    RSV_RA_partial(:,k,:,:) = varIxs.RSV_t_partial(:, :, :, :);
    RSVE_RA_partial(:,k,:,:) = bsxfun(@times, RSV_RA_partial(:,k,:,:), 100./totalVar_t);
%     ARSV_RA_partial(:,k,:,:) = varIxs.aRSV_t_partial(:, :, :, :);
%     ARSVE_RA_partial(:,k,:,:) = bsxfun(@times, ARSV_RA_partial(:,k,:,:), 100./totalVar_t);
    ISV_RA_partial(:,k,:,:) = varIxs.ISV_t_partial(:, :, :, :);
    ISVE_RA_partial(:,k,:,:) = bsxfun(@times, ISV_RA_partial(:,k,:,:), 100./totalVar_t);
    diffCorr_proj_signal(:,k,:,:) = varIxs.diffCorr_proj_signal(:, :, :, :);
    % compute p-values
    for kk = 1:size(RSV_RA_partial,3)
        % find signal position
%         pos = varIxs.signalN(kk,:);
        for j = 1:size(RSV_RA_partial,4)
            pRSV_RA_partial(:,k,kk,j) = sigTest(RSV_RA_partial(:, k, kk, j), RSVRand_partial_t(:,:,kk,j), 'upper', nullStatModel,'errThresh',gofThresh);
%             pARSV_RA_partial(:,k,kk,j) = sigTest(ARSV_RA_partial(:, k, kk, j), RSVRand_partial_t(:,:,kk,j), 'upper', nullStatModel);
            pISV_RA_partial(:,k,kk,j) = sigTest(ISV_RA_partial(:, k, kk, j), ISVRand_partial_t(:,:,kk,j), 'upper', nullStatModel,'errThresh',gofThresh);
%             pRSV_RA_partial(:,k,kk,j) = sigTest(RSV_RA_partial(:, k, kk, j), RSVRand_t(:,:,pos(j)), 'upper', nullStatModel);
%             pARSV_RA_partial(:,k,kk,j) = sigTest(ARSV_RA_partial(:, k, kk, j), RSVRand_t(:,:,pos(j)), 'upper', nullStatModel);
%             pISV_RA_partial(:,k,kk,j) = sigTest(ISV_RA_partial(:, k, kk, j), ISVRand_t(:,:,pos(j)), 'upper', nullStatModel);
            pDiffCorr_proj_signal(:,k,kk,j) = sigTest(diffCorr_proj_signal(:, k, kk, j), diffCorr_proj_signal_Rand(:,:,kk,j), 'upper', 'empirical');
        end
    end

end

%% compute variance related to common condition response

if ~isempty(CC)
    
    % compute common condition variance:    
    [CC_sum] = commonCondVar(RAs,CC);    
    
    % generate random vectors in the space of the common condition. Note
    % that here we do not have the option to use pre-existing random
    % vectors and instead have to generate them de novo.
    [CC_rand] = sampleRandSubspaces(1, cov(CC), 'commonCondVar', numSamples, CC);

    % extract random variance
    V_rand = NaN(1,numSamples);
    parfor s = 1:numSamples
        V_rand(s) = CC_rand{s}.V_RA;
    end
    
    % compute p-value for variance explained by each axis. For now, compute
    % EMPIRICAL p-value, since the distribution of variance does not match
    % a cannonical distribition
    CC_sum.pV_RA = NaN(size(CC_sum.V_RA));
    for k = 1:nRA
        CC_sum.pV_RA(k) = sigTest(CC_sum.V_RA(k), V_rand, 'upper', 'empirical');
    end
    
    % Save to Summary
    Summary.CC = CC_sum;
    Summary.CC.V_rand = V_rand;
end


%% save to summary
Summary.V_RA = V_RA;
Summary.VE_RA = VE_RA;
Summary.RSV_RA = RSV_RA;
Summary.RSVE_RA = RSVE_RA;
% Summary.ARSV_RA = ARSV_RA;
% Summary.ARSVE_RA = ARSVE_RA;
Summary.ISV_RA = ISV_RA;
Summary.ISVE_RA = ISVE_RA;
Summary.pV_RA = pV_RA;
Summary.pRSV_RA = pRSV_RA;
% Summary.pARSV_RA = pARSV_RA;
Summary.pISV_RA = pISV_RA;
Summary.corr_proj_signal = corr_proj_signal;

% partial correlations
Summary.signalN4Partial = varIxs.signalN;
Summary.RSV_RA_partial = RSV_RA_partial;
Summary.RSVE_RA_partial = RSVE_RA_partial;
Summary.pRSV_RA_partial = pRSV_RA_partial;
% Summary.ARSV_RA_partial = ARSV_RA_partial;
% Summary.ARSVE_RA_partial = ARSVE_RA_partial;
% Summary.pARSV_RA_partial = pARSV_RA_partial;
Summary.ISV_RA_partial = ISV_RA_partial;
Summary.ISVE_RA_partial = ISVE_RA_partial;
Summary.pISV_RA_partial = pISV_RA_partial;

% difference in sqaured correlation coefficients 
Summary.diffCorr_proj_signal = diffCorr_proj_signal;
Summary.pDiffCorr_proj_signal = pDiffCorr_proj_signal;

% temporal variance
if ~isempty(CC)
    % variance across time is only computed when common condition response
    % CC is provided
    Summary.V_RA_cond = V_RA_cond;
    Summary.VE_RA_cond = VE_RA_cond;
    Summary.pV_RA_cond = pV_RA_cond;
end


end

