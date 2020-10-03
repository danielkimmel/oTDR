function [BSSummary,oTDRSummaryAll] = bootstrapsRAs(N, cond2include, ...
    codedParams, dataset_times, codedParams4Dataset, spcConstraint, ...
    varargin)
% [BSSummary,oTDRSummaryAll] = bootstrapsRAs(N, cond2include, 
%   codedParams, dataset_times, codedParams4Dataset, spcConstraint, ...
%   [name,value])
%
% Args:
%   N: 1 x n struct where each element refers to individual neuron with fields:
%       .psth: matrix (t * r) of firing rate of neuron n on trial r
%           (columns) at time t (rows) relative to event-of-interest. 
%       .variableValues4Trial: matrix (r * p) of the values of the p
%           task variables (columns) for trial r (rows).
%       (If field ".cond" is present, it will be removed. Additional fields
%       will be ignored) 
%
%   cond2include: matrix (c * p) of the values of the p task variables
%       (columns) for the c condition (columns) for which you wish to
%       compute the trial-average. Values and number of columns p must
%       match those provided in N.conditionValues4Trial (see above). Note:
%       matrix likely already defined when making original call to oTDR().
%       Can also be found in oTDR output struct: oTDRSummary.Predictors
%
% THE FOLLOWING ARGS ARE REQUIRED TO PASS TO oTDR(). See oTDR.m for
% details:
%   codedParams: likely already defined when making original call to
%       oTDR(). Can also be found in oTDR output struct:
%       oTDRSummary.codedParams. Note, however, that this version is scaled
%       to range [0 1], whereas original variables may take on different
%       values. The scaling is performed automatically in oTDR() regardless
%       and thus should not affect output.
%   dataset_times: if not originally passed to oTDR(), can be found
%       in oTDR output struct: oTDRSummary.sRA.t4RA
%   codedParams4Dataset: likely already defined when making original call
%       to oTDR(). Can also be found in oTDR output struct:
%       oTDRSummary.meta.codedParams4Dataset
%   spcConstraint: this may have been defined when making original call to
%       oTDR(), but instead may have been computed within oTDR if oTDR was
%       called with numPCs4OTDR specified. In this case, the space
%       constraint can be found in oTDR output struct:
%       oTDRSummary.meta.spcConstraint
%       If no constraint was applied, pass an empty matrix [].
%       NOTE: we do not support constraining to top D PCs at the level of
%       bootstrapping because it would cause each resampled dataset to have
%       a different space constraint.
%
% THE FOLLOWING ARE OPTIONAL ARGS, SPECIFIED AS NAME, VALUE PAIRS:
%   param4nrm: struct with the following fields:
%       .mu: 1 x n vector of the mean of each neuron, which will be
%           subtracted from the resampled mean in each bootstrapped dataset.
%       .nrm: 1 x n vector of the std of each neuron, by which each
%           resampled mean will be normalized (i.e., divided by) after
%           subtraction by .mu (see above). 
%       NOTE: These normalizing terms are passed to each call to oTDR() and
%           used to z-score normalize the responses from each resampled
%           dataset. This ensures that all datasets share the same
%           normalization relative to the original data. The values for .mu
%           and .nrm can be conveniently found in the output struct from
%           oTDR() as computed for the full (i.e. non-resampled dataset) --
%           see oTDRSummary.meta.preprocessingSummary.
%       If not provided, each resampled dataset will be z-score normalized
%           per the mean and std of the given dataset. Thus each dataset
%           will have a different normalization. This is not advised.
%   nBootstraps: int specifying number of resampled bootstrapped datasets
%       to compute. Default = 700.
%   bitStats: logical on whether to compute summary stats and hypothesis
%       testing on the correlations within and between bootstrapped
%       datasets. (Note that the correlations are always computed, but the
%       summary stats and testing is optional.) Results are returned in
%       .stats substruct. Including bitStats 
%   alpha: scalar (range [0 1]) specifying the lower (*100) confidence
%       interval returned in the .stats substruct (see below).
%
% Returns:
%   BSSummary: struct with following fields:
%       .sRA: S x B x N tensor of the sRA coefficients for each of the S
%           task variables, B resampled datasets, and N neurons.
%       .optimGradFinal: 1 x B vector of the final gradient of the
%           oTDR optimization routine for each of B resampled datasets.
%           Unusually large values would suggest that oTDR for a given
%           dataset failed to converge. 
%       .sRA_p: N x S matrix of p-values for each of the N neurons and S
%           task variables. For a given neuron and task varible, the
%           p-value reflects the probability that a coefficient of zero
%           came from the distribution of sRA coefficients across resampled
%           datasets. Assumes that this distribution is normal. This is
%           distinct from a t-test of whether the distribution of samples
%           had mean=0, as it is instead measuring the probability that,
%           after many experiments, the oTDR  estimate of the coefficient
%           for a given neuron would be zero.
%       .corrWithin: S x Z matrix of the Pearson correlation between sRAs
%           for the same task variable (given S variables) for each of Z
%           pairs of resampled datasets (given B resampled datasets, Z =
%           (B.^2-B)/2).
%       .corrBetween: Y x B matrix of the Pearson correlation between each
%           of Y pairs of sRAs of different task variables (given S
%           variables, Y = S-choose-2) within the same dataset (given B
%           resampled datasets).
%       .corrBetweenInd: Y x 2 matrix where each row reports the indecies of the pair of task
%           variables whose correlation is reported in the corresponding
%           row of .corrBetween (see above). The indecies reference the
%           column numbers in codedParams (input arg, also given in output
%           substruct, .meta -- see below). For instance, given S=3 task
%           variables, .corrBetweenInd = 
%               [1     2
%                1     3
%                2     3]
%           indicating, for instance, that the first row of .corrBetween
%           reports the correlation between task variables 1 and 2, as
%           given by the first and second columns of codedParams matrix.           
%       .meta: for convenience, substruct of oTDR settings with following 
%           fields (see above for definitions): cond2include, codedParams,
%           dataset_times, codedParams4Dataset, spcConstraint, param4nrm
%       .stats: substruct of summary statistics and hypothesis testing with
%           fields:
%           .within: struct summarizing within-variable (i.e.,
%               between resampled dataset) correlations, with fields:
%               .r: 1 x S vector of the mean correlation coefficient (from
%                   BSSummary.corrWithin) between datasets for each of the
%                   S sRAs.
%               .r_sigma: 1 x S vector of the standard deviation of the
%                   correlation coefficient (from BSSummary.corrWithin)
%                   between datasets for each of the S sRAs.
%               .r_ci: S x 2 matrix of the lower (row 1) and upper (row 2)
%                   confidence intervals (as given by .ci) of the
%                   correlation coefficient between datasets for each of
%                   the S sRAs. CIs are computed empirically from the
%                   distribution BSSummary.corrWithin (see above). 
%               .ci: 1 x 2 vector of the percentile used for the lower and
%                   upper confidence intervals, respectively, given by the
%                   rows of .r_ci.
%               .p: 1 x S vector of p-values from the right-tailed t-test
%                   on the probability of the mean of the distribution of
%                   between-dataset correlation coefficients (from
%                   BSSummary.corrWithin) differing from zero.
%           .between: struct summarizing between-variable (i.e., within
%               resampled dataset) correlations. Fields are thes same as
%               .stats.within (see above), but refer to the distribution of
%               correlation coefficients measured within each resampled
%               dataset between pairs of variables (from
%               BSSummary.corrBetween). Given S variables, there are Y
%               (S-choose-2) pairs, the indecies of each are given in the Y
%               x 2 matrix, BSSummary.corrBetweenInd. In addition, are the
%               following fields: 
%               .r_allTrials: 1 x Y vector of the correlation coefficient
%                   between each of Y pairs of variables based on sRAs
%                   computed from all trials (i.e., full dataset), unlike
%                   .r, which is the mean across resampled datasets.
%               .p_allTrials: 1 x Y vector of the p-value (via Fisher
%                   z-transform) associated with the correlation
%                   coefficient (given by .r_allTrials) between each of Y
%                   pairs of variables based on sRAs computed from all
%                   trials (i.e., full dataset), unlike .p, which is the
%                   p-value associated with mean across resampled datasets
%                   (via t-test).
%           .nullModel: struct with data for null model of with fields:
%               .rPred_dist: Y x Z matrix of correlation coefficients
%                   between Y pairs of sRAs of different task variables A
%                   and B (given S variables, Y = S-choose-2) across Z
%                   pairs of resampled datasets (given B resampled
%                   datasets, Z = (B.^2-B)/2), such that the null
%                   distribution r_AB = sqrt(r_AA)*sqrt(r_BB), given that
%                   the null assumption that the true correlation r_AB = 1
%                   but was corrupted by independent noise in estimating A
%                   and B. See Spearman, C. The Proof and Measurement of
%                   Association between Two Things. Am. J. Psychol. 15,
%                   72-101 (1904).
%               .rPred_dist_mean, .rPred_dist_median, .rPred_dist_CI: mean,
%                   median, and confidence intervals (see .stats.within.ci)
%                   of distributions in .rPred_dist
%               .pVal_alt_ttest, pVal_alt_signrank: 1 x Y vector of p-value
%                   (via t-test or signrank test, respectively) that the
%                   distribution of null correlation coefficients
%                   (i.e., .rPred_dist) had a mean or median, respectively,
%                   equal to the observed correlation coefficient computed
%                   across all trials (i.e., .stat.between.r_allTrials) vs.
%                   the alternative hypothesis that the distribution's
%                   mean/median was less than the observed value (i.e.,
%                   left-tailed). (Note that "alt" refers to the fact that
%                   at one point this was an alternative approach, but it
%                   is now the primary approach, as used in Kimmel et al.,
%                   Nat Comm, 2020.)
%
%   oTDRSummaryAll: Optional. If second output arg is requested in call to
%       bootstrapsRAs(), then the complete oTDR output struct will be saved
%       as an element of oTDRSummaryAll for each resampled dataset. This is
%       useful for debugging purposes, but is memory intensive.

%% set default values for optional args

% number of bootstrap datasets
nBootstraps = 700;

% structure with normalizing terms
param4nrm = [];

% logical on whether to compute stats sub-struct
bitStats = true;

% alpha level
alpha = 0.05;

%% collect optionally provided parameters
warnopts(assignopts(who, varargin));

%% process inputs

% change name for convenience
B = nBootstraps;
clear nBootstraps

% warn
if isempty(param4nrm)
    warning('No normalizing terms provided in optional input "param4nrm". This will cause each resampled dataset to be normalized internally, and is not recommended.')
end

%% prepare data

% if trial-average data already present, remove it
if isfield(N,'cond')
    N = rmfield(N,'cond');
end

        
%% run bootstraps


% if output arg specified by user, will store entire oTDR summary for
% each resampled dataset.
if nargout > 1
    oTDRSummaryAll = cell(1,B);
else
    oTDRSummaryAll = [];
end

fprintf('COMPUTING oTDR FOR %d RESAMPLED DATASETS...\n',B)

% loop through remaining datasets
for s = 1:B
    
    fprintf('Bootstrapped dataset %d...\n',s)
    
    % run oTDR
    otdr_temp = runBS();

    % for first bootstrapped dataset, get various dimensions and initialize
    % vars
    if s==1
        % get dimensions
        D = size(otdr_temp.sRAStar.RA,1); % num of neurons
        S = size(otdr_temp.codedParams,2);  % num of task variables

        % initialize: 
        BSSummary.sRA = NaN(S,B,D); % task variables x datasets x neurons
        BSSummary.optimGradFinal = NaN(1,B); 
    end        
    
    % store sRAs and final gradient 
    BSSummary.sRA(:,s,:) = otdr_temp.sRAStar.RA';
    BSSummary.optimGradFinal(s) = otdr_temp.meta.optimGradFinal;
    
    % optionally store entire summary
    if nargout > 1
        oTDRSummaryAll{s} = otdr_temp;
    end
    
    clear otdr_temp
    
end

%% compute p-values 
% For each neuron and sRA, computes the probability that a coefficient of
% zero came from the distribution of sRA coefficients across samples (i.e.,
% datasets). Assumes that this distribution is normal. This is distinct
% from a t-test of whether the distribution of samples had mean=0, as it is
% instead measuring the probability that, after many experiments, our 
% estimate of the coefficient for a given neuron would be zero. 

% TODO(G): I do wonder if the above statistical test is the most
% reasonable. Why not just do a t-test?

% initialize
BSSummary.sRA_p = NaN(D,S);

% loop through sRAs
for s = 1:S
        
    % compute p-value using Gaussian fit (tests goodness of fit and returns
    % error if gaussian is a poor model). 
    BSSummary.sRA_p(:,s) = sigTest(zeros(1,D),squeeze(BSSummary.sRA(s,:,:)),'both','normal',...
        'bitIgnorePoorFit',true)';
    
    % warn if some were excluded
   if sum(isnan(BSSummary.sRA_p(:,s)))/D > 0.02
       warning('sRA P-values for %d neuron(s) were not computed because of error in fitting null distribution for sRA %d',...
           sum(isnan(BSSummary.sRA_p(:,s))),s);
   elseif any(isnan(BSSummary.sRA_p(:,s)))
       warning('sRA P-values for %d neuron(s)s were not computed because of error in fitting null distribution for sRA %d',...
           sum(isnan(BSSummary.sRA_p(:,s))),s);
   end
   
    % TODO: See above. Consider this ttest as an altnerative or remove. 
    % compute p-value based on t-test. Here the test stat is
    % mu/(sigma/sqrt(n)) -- p-value are much smaller
%     [~,BSSummary.sRA_p_ttest(:,s)] = ttest(squeeze(BSSummary.sRA(s,:,:)));
    
%     for d = 1:D
%         % compute p-value based on signrank test
%         [~,BSSummary.sRA_p_signrank(d,s)] = signrank(squeeze(BSSummary.sRA(s,:,d)));
%     end
end

%% compute correlations between sRA samples

% get dimensions
W = (B.^2-B)/2; % num of within signal pairwise correlations
BT = nchoosek(S,2); % num of between signal pairs

% initialize
BSSummary.corrWithin = NaN(S,W); % within signal correlation
BSSummary.corrBetween = NaN(BT,B); % between signal correlation
BSSummary.corrBetweenInd = NaN(BT,2); % ind of signal pairs

% only perform if more than one boostrapped dataset
if B > 1
    %%% compute correlations

    sPairN = 0;

    % loop through task variables
    for s = 1:S

        % compute within-variable correlation (all correlations between 
        % each pair of resampled datasets)
        foo = corr(squeeze(BSSummary.sRA(s,:,:))'); % transpose so that datasets are in columns
        % take upper triangle
        BSSummary.corrWithin(s,1:W) = foo(logical(triu(ones(size(foo)),1)));
        clear foo

        % loop through remaining task-variables to compute between-variable
        % correlations
        for s2=s+1:S

            % track task variable pair number
            sPairN = sPairN + 1;

            % store task variable pair index
            BSSummary.corrBetweenInd(sPairN,:) = [s s2];

            % compute between signal correlation, take just the diagonal (i.e.,
            % correlation between signals from the same resampled dataset)
            BSSummary.corrBetween(sPairN,1:B) = ...
                diag(corr(squeeze(BSSummary.sRA(s,:,:))',...
                squeeze(BSSummary.sRA(s2,:,:))')); % transpose so that samples are in columns
        end
    end
else
    warning('Cannot compute correlations between resampled datasets because %d dataset requested',B)
end

%% compute summary statistics and hypothesis testing

% initialize
stats = [];

% only proceed if requested
if bitStats
    %%%%%%%%%%%%%%%%%%%%%%%
    % Compute predicted r values given variability in within-signal r and
    % assuming perfect between signal correlation.
    
    % See original derivation in Spearman, C. The Proof and Measurement of
    % Association between Two Things. Am. J. Psychol. 15,72-101 (1904). 

    % The following is based on assumptions that within-signal variability
    % is due to noise that is: 
    % 1. independent between neurons
    % 2. independent between signals
    % As such, for signals A and B with true variability across neurons
    % given by sA and sB and with independent noise eA and eB, the expected
    % variability across neurons for a given sample is sA + eA and sB + eB.
    % We measure sA as the std across neurons of the mean beta measured
    % across bootstrapped resample datasets. This assumes that the mean
    % across bootstraps is the true beta for each neuron.
    %
    % There are two ways we can measure e: 
    % 1. (favored method) We assume that the within signal correlation that
    % we measure, rXX, is a function of the independent noise eA, by:
    % rXX = rAA*sA^2/(sA^2 + eA^2)
    % Assuming rAA = 1 (this is the within signal correlation, afterall):
    % eA = sqrt(sA^2 * (1-rXX)/rXX).
    % 2. (alternative -- not used) We measure e directly as the mean
    % (across neurons) of the std of betas (across resampled datasets).
    % This assumes that the varibility across datasets is driven by
    % independent sampling noise, and that it is independent across
    % neurons, as such, we can take the mean across neurons as the true
    % variability of the independent noise.  
    %
    % NB: it turns out these two methods are VERY similar.
    %
    % Finally, the effect of within signal noise eA and eB on the measured
    % correlation, rABHat, between A and B assuming a true correlation,
    % rAB, of 1 (that is, measured without within-varialbe variability) is
    % given by:  
    %    rABHat = rAB * sA * sB / (sqrt(sA + eA) * sqrt(sB + eB))
    % To arrive at this, let A and B be multivariate normal (MVN)
    % distributions given by the sum of MVM([muA, muB], [sA^2, rAB*sA*sB;
    % rAB*sA*sB, sB^2]) and the noise distribution MVN ([0 0],[eA^2, 0; 0,
    % eB^2]).
    % By setting rAB=1, our null model assumes that the decrement in rABHat
    % is entirely driven by independent sampling noise (and not by true
    % separation/independence of the signals A and B).
        
    % we need to compute oTDR on the full dataset so as to
    % compute the veridical between-variable correlations
    fprintf('COMPUTING oTDR ON FULL DATASET...\n')
	otdr_full = runBS(false); % set resampleFlg = false
    
    % COMPUTE basic stats
    nSignal = size(BSSummary.corrWithin,1);
    nSPair = size(BSSummary.corrBetween,1);

    % loop through variables for within-variables stats
    stats.within.ci = [alpha*100 100];
    for i = 1:nSignal    
        % stats
        stats.within.r(i) = mean(BSSummary.corrWithin(i,:));
        stats.within.r_sigma(i) = std(BSSummary.corrWithin(i,:));
        stats.within.r_ci(i,:) = prctile(BSSummary.corrWithin(i,:),stats.within.ci);
        [~,stats.within.p(i)] = ttest(BSSummary.corrWithin(i,:),0,'alpha',alpha,'tail','right');
    end

    % loop through pairs of variables for between-variable stats
    stats.between.ci = 100*[alpha/2 (1-alpha/2)];
    for i = 1:nSPair
        stats.between.r(i) = mean(BSSummary.corrBetween(i,:));
        stats.between.r_sigma(i) = std(BSSummary.corrBetween(i,:));
        stats.between.r_ci(i,:) = prctile(BSSummary.corrBetween(i,:),stats.between.ci);
        [~,stats.between.p(i)] = ttest(BSSummary.corrBetween(i,:),0,'alpha',alpha);

        % distribution r_AB = sqrt(r_AA)*sqrt(r_BB)
        stats.nullModel.rPred_dist(i,:) = prod(sqrt(BSSummary.corrWithin(BSSummary.corrBetweenInd(i,:),:)));
        stats.nullModel.rPred_dist_mean(i) = mean(stats.nullModel.rPred_dist(i,:));
        stats.nullModel.rPred_dist_median(i) = median(stats.nullModel.rPred_dist(i,:));
        stats.nullModel.rPred_dist_CI(i,:) = prctile(stats.nullModel.rPred_dist(i,:),stats.within.ci);

        % measure r_AB from all trials
        [foo,goo] = corr(otdr_full.sRAStar.RA(:,BSSummary.corrBetweenInd(i,:)));
        stats.between.r_allTrials(i) = diag(foo,1);
        stats.between.p_allTrials(i) = diag(goo,1);

    end


    % compute alternative p-value that relies on a predicted
    % DISTRIBUTION of r_AB = sqrt(r_AA)*sqrt(r_BB), instead of a scalar
    % value. In this case, the measured scalar value r_AB can be compared
    % directly against the distribution of predicted r_AB. Note that the
    % predicted r_AB can be positive or negative. Rather that test both null
    % hypotheses, we assume that the predicted r_AB is the same sign as the
    % observed r_AB (which we implement by taking abs(observed r_AB)).
    % We also have a choice for observed r_AB. We can use the value across all
    % trials from oTDRSummary.sRAStar.RA, or we can use the mean value across
    % bootstrapped samples from stats.between.r. Since in this case we only
    % need the point estimate (and not the bootstrapped variance), it makes
    % sense to use the most accurate point estimate, which is that from all
    % trials (computed de novo as otdr_full, see above).
    [~,stats.nullModel.pVal_alt_ttest] = ...
        ttest(bsxfun(@minus,abs(stats.between.r_allTrials),stats.nullModel.rPred_dist'),[],'tail','left');
    for i = 1:length(stats.between.r_allTrials)
        stats.nullModel.pVal_alt_signrank(i) = ...
            signrank(abs(stats.between.r_allTrials(i)) - stats.nullModel.rPred_dist(i,:)',[],'tail','left');
    end
    
    %%%% THE FOLLOWING APPROACH IS DEFUNCT
    % % compute sA, sB
    % stats.nullModel.trueSigma = std(squeeze(mean(s.sRA,2)));
    % % compute eA, eB by favored method
    % stats.nullModel.noiseSigma = ...
    %     sqrt(stats.nullModel.trueSigma.^2 .* ...
    %     (1 - stats.within.r)./stats.within.r);
    % % compute eA, eB by alt method
    % stats.nullModel.noiseSigma_AltMethod = mean(squeeze(std(s.sRA,[],2)));
    % % compute predicted value of sigma based on sqrt(sA^2 + eA^2)
    % stats.nullModel.sigmaPred = ...
    %     sqrt(stats.nullModel.trueSigma.^2 + stats.nullModel.noiseSigma.^2);
    % stats.nullModel.sigmaPred_AltMethod = ...
    %     sqrt(stats.nullModel.trueSigma.^2 + stats.nullModel.noiseSigma_AltMethod.^2);

    % stats.nullModel.rPred = NaN(1,nSPair);
    % stats.nullModel.pVal_alt = NaN(1,nSPair);

    % stats.nullModel.rPred_AltMethod = NaN(size(stats.nullModel.sigmaPred));

    % ind = 0;
    % for i = 1:nSignal
    %     for j = i+1:nSignal
    %         ind = ind+1;
    %         % compute predicted value of r, rABHat
    %         % NEW METHOD based on simplified derivation of rPred = sqrt(rAA) *
    %         % sqrt(rBB)  -- see proof in paper/methods.
    %         stats.nullModel.rPred(ind) = prod(sqrt(stats.within.r([i j])));
    %         
    % %         % compute predicted value of r, rABHat
    % %         stats.nullModel.rPred(ind) = 1 * prod(stats.nullModel.trueSigma([i j])) / ...
    % %             prod(stats.nullModel.sigmaPred([i j])); % "1" refers to the assumed underlying correlation
    % %         stats.nullModel.rPred_AltMethod(ind) = 1 * prod(stats.nullModel.trueSigma([i j])) / ...
    % %             prod(stats.nullModel.sigmaPred_AltMethod([i j])); % "1" refers to the assumed underlying correlation
    % 
    %     end
    % end

    % % compute t-stat that tests whether the predicted r-value for
    % % between-signal correlation is significantly greater than the distribution
    % % of r values measured.
    % [~,stats.nullModel.pVal_old] = ...
    %     ttest(bsxfun(@minus,BSSummary.corrBetween',stats.nullModel.rPred),[],'tail','left'); % transpose so that signals are in column    
    
end

% store
BSSummary.stats = stats;

%% store various meta data

meta.cond2include = cond2include;
meta.codedParams = codedParams;
meta.dataset_times = dataset_times;
meta.codedParams4Dataset = codedParams4Dataset;
meta.spcConstraint = spcConstraint;
meta.param4nrm = param4nrm;

BSSummary.meta = meta;

%% NESTED FUNCTION FOR CREATING RESAMPLED DATASETS ON RUNNING oTDR FOR EACH 
    function oTDRSummary = runBS(varargin)
        
        if length(varargin) > 0
            resampleFlg = varargin{1};
        else
            % by default, we resample
            resampleFlg = true;
        end
        
        tic
        
        % compute trial-average PTSHs, with trial-level resampling with
        % replacement. Only conditions with task variable values matching
        % those in cond2include will be included.
        N_trialavg = getTrialAvgPSTHs(N, cond2include, resampleFlg);
        
        % convert to "Data" struct (i.e., one struct element per condition)
        [Data] = genDataStruct(N_trialavg, cond2include); 
                
        % orthFlg is ignored by oTDR in this instance, since the bootstrap
        % analysis depends on the non-orthogonalized sRAStar. However, the
        % internal call to oTDR() expects this input arg and expects it in
        % a specific format, so we create a dummy orthFlg here:
        orthFlg = cell(size(codedParams4Dataset));
        for i = 1:numel(codedParams4Dataset)
            orthFlg{i} = ones(1,length(codedParams4Dataset{i}));
        end
        
        [oTDRSummary] = oTDR(Data, codedParams,...
            'dataset_times', dataset_times, ...
            'codedParams4Dataset', codedParams4Dataset, 'orthFlg', orthFlg,...
            'bit4Bootstrap', true, 'bitSerialOTDR', false,...
            'param4nrm', param4nrm, 'spcConstraint', spcConstraint);

        toc
        
    end % END nested function
end % end main functino