% function [Summary] = oTDR(Data, codedParams, [name,value])
%
% INPUTS:
% Data -- 1 x C struct with one element for each of C conditions as
%       returned by genDataStruct().
%
% codedParams -- C x P matrix where each row is a unique condition and each
%       column is a task variable (i.e. regression predictor) for which an
%       sRA will be computed. Number of conditions C must match number of
%       elements in Data (see above). (NOTE: scaling to [0, 1] occurs
%       within function.)
%
% OPTIONAL INPUTS (passed as name, value pairs):
% params_from_TDR -- struct output from TDR(). Serves as input to
%       getTimes4oTDR() to specify time bins in which to compute the sRAs.
%       If param_from_TDR contains field "dataset_times", these will be
%       used to define the time bins. Otherwise, magThresh (which depends on
%       values stored within param_from_TDR) or, if empty, tMinMax will be
%       used to define temporal periods (see below). If empty, and
%       dataset_times not provided as separate input (see below), then
%       TDR() will be called automatically using numPCs4TDR (see below).
% numPCs4TDR -- number of PCs into which to project the neural data prior
%       to finding dRAs (NOT sRAs). This functions as a denoising step.
%       Leave empty [] to include all PCs, i.e., no denoising. NOTE: this
%       applies only to the call to TDR (see above), which occurs only when
%       params_from_TDR and dataset_times are NOT provided.
% dataset_times -- 1 x E cell array, where each cell specifies the time
%       bins to include for each of E datasets in which to compute the
%       sRAs. Time bins are specified as a vector and reference the time
%       bin centers provided in Data.times. Overrides dataset_times
%       specified in params_from_TDR (see above). If empty, will attempt to
%       determine dataset times from params_from_TDR or getTimes4oTDR()
%       (see above).
% tMinMax -- R x 2 matrix of minimum start time (Col 1) and maximum end
%       time (Col 2) for inclusion of data in each of R datasets (rows).
%       Passed as input to getTimes4oTDR(). Leave empty to use predefined
%       temporal epochs, if any, defined in params_from_TDR (see above).
% magThresh -- scalar specifying threshold of signal magnitude use to
%       determine times to include when defining each of R datasets. Note,
%       when magThresh is supplied, a separate (possibly overlapping)
%       dataset is defined, and used for computing the sRA, for each coded
%       param in codedParams. When magThresh is provided, it supercedes any
%       values provided in tMinMax. 
% codedParams4Dataset -- R x 1 cell array specifying the index of the
%       parameters to include in the model for each of R datasets (as
%       defined by tMinMax, magThresh, or params_from_TDR, see above). Each
%       cell corresponds to 1 of R pre-defined datasets. Each cell contains
%       a 1 x P row vector specifying the column index of the coded
%       parameter (in codedParams) to include in the model (P < C, number
%       of coded parameters). To not include any parameters for a given
%       dataset, leave the cell empty. If codedParams4Dataset is empty,
%       will compute all parameters in the first dataset.
% orthFlg -- R x 1 cell array specifying the coded parameters for which the
%       corresponding static regression axes (sRA) should be
%       orthogonalized. Each cell corresponds to the cells in
%       codedParams4Dataset and contains a 1 x P row vector (must match
%       size of vector in corresponding cell of codedParams4Dataset) of 0s
%       and 1s specifying whether (1) or not (0) to orthogonalize the sRA
%       pertaining to the corresponding coded parameter. If orthFlg empty
%       or not provided, will orthogonalize all sRAs. 
% bit4Bootstrap -- logical on whether to run oTDR solely for the purpose of
%       generating bootstrapped data for purpose of calculating within- and
%       between-signal reliability of non-orthogonalized sRAs. When TRUE,
%       analyses are limited to computing sRAStar. All other RAs and
%       variance analyses are omitted. Note that subsampling must be
%       handled at a higher level and the subsampled data passed in DATA.
%       Default = FALSE.
% bitSerialOTDR -- logical on whether to run oTDR independently at each
%       time bin. Assumes single dataset at each time step including all
%       predictors and all predictors are orthogonalized. Only oTDR on
%       orthogonalized sRAs is returned. Default==0, in which case, input
%       datasets, predictors and orthogonalization are observed.
% regressBins -- scalar specifying number of consecutive time bins to
%       combine into a single regression time bin when computing serial
%       oTDR. Must be provided when performing serial oTDR. Ignored when
%       bitSerialOTDR == false.
% param4nrm -- structure with terms for z-score normalization +/-
%       common-condition subtraction. If any of the fields are not
%       provided, they will be computed from the input data. 
%           .mu -- 1 x N vector of mean response for N neurons
%           .nrm -- 1 x N vector of standard deviation of response for N
%               neurons
%           .commonConds -- T x N matrix of mean response across conditions
%           for T time points and N neurons.
% vect4Angle -- N x P matrix specifying the N dimensions of P vectors. All
%       angles between the provided vectors and the discovered
%       non-orthogonalised sRAs (i.e., sRAStar) are measured with
%       associated p-values using the RARand (biased vector) null model.
%       (N dimenions must match the N dimensions provided in R.)
% bitExcludeCondition -- C x 1 logical vector on whether to exclude given
%       condition (referencing conditions defined in Data) when computing
%       the sRAs. Projection of the excluded condition onto the sRAs will
%       nonetheless be performed and reported separately in fields with
%       suffix "_excluded". Likewise, the task variable values for the
%       excluded conditions are returned in field "codedParams_excluded"
%       and scaled/centered by parameters used for the included conditions
%       (which did not take excluded conditions into account). For this
%       reason, the noramlized predictor values for an excluded condition
%       may be outside of the range [0, 1]. 
% numPCs4OTDR -- integer specifying number D of top PCs: the sRA will be
%       constrained to live in the space spanned by the top D PCs. The PCs
%       are computed via PCA on the firing rate response matrix of neurons
%       x times*conditions. When empty, no constraint is applied.
% spcConstraint -- matrix (N x V) defining V orthogonal, N-dimensional
%       vectors, the space spanned by which will constrain the sRAs. Leave
%       empty to apply no constraint. When numPCs4OTDR is provided (see
%       above), spcConstraint will be computed internally and an error will
%       be throw if spcConstraint is also provided as an input.
%
%%% TODO: NEED TO REVIEW AND UPDATE OUTPUT DOCUMENTATION
% oTDR Output:
% -Summary
%   .data: data matrix (num conditions and times x number of neurons)
%   .sRA: static vectors of interest (from oTDR)
%       .sRA1: Benefit vector of interest.
%       .sRA2: Choice vector of interest.
%       .sRA3: Expected Reward vector of interest.
%   .Times:
%       .all_times: times w r t offer of the data
%       .t1: Benefit vector (sRA1) time
%       .t2: Choice vector (sRA2) time
%       .t3: Expected reward vector (sRA3) time
%   .varAnalysis: variance significance analysis
%   .varAnalysis: variance significance analysis
%       .V_sRA1: variance across conditions at all times of projections on dRA1 vectors of size (all_times x 1)
%       .V_sRA2: variance across conditions at all times of projections on dRA2 vectors of size (all_times x 1)
%       .V_sRA3: variance across conditions at all times of projections on dRA3 vectors of size (all_times x 1)
%       .VE_sRA1: percent variance explained across conditions at all times of projections on dRA1 vectors of size (all_times x 1)
%       .VE_sRA2: percent variance explained across conditions at all times of projections on dRA2 vectors of size (all_times x 1)
%       .VE_sRA3: percent variance explained across conditions at all times of projections on dRA3 vectors of size (all_times x 1)
%       .RSV_sRA1: relevant benefit variance across conditions at all times of projections on dRA1 vectors of size (all_times x 1)
%       .RSV_sRA2: relevant choice variance across conditions at all times of projections on dRA2 vectors of size (all_times x 1)
%       .RSV_sRA3:  relevant expected reward variance across conditions at all times of projections on dRA3 vectors of size (all_times x 1)
%       .RSVE_sRA1: percent relevant benefit variance across conditions at all times of projections on dRA1 vectors of size (all_times x 1)
%       .RSVE_sRA2: percent relevant choice variance across conditions at all times of projections on dRA2 vectors of size (all_times x 1)
%       .RSVE_sRA3: percent relevant expected reward variance across conditions at all times of projections on dRA3 vectors of size (all_times x 1)
%       .ISV_sRA1: irrelevant benefit variance across conditions at all times of projections on dRA1 vectors of size (all_times x 1)
%       .ISV_sRA2: irrelevant choice variance across conditions at all times of projections on dRA2 vectors of size (all_times x 1)
%       .ISV_sRA3:  irrelevant expected reward variance across conditions at all times of projections on dRA3 vectors of size (all_times x 1)
%       .ISVE_sRA1: percent irrelevant benefit variance across conditions at all times of projections on dRA1 vectors of size (all_times x 1)
%       .ISVE_sRA2: percent irrelevant choice variance across conditions at all times of projections on dRA2 vectors of size (all_times x 1)
%       .ISVE_sRA3: percent irrelevant expected reward variance across conditions at all times of projections on dRA3 vectors of size (all_times x 1)
%       .zV_sRA1: z-scored deviation from chance of variance across conditions at all times of projections on dRA1 vectors of size (all_times x 1)
%       .zV_sRA2: z-scored deviation from chance of variance across conditions at all times of projections on dRA2 vectors of size (all_times x 1)
%       .zV_sRA3: z-scored deviation from chance of variance across conditions at all times of projections on dRA3 vectors of size (all_times x 1)
%       .zRSV_sRA1: z-scored deviation from chance of relevant benefit variance across conditions at all times of projections on dRA1 vectors of size (all_times x 1)
%       .zRSV_sRA2: z-scored deviation from chance of relevant choice variance across conditions at all times of projections on dRA2 vectors of size (all_times x 1)
%       .zRSV_sRA3:  z-scored deviation from chance of relevant expected reward variance across conditions at all times of projections on dRA3 vectors of size (all_times x 1)
%       .zISV_sRA1: z-scored deviation from chance of irrelevant benefit variance across conditions at all times of projections on dRA1 vectors of size (all_times x 1)
%       .zISV_sRA2: z-scored deviation from chance of irrelevant choice variance across conditions at all times of projections on dRA2 vectors of size (all_times x 1)
%       .zISV_sRA3:  z-scored deviation from chance of irrelevant expected reward variance across conditions at all times of projections on dRA3 vectors of size (all_times x 1)
%
% Copyright Gamal Elsayed and Daniel Kimmel, 2016

function [Summary] = oTDR(Data, codedParams, varargin)

%% define default values for optional input params

params_from_TDR = [];
dataset_times = [];
tMinMax = [];
magThresh = [];
codedParams4Dataset = [];
orthFlg = [];
bit4Bootstrap = false;
bitSerialOTDR = false;
regressBins = [];
param4nrm = [];
vect4Angle = [];
bitExcludeCondition = false(length(Data),1);
numPCs4TDR = [];
numPCs4OTDR = [];
spcConstraint = [];

%% collect optional name,value pairs
warnopts(assignopts(who, varargin));


%% checks

% number of coded parameter
nCodedParamsMax = size(codedParams,2);

if isempty(codedParams4Dataset)
    warning('codedParams4Dataset was empty or not provided. Will compute sRAs for all codedParams in the first dataset.');
    codedParams4Dataset = {1:size(codedParams,2)};
end

if isempty(orthFlg)
    warning('orthFlg was empty or not provided. Will orthgonalize all sRAs across all datasets.');
    orthFlg = cell(size(codedParams4Dataset));
    for i = 1:numel(codedParams4Dataset)
        orthFlg{i} = ones(1,length(codedParams4Dataset{i}));
    end
end

if ~isequal(length(codedParams4Dataset),length(orthFlg))
    error('Number of datasets referenced by codedParams4Dataset and orthFlg must be equal')
end

for i = 1:length(codedParams4Dataset)
    if length(codedParams4Dataset{i}) > nCodedParamsMax || any(~ismember(codedParams4Dataset{i},[1:nCodedParamsMax]))...
            || length(codedParams4Dataset{i}) ~= length(unique(codedParams4Dataset{i}))
        error('For a given dataset (i.e., cell element), codedParams4Dataset can only reference coded parameters from 1 to %d, and can reference each parameter at most once', nCodedParamsMax)
    end
    
    if length(codedParams4Dataset{i}) ~= length(orthFlg{i})
        error('For each dataset (i.e., cell element), number of predictors included by codedParams4Dataset must match number of flags in orthFlg')
    end
    
    if size(codedParams4Dataset{i},1) > size(codedParams4Dataset{i},2) || ...
        size(orthFlg{i},1) > size(orthFlg{i},2)
        error('Vectors in each cell of codedParams4Dataset and orthFlg must be ROW vectors')
    end

end

%% preprocessing 

% separate excluded conditions, if any
DataExcluded = Data(bitExcludeCondition);
Data = Data(~bitExcludeCondition);
codedParamsExcluded = codedParams(bitExcludeCondition,:);
codedParams = codedParams(~bitExcludeCondition,:);

% Normalize neural responses and remove common-condition response. When
% param4nrm is not empty, data will be normalized by predefined values
% within param4nrm; otherwise values are computed directly from included
% data.
[dataTensor, preprocessingSummary] = preprocess4TDR(data2Tensor(Data), param4nrm);
[T, N, C] = size(dataTensor);

% process data for excluded conditions (apply normalization used above)
if ~isempty(DataExcluded)
    dataTensorExcluded = preprocess4TDR(data2Tensor(DataExcluded), preprocessingSummary);
else
    dataTensorExcluded = [];
end

% scale and center coded params to range [0, 1]
codedParams_min = min(codedParams);
codedParams_range = range(codedParams);
codedParams = bsxfun(@minus, codedParams, codedParams_min);
codedParams = bsxfun(@times, codedParams, 1./codedParams_range);
% scale and center coded params for the excluded conditions by same factors
% as the included conditions. Note this may cause values for the excluded
% conditions to be outside of range [0, 1]
codedParamsExcluded = bsxfun(@minus, codedParamsExcluded, codedParams_min);
codedParamsExcluded = bsxfun(@times, codedParamsExcluded, 1./codedParams_range);

% store trial counts by condition by unit
trCountmtx = vertcat(Data.trialRept)';

% store matrix of unnormalized values of the task variables for each
% condition (NOTE: may not include same parameters as defined in
% codedParams)
Predictors = vertcat(Data.Predictors);
PredictorsExcluded = vertcat(DataExcluded.Predictors);

% store times for each time bin 
all_times = unique(vertcat(Data.times),'rows');
% confirm that all conditions have same time points
if size(all_times,1) ~= 1
    error('Time points must be the same across all conditions passed in Data');
end 

%% Constrain solution to lie on top PCs.

% only compute space constraint if not provided as input arg
if isempty(spcConstraint)
    if isempty(numPCs4OTDR)
        spcConstraint = [];     % Pass in top PCs here if the desired projections are constrained to be within the top PCs
    else    
        XN = reshape(permute(dataTensor,[1 3 2]), [], N);
        PCs = pca(XN);
        numPCs4OTDR = min(numPCs4OTDR, N);
        spcConstraint = PCs(:, 1:numPCs4OTDR);  % Pass in top PCs here if the desired projections are constrained to be within the top PCs
    end
elseif ~isempty(numPCs4OTDR)
    % check that numPCs4OTDR not also specified
    error('Cannot provide both a specific space contraint (spcConstraint) and a request to contrain the solution to a given number of PCs (numPCs4OTDR). Pick one.')
end

%% decide to run once for input param set, or run independently for each time point

if bitSerialOTDR

    % set all coded params for each dataset
    codedParams4Dataset = {[1:size(codedParams,2)]};
        
    % set orthogonalization flags
    orthFlg = {true(1,size(codedParams,2))};

    % reset certain values 
    tMinMax = NaN;
    magThresh = NaN;

    % combine time bins into regression bins so as to have more accurate
    % estimates of firing rate (same as in TDR.m)
    regressTimes = all_times(1:floor(T/regressBins)*regressBins);
    regressTimes = round(mean(reshape(regressTimes, regressBins, floor(T/regressBins)),1), 2);
    % num of regression times
    rT = length(regressTimes);
    
    % find corresponding regression time for each all_time
    posRT = interp1(regressTimes,[1:length(regressTimes)],all_times,'nearest','extrap');
    
    
    % generate random vectors based on dataspace across all times and pass
    % these to the function
    bigA = reshape(permute(dataTensor,[1 3 2]), [], N);
    numSamples = 10000;
    [~,RA_rand] = sampleRandSubspaces(1, cov(bigA), [], numSamples);     
    RARand = NaN(N,numSamples);
    for s = 1:numSamples
        RARand(:,s) = RA_rand{s}{:};
    end
    clear RA_rand
    
    % initialize
    clear srl 
    errPerIter = cell(rT,1);
    bitErrPerTime = false(rT,1);
    optimGradFinal = zeros(rT,1);
        
    % repeat oTDR independently for each time point
    for s = 1:rT
        tic
        fprintf('PROCESSING SERIAL OTDR TIME NUM %d \n',s);
        
        % make single dataset
        R = {dataTensor(posRT==s,:, :)};
        % force single times
        dataset_times = {all_times(posRT==s)};
        
        % special common condition input for reversing common-condition
        % subtraction for projections
        commCond4Proj = preprocessingSummary.commonConds(posRT==s,:);
        
        % call oTDR
        try
            [srl(s),errPerIter{s},optimGradFinal(s)] = runOTDR(R,codedParams,dataset_times,codedParams4Dataset,...
                orthFlg,dataTensor(posRT==s,:,:),dataTensorExcluded(posRT==s,:),[],commCond4Proj,...
                trCountmtx,bitSerialOTDR,bit4Bootstrap,RARand,vect4Angle,spcConstraint);
        catch
            bitErrPerTime(s) = true;
            warning('Error running serial oTDR for regression time num %d (%0.2f s)',s,mean(dataset_times{:}));
%             disp(lasterr);
        end
        toc
    end

    %%%
    % Compile certain fields across serial runs
    %%
    %%% RAs
    % instantiate
    nSRA = size(srl(1).sRA.RA,2);
    Summary.sRA.RA = NaN(rT,N,nSRA);
    % loop through times
    for s = 1:length(srl)
        % skip if error
        if bitErrPerTime(s)
            continue
        end
        Summary.sRA.RA(s,:,:) = srl(s).sRA.RA; 
    end
    
    %%% VARANALYSIS
    % loop through all varAnalysis fields
    fn = fieldnames(srl(1).sRA.varAnalysis);
    for f = 1:length(fn)
        % size of field
        sz = size(srl(1).sRA.varAnalysis.(fn{f}));
        
        if sz(1) == regressBins
            % if first dimension is size regressBins, assume field varies with time
            % instantiate
            Summary.sRA.varAnalysis.(fn{f}) = NaN([T,sz(2:end)]);
            
            % loop through times
            for s = 1:length(srl)
                % skip if error
                if bitErrPerTime(s)
                    continue
                end
                % store data
                foo = repmat({':'},1,ndims(srl(1).sRA.varAnalysis.(fn{f}))-1);
                Summary.sRA.varAnalysis.(fn{f})((s-1)*regressBins+[1:regressBins],foo{:}) = ...
                    srl(s).sRA.varAnalysis.(fn{f});
            end
   
        else
            % otherwise assume it is a repeating field and store only once
            Summary.sRA.varAnalysis.(fn{f}) = srl(1).sRA.varAnalysis.(fn{f});
        end
    end
    
    %%% projections
    % get highest level field names
    fn = fieldnames(srl(1).sRA);
    % loop through fields
    for f = 1:length(fn)
        if strncmp('proj',fn{f},4)
            % size of field
            sz = size(srl(1).sRA.(fn{f}));
            
            % if projection, instantiate field
            Summary.sRA.(fn{f}) = NaN([T,sz(2:end)]);

            % loop through times
            for s = 1:length(srl)
                % skip if error
                if bitErrPerTime(s)
                    continue
                end
                % store data
                foo = repmat({':'},1,ndims(srl(1).sRA.(fn{f}))-1);
                Summary.sRA.(fn{f})((s-1)*regressBins+[1:regressBins],foo{:}) = ...
                    srl(s).sRA.(fn{f});
            end
        elseif strcmp('t4RA',fn{f})
            % special case of times for RA

            % instantiate
            Summary.sRA.(fn{f}) = cell(1,nSRA);
            % loop through RAs
            for r = 1:nSRA
                Summary.sRA.(fn{f}){r} = [];
                
                % loop through times
                for s = 1:length(srl)
                    % skip if error
                    if bitErrPerTime(s)
                        continue
                    end
                    % store times
                    Summary.sRA.(fn{f}){r} = horzcat(Summary.sRA.(fn{f}){r},...
                        srl(s).sRA.(fn{f}){r});
                end

            end
            
        elseif ~ismember(fn(f),fieldnames(Summary.sRA))
            % Provided field has not been collected already, assume it is a
            % repeating field and store only once
            Summary.sRA.(fn{f}) = srl(1).sRA.(fn{f});
%         else
%             % skip field
%             fprintf('Field %s will be skipped \n',fn{f})
        end
        
    end
 
    
else
    % run oTDR once based on input set of datasets and predictors
    
    % pick times for oTDR and select data based on these times.
    %
    % determine what user input was provided, if any
    if isempty(dataset_times)
        % if dataset_times not provided, move onto extracting times from
        % params_from_TDR
        params_temp = [];
        if isstruct(params_from_TDR)
            % if struct, assume it is the TDR output.
            params_temp = params_from_TDR;
        else
            % if not user provided, run TDR to get "params" struct        
            params_temp = TDR(Data, codedParams, numPCs4TDR, regressBins, false, false);
        end
        
        % confirm that times in params_temp matches those in Data
        if ~isequal(round(all_times,10), round(params_temp.Times.all_times,10))
            error('When using TDR output to determine dataset times, the time bin centers in TDR must match those provided in Data')
        end
        % confirm that number of original time bins to combine for
        % regression matches between TDR and oTDR
        if regressBins ~= params_temp.Times.regressBins
            error('When using TDR output to determine dataset times, the number of time bins to combine for regression must match between TDR and oTDR');
        end

        % determine time bins that define datasets in which to compute the sRAs
        dataset_times = getTimes4oTDR(params_temp,tMinMax,magThresh);
    end
    
    % preallocate cell array for each dataset
    R = cell(1,length(dataset_times));
    
    % extract data for each dataset
    for i = 1:length(dataset_times)
        R{i} = dataTensor(ismember(round(all_times(:),2),dataset_times{i}),:, :);
    end
        
    % dummy
    RARand = [];
    
    % call oTDR
    [Summary,errPerIter,optimGradFinal] = runOTDR(R,codedParams,dataset_times,...
        codedParams4Dataset,orthFlg,dataTensor,dataTensorExcluded,...
        preprocessingSummary.commonConds,preprocessingSummary.commonConds,...
        trCountmtx,bitSerialOTDR,bit4Bootstrap,RARand,vect4Angle,spcConstraint);
end

%% Store high level info:

if bitSerialOTDR
    Summary.RARand = RARand;
end

%%%% data and predictors
Summary.dataTensor = dataTensor;
Summary.Predictors = Predictors;
Summary.PredictorsExcluded = PredictorsExcluded;
Summary.codedParams = codedParams;
Summary.codedParamsExcluded = codedParamsExcluded;

%%%% meta info
Summary.meta.preprocessingSummary = preprocessingSummary;
% store input vars
Summary.meta.magThresh = magThresh; % will be empty if no threshold used
if isempty(magThresh)
    Summary.meta.tMinMax = tMinMax; % not all datasets are necessarily used.
else
    Summary.meta.tMinMax = []; % empty when threshold is used
end
Summary.meta.codedParams4Dataset = codedParams4Dataset;
Summary.meta.orthFlg = orthFlg;
if ~bit4Bootstrap
    if iscell(errPerIter)
        errPerIterCell = errPerIter;
        % find max number of iterations
        maxIter = 0;
        for i = 1:length(errPerIterCell)
            maxIter = max(maxIter,length(errPerIterCell{i}));
        end
        % instantiate
        errPerIter = NaN(length(errPerIterCell),maxIter);
        % populate
        for i = 1:length(errPerIterCell)
            errPerIter(i,1:length(errPerIterCell{i})) = errPerIterCell{i};
        end
    end
    Summary.meta.errPerIter4SRA = errPerIter;
end
Summary.meta.optimGradFinal = optimGradFinal;
Summary.meta.bitSerialOTDR = bitSerialOTDR;
Summary.meta.numPCs4OTDR = numPCs4OTDR;
Summary.meta.spcConstraint = spcConstraint;

%%%% TIMES
Summary.Times.all_times = all_times;
if bitSerialOTDR
    Summary.Times.regressTimes = regressTimes;
    Summary.Times.regressBins = regressBins;
end

end % END main function
