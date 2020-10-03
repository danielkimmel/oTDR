%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Daniel L. Kimmel and Gamaleldin F. Elsayed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script runs Optimal Targeted Dimensionality Reduction (oTDR) against
% the data from Kimmel et al., Nature Communications, 2020. 
%
% The analyses are run independently for each animal.
%
% Set the desired parameters in first code block 
%
% Analysis code is found at:
% https://github.com/danielkimmel/oTDR
%
% Test data are found at:
% https://github.com/danielkimmel/Kimmel_NatComm_2020
%%%%%%%%%

%% add oTDR to path 
%%%
%%% For this block to run, file must be executed from command line.
%%% Once path is loaded by calling startup.m in root directory of package,
%%% this block can be skipped
%%%

% name of current file (for later testing purposes)
curfilename = 'test_oTDR.m';

% path to oTDR package
otdr_package_path = '';
if isempty(otdr_package_path)
    otdr_package_path = input('Enter absolute path to oTDR package: ','s');
end
if isempty(otdr_package_path)
    error('Must provide path to oTDR package');
end

% get current directory 
curdir = pwd;

% get directory of current file (only works if file is being executed from
% command line, not if executing by line or code block)
filedir = fileparts(mfilename('fullpath'));

try
    % navigate to top-level directory of oTDR package
    cd(otdr_package_path);
    % run startup file there
    startup;
catch
    error('Run startup.m in root directory of oTDR package. Alternatively, execute %s from command line.',curfilename);
end
% return to original directory
cd(curdir);

% confirm that path was appropriately updated
if ~exist(curfilename,'file')
    error('MATLAB path was not successfully updated. Run startup.m in root directory of oTDR package. Alternatively, execute %s from command line.',curfilename);
end

%% set overall parameters

% define path to data files for public release. 
dataPath = '';
if isempty(dataPath)
    dataPath = input('Enter absolute path to the data repository (must be fetched separately): ','s');
end
if isempty(dataPath)
    error('Path to data repository (i.e. "dataPath") must be provided')
end
    
% specify monkey name: 'N' or 'K'
monkeyName = 'N'; 

% vector specifying the types of units to include in analysis: 1=single
% unit, 2=multi unit (i.e., poorly isolated single unit, 3=hash (cluster of
% many indistinguishable units).  
unitTypes2include = [1 2 3];  

% bin width in secs
binWidth = 0.100; 

% logical on whether to limit analysis to cells with singleton condition
bitSingleton = false;

% define whether "choice" is encoded symmetrically as [-0.5 0.5] (TRUE) or
% assymetrically as [0, 1] (FALSE), as in the main text. This determines
% whether expected reward (benefit * choice) varies with benefit for both
% accept and reject choices, or only for accept choices, respectively.
bitSymmetricChoiceEncoding = false; 

% set animal-specific parameters
switch monkeyName
    case 'N' 
        % tMinMax specifies the temporal epochs in which the sRAs are
        % computed.  
        % THESE ARE TIMES THAT WERE USED IN MAIN TEXT.
        tMinMax = [0 0.5; 0.5 Inf];
        
        % threshold time for 'late' rejections (rejections before/after
        % this time will be classified as 'early'/'late'). Currently using
        % median rejection time across all neurons and conditions.
        lateRejectThresh = 0.92; 
        
        % number of PCs spanning the subspace into which the dRAs are
        % restricted. 
        numPCs = 8;
        
        % timing parameters for present offer epoch
        t_bef = 0.6; % time (s) to include before offer onset
        t_after = 4.6; % time (s) to include after offer onset
        workDur = 0.5 + 4; % duration (s) of work period
        
        % timing parameters for previous offer epoch 
        t_befPP = 10.4; % time (s) to include before fixation on present trial
        t_afterPP = 0.8; % time (s) to include after fixation on present trial

    case 'K'
        % tMinMax specifies the temporal epochs in which the sRAs are
        % computed.  
        % THESE ARE TIMES THAT WERE USED IN MAIN TEXT.
        tMinMax = [0 0.5; 0.5 Inf];
        
        % threshold time for 'late' rejections (rejections before/after
        % this time will be classified as 'early'/'late'). Currently using
        % median rejection time across all neurons and conditions.
        lateRejectThresh = 1.3;

        % number of PCs spanning the subspace into which the dRAs are
        % restricted. 
        numPCs = 20;
        
        % timing parameters for present offer epoch
        t_bef = 0.6; % time (s) to include before offer onset
        t_after = 7; % time (s) to include after offer onset
        workDur = 0.5 + 6.5; % duration (s) of work period
        
        % timing parameters for previous offer epoch 
        t_befPP = 12.8; % time (s) to include before fixation on present trial
        t_afterPP = 0.8; % time (s) to include after fixation on present trial
        
    otherwise
        error('Monkey name %s not recogized',monkeyName)
end

% for previous offer epoch, animals share the same window
tMinMaxPP = [0 0.5];


%% LOAD data

% check dataPath exists
if isempty(dataPath)
    error('Must define dataPath in parameters section above');
end

% update with path to raw data
dataPath_raw = fullfile(dataPath,'raw_data');

if ~exist(dataPath_raw,'dir')
    error('dataPath "%s" was not found',dataPath_raw);
end

% load raw data file (contains struct "data_raw")
foo = fullfile(dataPath_raw,sprintf('data_monkey%s.mat',monkeyName));
if ~exist(foo,'file')
    error('Raw data file not found at "%s"',foo);
end
load(foo);
clear foo

% determine which cells are included
bitCellIncluded = data_raw.bitCell2Include; % predefined in raw data file
% further select cells based on unit type (specified above)
bitCellIncluded = bitCellIncluded & ismember(vertcat(data_raw.byUnit.unitType), unitTypes2include(:).'); 
% initialize
bitCellWithSingleton = false(size(bitCellIncluded));
% loop through raw data for each cell
for i = 1:length(data_raw.byUnit)
   % Some cells are excluded from the previous-trial analysis (i.e., Figure
   % 6 in main text), even though they may have beenincluded in the
   % present-trial analysis, because of animal-imposed long delays between
   % a significant number of trials. Here we exclude these cells from all
   % analyses so that all analyses share the same cell population. Note,
   % that in the present paper, this additional exclusion did not impact
   % the cells included, as all cells included in the present-trial
   % analysis were also included in the previous-trial analysis.
   bitCellIncluded(i) = bitCellIncluded(i) & (~data_raw.byUnit(i).bitExcludeUnitFromPrevTrialAnal); 
   
   % determine cells for which singleton condition was collected
   bitCellWithSingleton(i) = any(data_raw.byUnit(i).trialInfo.bitSingleton);
   
end

% optionally limit cells by those with singleton condition
if bitSingleton
    bitCellIncluded = bitCellIncluded & bitCellWithSingleton;
end

%%% define conditions to include in the analysis

% Each row specifies a unique condition. Columns are: 
% 1) number of rewards offered (Benefit),
% 2) choice to accept (1) or reject (0) offer (Choice)
% 3) whether the 8-reward offer was signaled by (1) single purple icon,
%    i.e., singleton condition, or instead (0) by 8 yellow icons.
cond2include = [];
cond2include(:,1) = data_raw.condIndDic(:,2); 
cond2include(:,2) = data_raw.condIndDic(:,1); 
cond2include(:,3) = zeros(size(cond2include(:,1))); % singleton condition not included

% include additional conditions for separate analysis of late rejections
cond2include_withLateReject = [cond2include(cond2include(:,2)==0,:) false(sum(cond2include(:,2)==0),1);
    cond2include(cond2include(:,2)==0,:) true(sum(cond2include(:,2)==0),1);
    cond2include(cond2include(:,2)==1,:) false(sum(cond2include(:,2)==1),1)];

%% calculate firing rates

% Generate PTSHs for each neuron on each trial, aligned to various trial
% events (each alignment produces separate struct output argument)
[N, NprevPres, NpresRespPrevCond, N_withLateReject, timeIntervals] = genNstructFromSpks(...
    data_raw.byUnit(bitCellIncluded), binWidth, round(t_bef./binWidth), ...
    round(t_after./binWidth), 'offer', round(t_befPP./binWidth), ...
    round(t_afterPP./binWidth),lateRejectThresh);

% Generate trial-average responses for each neuron and condition. For each
% output struct (e.g., "Data"), each element is a unique condition and
% trial-average responses are stored in field ".A" as time points X
% neurons.
%
% Present-trial epoch, aligned to offer. 
[Data] = genDataStruct(N, cond2include); 
% Same as above, but include condition when the 8-reward offer was signaled
% by single purple icon. 
[DataSingle] = genDataStruct(N, [8 1 1]); 
% Same as present-trial epoch, but with separate conditions for early and
% late rejections
[Data_withLateReject] = genDataStruct(N_withLateReject, cond2include_withLateReject); 
% Previous-trial epoch: aligned to fixation on new trial and
% extending retrospectively into previous trial.
[DataPP] = genDataStruct(NprevPres, cond2include); 
% Present trial epoch, but with conditions coded by previous offer
[DataPresRespPrevCond] = genDataStruct(NpresRespPrevCond, cond2include); 
% Sanity check: replace earlier logical vector on whether singleton
% condition was run for given cell with new logical vector computed based
% on trial count
bitCellWithSingleton = logical(DataSingle.trialRept);

%% Define targeted variables (i.e., "codedParams")

% must define coded params separately for each set of trial-average
% responses (as generated above):
[codedParams,~,paramLabel] = define_codedParams(Data, bitSymmetricChoiceEncoding);
codedParamsSingle = define_codedParams(DataSingle, bitSymmetricChoiceEncoding);
codedParamsPP = define_codedParams(DataPP, bitSymmetricChoiceEncoding);
codedParamsPresRespPrevCond = define_codedParams(DataPresRespPrevCond, bitSymmetricChoiceEncoding);

% generate labels for previous-trial variables
paramLabel_PP = paramLabel;
for i = 1:length(paramLabel_PP)
    if contains(paramLabel_PP{i},{'Benefit','Choice'},'IgnoreCase',true)
        paramLabel_PP{i} = ['Prev ',paramLabel_PP{i}];
    elseif strcmpi(paramLabel_PP{i},'Expected Reward')
        paramLabel_PP{i} = 'Experienced Reward';
    end
end

%% Perform TDR on present-trial epoch
% Compute dynamic regression axes (dRAs), as in Figure 5 of main text.
%
% NOTE: we refer to analyses of dRAs as "TDR" because they did not require
% a numerical optimization step. In constrast, we refer to analyses of the
% static regression axes (sRA) as "oTDR", since these did require numerical
% optimization. Nonetheless, all analyses are based on our method of
% "optimal targeted dimensionality reduction", and should not be
% confused with the earlier "targeted dimensionality reduction" method
% (Mante et al., Nature, 2013), on which our method expands.

% define the number of time bins to average together when computing the
% dRAs. E.g., if time bins are 100ms and regressBins = 2, then each dRA
% will be computed in consecutive, non-overlapping time bins of 200ms.
regressBins = 2; % as used in main paper
% regressBins = 1; % as used in Supp Fig 22
% regressBins = 5; % as used in Supp Fig 22

% define number of surrogate datasets to generate for hypothesis testing.
numSamples4Surrogates = 1000;

% logical on whether to do ridge regression (TRUE), as in main text, or not
% (FALSE), as in Supp Fig 22.  
loocvFlg = true; % logical on whether to do ridge regression (true), or no regularization (false) -- MAIN PAPER
% loocvFlg = false; % as in Supp Fig 22

% number of top PCs onto which to project data for purposes of noise
% reduction. Set to numPCs to inherit the parameter value specified above
% and used in main text. Set to empty [] to skip the noise reduction step
% (equivalent to using all PCs), as in Supp Fig 22.
numPCsTemp = numPCs; % MAIN PAPER
% numPCsTemp = []; % no denoising, as in Supp Fig 22

% compute dRAs
[TDRSummary] = TDR(Data, codedParams, numPCsTemp, regressBins, true, true,...
    'numSamples4Surrogates',numSamples4Surrogates,...
    'loocvFlg',loocvFlg);

% store additional info
TDRSummary.bitSymmetricChoiceEncoding = bitSymmetricChoiceEncoding;
TDRSummary.bitCellIncluded = bitCellIncluded;
TDRSummary.siteInfo = data_raw.siteInfo;

%% Perform TDR on previous-trial epoch
% Same as above but for previous-trial epoch. Not shown in main paper.

regressBins = 2;
numSamples4Surrogates = 1000;
loocvFlg = true; % logical on whether to do ridge regression (true), or no regularization (false)

[TDRSummary_PP] = TDR(DataPP, codedParamsPP, numPCs, regressBins, true, true,...
    'numSamples4Surrogates',numSamples4Surrogates,...
    'loocvFlg',loocvFlg);
TDRSummary_PP.bitSymmetricChoiceEncoding = bitSymmetricChoiceEncoding;
TDRSummary_PP.bitCellIncluded = bitCellIncluded;
TDRSummary_PP.siteInfo = data_raw.siteInfo;

%% Perform TDR on present-trial epoch using previous offer coding
% As in above two analyses, but here the neural responses from the present
% trial are regressed against the variables from the previous trial,
% thereby testing for encoding of the previous trial during the present
% trial. Not shown in main paper.

regressBins = 2;
numSamples4Surrogates = 1000;
loocvFlg = true; % logical on whether to do ridge regression (true), or no regularization (false)

[TDRSummary_presRespPrevCond] = TDR(DataPresRespPrevCond, codedParamsPresRespPrevCond, ...
    numPCs, regressBins, true, true,...
    'numSamples4Surrogates',numSamples4Surrogates,...
    'loocvFlg',loocvFlg);
TDRSummary_presRespPrevCond.bitSymmetricChoiceEncoding = bitSymmetricChoiceEncoding;
TDRSummary_presRespPrevCond.bitCellIncluded = bitCellIncluded;
TDRSummary_presRespPrevCond.siteInfo = data_raw.siteInfo;
%saveFig2Directory('results/TDRanalysis')
%[Summary] = bootstrapdRAs(data_raw, bitCellIncluded, binWidth, t_bef, t_after, cond2include, numPCs, regressBins);


%% Perform oTDR 
% Compute static regression axes (sRAs) for the present trial, as in
% Figures 3 and 4 of main text.
%
% NOTE: we refer to analyses of sRAs as "oTDR" because they required a
% numerical optimization step. In constrast, we refer to analyses of the
% dynamic regression axes (dRA) as "TDR", since these did not require
% numerical optimization. Nonetheless, all analyses are based on our method
% of "optimal targeted dimensionality reduction", and should not be
% confused with the earlier "targeted dimensionality reduction" method
% (Mante et al., Nature, 2013), on which our method expands.

% define the number of time bins to average together when computing the
% dRAs. E.g., if time bins are 100ms and regressBins = 2, then each dRA
% will be computed in consecutive, non-overlapping time bins of 200ms. This
% value is only used when performing Serial oTDR, i.e., computing sRAs
% independently in each time bin (see below). These serial sRAs are for
% comparison to the dRAs (see above), and thus the same value for
% regressBins should be used for both analyses.
regressBins = 2; 

% Generally, the time period(s) in which the sRAs are computed is defined
% by tMinMax (see above). However, one can optionally select for time(s)
% based on when the normalized magnitude of the corresponding dRA exceeds threshold
% magThresh, supported on [0, 1]. See getTimes4oTDR(). When empty, time
% periods from tMinMax will be used.
magThresh = []; 

% Cell array specifying which coded params (i.e., task-relevant variables)
% to include for which dataset. Each cell corresponds to a dataset (i.e.
% temporal epoch) as specified by tMinMax or magThresh (see above).
% Therefore cannot have more cells than datasets. However, one can specify
% fewer coded params than datasets, in which case the extra datasets are
% ignored.
% USED IN MAIN PAPER:
codedParams4Dataset = {[1];[2 3]}; 
% USED IN SUPP FIG 12:
% codedParams4Dataset = {[1];[1 2 3]}; % BENEFIT IN EARLY AND LATE EPOCHS
% codedParams4Dataset = {[1 2 3];[1 2 3]}; % ALL PREDICTORS IN EARLY AND LATE EPOCHS
% codedParams4Dataset = {[1 2 3]}; % SINGLE EPOCH FOR ALL PREDICTORS
% NOT USED IN PAPER:
% codedParams4Dataset = {[1];[2];[3]}; % SEPARATE EPOCHS FOR EACH PREDICTOR

% specify which coded params to orthogonalize across datasets. Cell array
% should mirror array in codedParams4Dataset (see above). For each
% predictor specified in codedParams4Dataset, included a 1 or a 0 in
% orthFlg if the sRA for the predictor is to be orthogonalized or not,
% respectively, with respect to ALL other sRAs (not just the sRAs within
% the same dataset/epoch).
% USED IN MAIN PAPER
orthFlg = {[1];[1 1]}; % ALL ORTHOGONALIZED
% USED IN SUPP FIG 12:
% orthFlg = {[1];[1 1 1]}; % BENEFIT EARLY AND LATE -- ALL ORTHOGONALIZED
% orthFlg = {[1];[0 1 1]}; % BENEFIT EARLY AND LATE -- ONLY SRAS OF INTEREST ORTHOGONALIZED
% orthFlg = {[1 1 1];[1 1 1]}; % ALL PREDICTORS EARLY AND LATE -- ALL ORTHOGONALIZED
% orthFlg = {[1 1 1]}; % SINGLE EPOCH FOR ALL PREDICTORS
% NOT USED IN PAPER:
% orthFlg = {[1];[1];[1]}; % SEPARATE EPOCHS FOR EACH PREDICTOR -- ALL ORTHOGONALIZED

bitSerialOTDR = false; % logical to run oTDR independently at each time bin

%%% COMPUTE
% define vector to exclude the singleton condition when computing the sRAs
bitExcludeCondition = [false(length(Data),1); true];

% call oTDR. Note that we concatenate the non-singleton and singleton data
% and associated codedParams. However, the singleton condition is excluded
% from analysis by the vector bitExcludeCondition
[oTDRTemp] = oTDR([Data DataSingle], [codedParams; codedParamsSingle],...
    'params_from_TDR', [], 'numPCs4TDR', numPCs, 'dataset_times', [],...    
    'tMinMax', tMinMax, 'magThresh', magThresh,...
    'codedParams4Dataset', codedParams4Dataset, 'orthFlg', orthFlg,...
    'bit4Bootstrap', false, 'bitSerialOTDR', bitSerialOTDR,...
    'regressBins', regressBins, 'param4nrm', [], 'vect4Angle', [],...
    'bitExcludeCondition', bitExcludeCondition, 'numPCs4OTDR', []);

% store additional info
oTDRTemp.bitCellIncluded = bitCellIncluded;
oTDRTemp.siteInfo = data_raw.siteInfo;
oTDRTemp.Times.timeIntervals = timeIntervals;
oTDRTemp.meta.bitSymmetricChoiceEncoding = bitSymmetricChoiceEncoding;

% rename summary struct
if bitSerialOTDR
    oTDRSummary_serial = oTDRTemp;
else
    oTDRSummary = oTDRTemp;
end

clear oTDRTemp  

%% project rejection-aligned responses onto sRAs
% Here we recompile the trial-average data alinged to the rejection
% (instead of the offer, as above) and project these responses onto the
% sRAs computed above.

% initialize
rejectionTime = [];

% find median rejection time, defined as rejection time - offer time.
% NOTE this median rejection time is only used to determine start and end
% times of epoch around the rejection. The sorting of early vs. late
% rejections is determined by lateRejectThresh, as set above.
for i = 1:length(data_raw.byUnit)
    if ~bitCellIncluded(i)
        continue
    end
    
    rejectionTime = [rejectionTime, ...
        data_raw.byUnit(i).trialInfo.rejectTime ...
        - data_raw.byUnit(i).trialInfo.offerTime];
end
medRejectionTime = nanmedian(rejectionTime);
clear rejectionTime

% calculate firing rates aligned to rejection. Window extends
% retrospectively by median rejection time + minimum pre-offer fixation
% duration (0.5 s), and extends prospectively by trial duration - median
% rejection time. Note this returns 1xN structures (one element
% per unit) with trials sorted by condition, N_reject has the standard
% conditions as above. N_reject_withLateReject further separates reject
% choices into early vs late rejections based on medRejectionTime.
[N_reject,~,~,N_reject_withLateReject] = genNstructFromSpks(data_raw.byUnit(bitCellIncluded), binWidth, ...
    ceil((medRejectionTime+0.5)./binWidth), ...
    ceil((workDur - medRejectionTime)./binWidth),'rejection', [], [], lateRejectThresh);

% compile individual units into "Data_reject" struct by condition,
% excluding "accept" conditions. 
[Data_reject] = genDataStruct(N_reject, cond2include(cond2include(:,2)==0,:)); 
% repeat the above step, except now include separate conditions for for
% early vs. late rejections (given by new fourth column of
% Data_reject_withLateReject.Predictors). We use this output in the next
% code block.
[Data_reject_withLateReject] = genDataStruct(N_reject_withLateReject, ...
    cond2include_withLateReject(cond2include_withLateReject(:,2)==0,:)); 

% project rejection aligned data onto sRAs. Note that we z-score normallize the
% data using the parameters for the entire dataset (all
% conditions, not just rejections).  We do not remove the common condition
% response since we are interested in the overall temporal course
% peri-rejection.
oTDRSummary.sRA.projRejectAlign = projData(preprocess4TDR(data2Tensor(Data_reject), oTDRSummary.meta.preprocessingSummary,false),oTDRSummary.sRA.RA);
% also must store times relative to rejection
oTDRSummary.Times.all_times_rejectAlign = N_reject(1).all_times;

%% project early and late rejection responses onto sRAs
% Here we project trial-average responses onto the sRAs, as above. However,
% now we separate rejection trials into early and late rejections (based on
% median rejection time) and compute these trial-averages before
% projection. Because of the further segregation, the trials available per
% condition is decreased, and we must further pare-down the data to ensure
% adequate trial counts

% Set params to determine which conditions to exclude based on trial
% counts. The selection process works identically as for the main analysis
% above, but here we can use different thresholds, as desired.
%
% Set min number of trials per condition per unit. Units-Conditions with
% fewer trials, are subject to exclusion depending on minPropNeuronPerCond
% below.
minTrialPerCond = 5; 
% Set min proportion of neurons meeting minTrialPerCond trial-count
% threshold (above) required to include condition. Conditions with less
% than minPropNeuronPerCond proportion of neurons meeting the trial
% count threshold will be excluded
minPropNeuronPerCond = 0.85; 

% store predictors that include early and late rejection conditions (early
% vs. late status given by 0 or 1, respectively, in new fourth column of
% Predictors_withLateReject):
oTDRSummary.Predictors_withLateReject = vertcat(Data_withLateReject.Predictors);

% summarize trial count per condition
oTDRSummary.trCount_withLateReject = vertcat(Data_withLateReject.trialRept);
  
%%% Here we apply the trial-count and unit-per-condition criteria defined
%%% above
%
% since fewer trials are available for early/late rejections, compute
% proportion of neurons meeting min trial count per condition for each
% condition and eliminate conditions below threshold
bitInclCond = 1 - sum(oTDRSummary.trCount_withLateReject<minTrialPerCond,2) ...
    ./ sum(isfinite(oTDRSummary.trCount_withLateReject),2) ...
    >= minPropNeuronPerCond;

if any(~bitInclCond)
    disp('The following conditions will be eliminated:')
    disp(oTDRSummary.Predictors_withLateReject(~bitInclCond,:));
end

% eliminate neurons not meeting min trial count per condition
bitInclUnit = ~any(...
    oTDRSummary.trCount_withLateReject(bitInclCond,:)...
    <minTrialPerCond,1);
if any(~bitInclUnit)
    sprintf('The following %d units will be eliminated:',sum(~bitInclUnit))
    disp(find(~bitInclUnit));
end    

% store temp data struct based on Data_withLateReject after removing
% eliminated conditions
DataTemp = Data_withLateReject(bitInclCond);

% loop through conditions in DataTemp, removing eliminated units
for c = 1:length(DataTemp)
    DataTemp(c).A = DataTemp(c).A(:,bitInclUnit);
    DataTemp(c).trialRept = DataTemp(c).trialRept(bitInclUnit);
end

% remove eliminated units from the original preprocessingSummary. This will
% be used when normalizing the responses below. 
preprocSummTemp = oTDRSummary.meta.preprocessingSummary;
preprocSummTemp.var_acrossCond = preprocSummTemp.var_acrossCond(:,bitInclUnit);
preprocSummTemp.var_acrossTime = preprocSummTemp.var_acrossTime(bitInclUnit,:);
preprocSummTemp.mu = preprocSummTemp.mu(bitInclUnit);
preprocSummTemp.nrm = preprocSummTemp.nrm(bitInclUnit);
preprocSummTemp.commonConds = preprocSummTemp.commonConds(:,bitInclUnit);

% project trial-average data for the newly defined conditions: early
% rejections, late rejections, accepts.
oTDRSummary.sRA.projLateReject = projData(preprocess4TDR(data2Tensor(DataTemp), ...
    preprocSummTemp,true),oTDRSummary.sRA.RA(bitInclUnit,:));

% update predictor and trial-count matrices to reflect pared-down
% conditions
oTDRSummary.Predictors_withLateReject = oTDRSummary.Predictors_withLateReject(bitInclCond,:);
oTDRSummary.trCount_withLateReject = oTDRSummary.trCount_withLateReject(bitInclCond,bitInclUnit);

% Here we perform the rejection-aligned analysis, as in the prior code
% block, but now with separate conditions for early vs. late rejections. 
% We begin with the Data_reject_withLateReject struct as computed in the
% earlier code block. First we exclude the same conditions that were just
% excluded for the offer-aligned analysis.
bitInclCond_reject = ismember(vertcat(Data_reject_withLateReject.Predictors),...
    oTDRSummary.Predictors_withLateReject,'rows');
DataTemp_rejectAlign = Data_reject_withLateReject(bitInclCond_reject);
% Likewise, we exclude the units that were excluded in the offer-aligned
% analysis
for c = 1:length(DataTemp_rejectAlign)
    DataTemp_rejectAlign(c).A = DataTemp_rejectAlign(c).A(:,bitInclUnit);
    DataTemp_rejectAlign(c).trialRept = DataTemp_rejectAlign(c).trialRept(bitInclUnit);
end

% project REJECTION-ALIGNED trial-average data from early and late
% rejection trials. Note that we z-score normalize the data using the
% parameters for all trials (not just the rejections) of the included units
% (per bitInclUnit) of the pared-down dataset.  Also, we do not remove the
% common condition response since we are interested in the overall temporal
% dynamics peri-rejection.
oTDRSummary.sRA.projRejectAlign_wLateReject = projData(preprocess4TDR(data2Tensor(DataTemp_rejectAlign), ...
    preprocSummTemp,false),oTDRSummary.sRA.RA(bitInclUnit,:));

clear Data_reject_withLateReject preprocSummTemp DataTemp

%% oTDR for Previous offer epoch
% Compute static regression axes (sRAs) for the previous trial, as in
% Figure 6 of main text.

% specify which coded params to include for which dataset. Note that for
% the previous trial analysis, all sRAs are computed within a single
% dataset/epoch, defined as the first 0.5s of fixation on the present
% trial. 
codedParams4DatasetPP = {[1 2 3]};
% specify which coded params to orthogonalize across datasets (all params
% are orthogonalized)
orthFlgPP = {[1 1 1]};

% logical to run oTDR independently at each time bin
bitSerialOTDR = false; 

%%% COMPUTATION
% For z-score normalization, apply normalization terms from present-trial
% analysis. This requires oTDRSummary from above.
param4nrm.mu = oTDRSummary.meta.preprocessingSummary.mu;
param4nrm.nrm = oTDRSummary.meta.preprocessingSummary.nrm;

% run oTDR
[oTDRTemp] = oTDR(DataPP, codedParamsPP, ...
    'numPCs4TDR', numPCs, ...    
    'tMinMax', tMinMaxPP, ...
    'codedParams4Dataset', codedParams4DatasetPP, 'orthFlg', orthFlgPP,...
    'bit4Bootstrap', false, 'bitSerialOTDR', bitSerialOTDR,...
    'regressBins', regressBins, 'param4nrm', param4nrm, ...
    'vect4Angle', oTDRSummary.sRAStar.RA, 'numPCs4OTDR', []);    

% save additional params
oTDRTemp.bitCellIncluded = bitCellIncluded;
oTDRTemp.siteInfo = data_raw.siteInfo;
oTDRTemp.meta.bitSymmetricChoiceEncoding = bitSymmetricChoiceEncoding;

% rename summary struct
if bitSerialOTDR
    oTDRSummary_PP_serial = oTDRTemp;
else
    oTDRSummary_PP = oTDRTemp;
end

clear oTDRTemp  param4nrm


%%% projections and var analysis between present and previous trial epochs
%%% (i.e., comparing sRAs in oTDRSummary (present trial) and oTDRSummary_PP
%%% (previosu trial))
if ~bitSerialOTDR

   % store temporal intervals of trial events
    oTDRSummary_PP.Times.timeIntervals = timeIntervals;
 
    % extract tensor of present-trial trial-average data, sorted by
    % conditions on present trial, normalized by mean and variance
    % estimated from present-trial responses, and with the common-condition
    % (CC) response subtracted.
    dataTensor = preprocess4TDR(data2Tensor(Data), oTDRSummary.meta.preprocessingSummary,true);
    
    % extract tensor of present-trial data BUT where the trial-average
    % responses are sorted by the conditions on the previous trial. Also,
    % as above, responses are normalized by mean and variance estimated
    % from present-trial responses, andthe common-condition (CC) response
    % is subtracted.
    dataTensorPresRespPrevCond = preprocess4TDR(data2Tensor(DataPresRespPrevCond), oTDRSummary.meta.preprocessingSummary,true);
    
    % extract tensor of responses from the previous trial (sorted by
    % previous-trial conditions), normalized by mean and variance estimated
    % from PREVIOUS-trial responses, and with the common-condition (CC)
    % response (also estimated from the previous trial) subtracted.
    dataTensorPP = preprocess4TDR(data2Tensor(DataPP), oTDRSummary_PP.meta.preprocessingSummary,true);

    %%%% project present-trial responses on previous-trial sRAs 
    % (we save the projections in the previous-trial struct
    % oTDRSummary_PP.)
    %
    % project present-trial data (sorted by present-trial conditions) onto
    % previous-trial sRAs. This tests whether the previous-trial sRAs
    % explain variance in the present-trial responses related to the
    % present-trial conditions. We expect they would not, and so this is in
    % effect a sanity check. (Note: no need to apply mean subtraction here,
    % as already applied in constructing data tensor)
    oTDRSummary_PP.sRA.projPresent = projData(dataTensor,oTDRSummary_PP.sRA.RA,true);
    % perform variance analysis of projections (Note that "codedParams"
    % input includes the scaled predictors used to compute the sRAs, not
    % the original predictors in the "Data" struct.)
    [oTDRSummary_PP.sRA.varAnalysis_present] = TDRvarAnalysis(oTDRSummary_PP.sRA.RA, dataTensor, oTDRSummary.codedParams, [], [], []);
    
    % Again we project present-trial responses onto the previous-trial
    % sRAs, but now the responses are sorted by the previous-trial
    % conditions. This tests whether the previous-trial sRAs continue to
    % represent information related to the previous-trial variables into
    % the present-trial (as aligned to the present-trial offer). This
    % analysis appears in the main text in Figure 6e-h, right-hand panels.
    oTDRSummary_PP.sRA.projPresRespPrevCond = projData(dataTensorPresRespPrevCond,oTDRSummary_PP.sRA.RA,true);
    % And perform variance analysis. See above note re: codedParams.
    [oTDRSummary_PP.sRA.varAnalysis_presRespPrevCond] = ...
        TDRvarAnalysis(oTDRSummary_PP.sRA.RA, dataTensorPresRespPrevCond, oTDRSummary_PP.codedParams, [], [], []);
    % store times of the present-trial epoch within the previous-trial
    % summary struct.
    oTDRSummary_PP.Times.all_times_present = oTDRSummary.Times.all_times;
    
    %%%% project previous-trial responses onto present-trial sRA 
    % (we save the projections with in the present-trial struct
    % oTDRSummary.)
    %
    % project previous-trial responses (sorted by previous-trial
    % conditions) onto present-trial sRAs. This tests whether the
    % present-trial sRAs explain variance related to the previous-trial
    % conditions in responses from the previous-trial. Note that during the
    % previous trial, the previous-trial conditions are in effect "present-
    % trial" conditions as they refer to trial events coincident with the
    % neural response. As such, we would expect the present-trial sRAs to
    % encode information about the coincident trial, regardless of whether
    % it was designated ?present? or ?previous.? (Note: no need to apply
    % mean subtraction, as already applied in constructing data tensor.)
    % This analysis appears in the main text in Figure 6a-d, left-hand
    % panels.
    oTDRSummary.sRA.projPP = projData(dataTensorPP,oTDRSummary.sRA.RA,true);
    % And perform variance analysis. See above note re: codedParams.
    [oTDRSummary.sRA.varAnalysis_PP] = TDRvarAnalysis(oTDRSummary.sRA.RA, dataTensorPP, oTDRSummary_PP.codedParams, [], [], []);
    % store times of the previous-trial epoch within the present-trial
    % summary struct.
    oTDRSummary.Times.all_times_PP = oTDRSummary_PP.Times.all_times;
    
    clear dataTensor dataTensorPP
end


%% Bootstrap resampling -- present-trial epoch

% This analysis computes non-orthogonalized present-trial sRAs (i.e.,
% "sRAStar") in separate datasets that have been constructed by randomly
% selecting trials with replacement from the original data. This
% distribution of sRAs across resampled (bootstrapped) datasets is used in
% the main text to: estimate reliability of sRA coefficients for individual
% units (i.e., error bars in Figure 3a,c), estimate the reliability of a
% given sRA across datasets (i.e., "Within-sRA reliability", Figure 3b,d),
% and estimate the null distribution of correlation between sRAs for
% different variables (i.e., "Between-sRA correlation", Figure 3b,d).

% Note that this analysis inherits the analysis parameters defined in the
% earlier oTDR code block above. 

% Set number of resampled datasets to generate. Set to 700 since this
% represents the max number of trials per experiment across both animals.
numBS = 700;

% Extract the per unit mean and variance from the present-trial epoch,
% which will be used to z-score normalize the responses from each resampled
% dataset. This ensures that all datasets share the same normalization
% relative to the original data. This step requires the oTDRSummary struct
% from the earlier code block.
param4nrm.mu = oTDRSummary.meta.preprocessingSummary.mu;
param4nrm.nrm = oTDRSummary.meta.preprocessingSummary.nrm;

% Keep N structs from above. Note that trial-average substructure ".cond"
% will be removed by bootstrapsRAs().

% Run boostrapped analysis.
[BSSummary] = bootstrapsRAs(N, cond2include, ...
    codedParams, oTDRSummary.sRA.t4RA, codedParams4Dataset, [],...
    'nBootstraps',numBS,'param4nrm',param4nrm);
  
%% Bootstrap resampling -- previous-trial epoch

% This analysis is identical to the above code block except computes
% non-orthogonalized PREVIOUS-trial sRAs. The results appear in
% Supplementary Figures 28 and 29.

% Note that this analysis inherits the analysis parameters defined in the
% earlier oTDR (previous trial) code block above. 

% Set number of resampled datasets to generate. Set to 700 since this
% represents the max number of trials per experiment across both animals.
numBS = 700;

% Extract the per unit mean and variance from the present-trial epoch,
% which will be used to z-score normalize the responses from each resampled
% dataset. This ensures that all datasets share the same normalization
% relative to the original data. This step requires the oTDRSummary struct
% from the earlier code block.
param4nrm.mu = oTDRSummary_PP.meta.preprocessingSummary.mu;
param4nrm.nrm = oTDRSummary_PP.meta.preprocessingSummary.nrm;

% Run boostrapped analysis.
[BSSummary_PP] = bootstrapsRAs(NprevPres, cond2include, ...
    codedParamsPP, oTDRSummary_PP.sRA.t4RA, codedParams4DatasetPP, [],...
    'nBootstraps',numBS,'param4nrm',param4nrm);

%% ANGLE PROFILE -- Stability analysis
% This block fits the "angle profile", specifically a boxcar step function,
% to the time-series of angles between pairs of dRAs. It generates a
% summary struct angleProfile that will be used in later code blocks. Each
% element of angleProfile contains the comparison between a single pair of
% dRAs, either of the same variable (e.g., benefit at time t1 vs benefit at
% time t2) or between different variables (e.g., benefit at time t1 vs
% choice at time t2). Because the first (reference) variable for pairs of
% different variables is not symmetric, there are actually two separate
% entries for each comparison of different variables: the non-transposed
% and transposed entries. As such, for 3 variables, there are a total of 9
% comparisons (3 same, 3 different non-transposed, 3 different transposed).
%
% Within each entry of angleProfile, there are separate fields for the
% different analyses. The fields starting with "main_step" refer to the
% standard analysis in Figure 5 (aka, "CIS" comparison) where we compare
% the "folded" or "reflected" angle between dRAs, such that the angles
% equidistant from 90 deg are treated identically, e.g. an aboslute angle
% of 80 or 110 are both transformed to 80 deg.
%
% Fields starting with "main_step_trans" refer to a special analysis in
% Kimmel et al., 2020, Supp Figure 23 (aka, "TRANS" comparison) where we
% compare the "unfolded" (non-transformed) angle between dRAs, and limit
% analysis to angles to range [90 180] deg. This effectively assess sign
% reversals, i.e., when the direction of encoding flips mid-trial from
% positive to negative, or vice versa. 
%
% Fields ending in "_PP" refer to dRAs computed during the PREVIOUS-trial
% epoch with regards to the previous-trial conditions. Fields ending in
% "_presRespPrevCond" refer to dRAs computed during the PRESENT-trial
% epoch with regards to the previous-trial conditions. (Note: neither of
% these analyses are presented in Kimmel et al., 2020.). Fields without
% either suffix refer to dRAs computed during the present-trial epoch with
% regards to the present-trial conditions, as shown in Kimmel et al., 2020.
% 
% Within each analysis, after the angle profile is fit to the timeseries of
% angles, the angle profile is applied to the angles computed in each of
% the surrogate datasets, and the mean angle (i.e., "similarity") during
% the period defined by the profile is extracted from each surrogate. In
% fields prefixed by "step_", various summary statistics are reported
% (e.g., mean, percentiles, etc). When field specifies "theta", it refers
% to the veridical data. When field specifies "thetaS", it refers to the
% null distribution of similarities from the surrogate data. In addition,
% the veridical similarity is compared to this null distribution and a
% p-value is reported.

%%%%%%%%% USER PARAMETERS

% cell array specifying temporal epochs to use for analysis. All entries
% will be computed. Options include
% 'present', 'previous', and 'present_response_previous_condition', which
% refer to fields with no suffix, "_PP", and "_presRespPrevCond",
% respectively (see above).
epoch2Use = {'present','previous','present_response_previous_condition'};

% logical on whether to fit angle profiles to each surrogate dataset,
% stored in field ".bS" (times x fit parameters x surrogate dataset). In
% standard analysis, this is omited (=0), since fits from veridical data
% are applied to surrogate data. Applying fits to surrogate data (=1) does
% not commit one to a particular analysis later, but is time consuming.
% Note that analysis in Kimmel et al., 2020, did NOT fit angle profiles to
% the surrogate data. Likewise, the analysis code to summarize these fits
% has NOT been included with the public release.
bitFitSurr = 0; 

% logical on whether to fit a 2-step (dual) boxcar function (instead of a
% single boxcar). This is not used in the standard analysis and extends the
% computation time.
bit2Step = false; 

% tolerance for error in gamma fits. (Gamma distributions are fit to the
% null distributions so as to compute smooth CDFs. We tolerate a certain
% amount of error in these gamma fits, as specified below.)
if strcmp('K',monkeyName)
    % Monkey K requires greater tolerance for step fits, but this is OK
    % because the higher error occurs only on 1 or 2 time bins and only in
    % the cross-correlograms
    errThreshSigTestGamma = 0.165;
    
    % separate error thresh for previous-trial epoch, which happens to be
    % the same for monkey K.
    errThreshSigTestGamma_PP = errThreshSigTestGamma; 
else
    % error thresh for present- and previous-trial epochs
    errThreshSigTestGamma = 0.11;
    errThreshSigTestGamma_PP = 0.35;
end
% need more tolerance for 2-step fits
if bit2Step
    errThreshSigTestGamma = errThreshSigTestGamma + 0.03;
    errThreshSigTestGamma_PP = errThreshSigTestGamma_PP + 0.03;
end



%%%%%%%%%% COMPUTATION

% set angle offset. This is the amount by which theta is offset as theta* =
% angleOffset - theta. 
% Set angleOffset to EMPTY to use the max angle for any given row of theta
% as angleOffset. Note in this case, the offset may be <90deg.
angleOffset = 90;

% (The following flag only matters for gradient descent fit method. Now we
% use the complete parameter search method, for which bitContrainOrder does
% not matter).
% logical on fit formula. We can either constrain the horizontal order of
% two terms (TRUE), in which case "c" refers to the width:
%     a*(1/(1+exp(-%f*(x-b))) - 1/(1+exp(-%f*(x-(b+c)))))
% Or we can constrain the endpoint of the pulse (FALSE), in which case "c"
% refers to the endpoint:
%     a*(1/(1+exp(-%f*(x-b))) - 1/(1+exp(-%f*(x-(c)))))
% In the future, if we really cared, we could use FMINCON() to constrain
% both "c" and "b+c".
bitConstrainOrder = false;

% number of repetitions of fit before determining best fit. Only applies to
% gradient decent method, which is not used presently.
nIter = 100; 

% loop through epochs
for e = epoch2Use
    
    % determine flags for epoch
    switch lower(e{:})
        case 'present'
            % whether to use previous offer TDR
            bitPP = false;

            % whether to use previous offer conditions on present trial responses
            % (superceded by bitPP)
            bitPRPC = false;
            
        case 'previous'
            bitPP = true;
            bitPRPC = false;
        case 'present_response_previous_condition'
            bitPP = false;
            bitPRPC = true;
        otherwise
            error('epoch "%s" was not recognized.', e{:});
    end
    
    % loop through "CIS" and "TRANS" analyses
    for ct = {'cis','trans'}

        switch ct{:}
            case 'cis'
                % logical on whether to reflect theta about 90deg, forcing
                % theta into range [0 90] 
                bitReflectTheta = 1; 

                % define range to include, inclusive. data OUTSIDE these
                % bounds will be considered outliers and excluded. Refer to
                % theta in original coordinates, before any transformation
                % applied below. If you wish to include all angles, set to
                % [-Inf Inf];
                dataRange2Incl = [0 90];

                % set scalar to apply to angles such that theta* =
                % angleScale*(angleOffset - theta). Generally, scalar will
                % be 1 (default) or -1 to invert the angles. This is useful
                % when attempting to fit angles > 90deg.
                angleScale = 1;
                
            case 'trans'
                % logical on whether to reflect theta about 90deg, forcing
                % theta into range [0 90] 
                bitReflectTheta = 0; 

                % define range to include, inclusive. data OUTSIDE these
                % bounds will be considered outliers and excluded. Refer to
                % theta in original coordinates, before any transformation
                % applied below. If you wish to include all angles, set to
                % [-Inf Inf];
                dataRange2Incl = [90 180];

                % set scalar to apply to angles such that theta* =
                % angleScale*(angleOffset - theta). Generally, scalar will
                % be 1 (default) or -1 to invert the angles. This is useful
                % when attempting to fit angles > 90deg.
                angleScale = -1;
        end
        
        %%%%%%%%%%%%%

        % gather TDR summary struct to use depending of whether using present or
        % previous offer 
        if bitPP
            tdrTemp = TDRSummary_PP;
            paramLabelTemp = paramLabel_PP;
            errThreshSigTestGamma_temp = errThreshSigTestGamma_PP;
        elseif bitPRPC
            tdrTemp = TDRSummary_presRespPrevCond;
            paramLabelTemp = paramLabel_PP;
            errThreshSigTestGamma_temp = errThreshSigTestGamma;
        else
            tdrTemp = TDRSummary;
            paramLabelTemp = paramLabel;
            errThreshSigTestGamma_temp = errThreshSigTestGamma;
        end
        
        % TODO: The following is a work-around for fact that
        % TDRstabilityAnalysis() is hard-coded to generate surrogate data
        % for fixed number of signal (i.e., dRA) pairs. In future release
        % that includes more general form of TDRstabilityAnalysis(), the
        % following code will need to be updated
        %
        % Find number of fields referring to angles. should include 3 auto
        % angles and 3 cross angles AND a repeat of the 3 cross angles so
        % we can measure them in both directions
        fn = fieldnames(tdrTemp.angleAnalysis);
        fieldI = strncmp('angle',fn,5);
        % cutdown fieldnames:
        fn = fn(fieldI);
        % extract numbers from field names:
        foo = regexp(fn,'(\d{1})(\d{1})','tokens');
        sigN = [str2num(cell2mat(cellfun(@(c) (c{1}{1}),foo,'UniformOutput',0))), ...
            str2num(cell2mat(cellfun(@(c) (c{1}{2}),foo,'UniformOutput',0)))];
        % sort sigN by auto then cross angles
        bitAuto = diff(sigN,1,2) == 0;
        sigN = [sortrows(sigN(bitAuto,:)); sortrows(sigN(~bitAuto,:)); sortrows(sigN(~bitAuto,:))];
        nSignal = size(sigN,1);
        % store key on whether signal pair is/should be transposed
        bitTranspose = [false(sum(bitAuto),1); false(sum(~bitAuto),1); true(sum(~bitAuto),1)];

        s2u = 'main_step';
        if strcmpi(ct,'trans')
            % append "trans" suffix"
            s2u = [s2u,'_trans'];
        end
        if bit2Step
            s2u = [s2u,'_2'];
        end
        % append for using previous offer
        if bitPP
            s2u = [s2u,'_PP'];
        elseif bitPRPC
            s2u = [s2u,'_presRespPrevCond'];
        end

        % define time range to include in fits as [min max] inclusive. 
        if bitPP || bitPRPC
            % previous trial -- include all times
            tRange2Incl = [-Inf Inf];
        else
            % present trial -- limit to after offer
            tRange2Incl = [0 Inf];
        end

        % instantiate 
        [angleProfile(1:nSignal).(s2u)] = deal(struct([]));

        tic
        for i = 1:nSignal
            
            %%%%%%%%%%%%%%%%%
            % FIT ANGLE PROFILES TO VERIDICAL DATA            

            % extract surrogate angles. 
            thetaS = tdrTemp.angleAnalysis.(['surrAngle',num2str(sigN(i,1)),num2str(sigN(i,2))]);
            
            % if not fitting profile to surrogates, set temporarily to
            % empty, which tells stabilityFitStep() not to fit profiles to
            % them.
            if bitFitSurr
                thetaS_temp = thetaS;                
            else
                thetaS_temp = [];
            end

            if ~bitTranspose(i)

                [angleProfile(i).(s2u)] = stabilityFitStep(...
                    tdrTemp.angleAnalysis.(['angle',num2str(sigN(i,1)),num2str(sigN(i,2))]),...
                    tdrTemp.angleAnalysis.regressTimes,...
                    'thetaS',thetaS_temp,...
                    'nIter',nIter,...
                    'dataRange2Incl',dataRange2Incl,...
                    'angleOffset',angleOffset,...
                    'angleOffsetMin',90,...
                    'angleScale',angleScale,...
                    'bitConstrainOrder',bitConstrainOrder,...
                    'bitReflectTheta',bitReflectTheta,...
                    'bit2Step',bit2Step,...
                    'tRange2Incl',tRange2Incl);

                % extract and store which signals are being compared:
                angleProfile(i).(s2u).signalN = sigN(i,:);

            else
            % for cross-signal angles, we have to test the stability both with
            % signal 1 being in the rows, and transposed such that signal 2 is in
            % the rows
                [angleProfile(i).(s2u)] = stabilityFitStep(...
                    tdrTemp.angleAnalysis.(['angle',num2str(sigN(i,1)),num2str(sigN(i,2))])',... % NOTE TRANSPOSE
                    tdrTemp.angleAnalysis.regressTimes,...
                    'thetaS',permute(thetaS_temp, [2 1 3]),... % NOTE PERMUTE
                    'nIter',nIter,...
                    'dataRange2Incl',dataRange2Incl,...
                    'angleOffset',angleOffset,...
                    'angleOffsetMin',90,...
                    'angleScale',angleScale,...
                    'bitConstrainOrder',bitConstrainOrder,...
                    'bitReflectTheta',bitReflectTheta,...
                    'bit2Step',bit2Step,...
                    'tRange2Incl',tRange2Incl);

                % extract and store which signals are being compared:
                angleProfile(i).(s2u).signalN = fliplr(sigN(i,:));
        %         % name it
        %         angleProfile(i).(s2u).signalDesc = sprintf('%s X %s',paramLabelTemp{sigN(i,2)},paramLabelTemp{sigN(i,1)});

            end
            
            clear thetaS_temp

            angleProfile(i).(s2u).bitTranspose = bitTranspose(i);

            % store description of signals
            % if using previous offer predictors, add prefix
            if bitPP || bitPRPC
                foo = 'Prev ';
            else
                foo = '';
            end
            angleProfile(i).(s2u).signalDesc = sprintf('%s%s',foo,paramLabel{angleProfile(i).(s2u).signalN(1)});
            % add second signal if cross signals
            if length(unique(angleProfile(i).(s2u).signalN)) > 1
                angleProfile(i).(s2u).signalDesc = [angleProfile(i).(s2u).signalDesc,' X ',sprintf('%s%s',foo,paramLabel{angleProfile(i).(s2u).signalN(2)})];
            end
            
            
            %%%%%%%%%%%%%%%%%
            % COMPUTE SIMILARITY IN SURROGATE DATA 
            % using profiles fit to veridical data above.

            % extract from above output
            theta = angleProfile(i).(s2u).theta; % veridical angles
            tSet = angleProfile(i).(s2u).tMaster; % still use tMaster (not t), since no special coversion happens here
            nTimes = length(tSet);
            nStep = size(angleProfile(i).(s2u).b,3);
            

            % NaN-out timepoints from surrogates that were removed from
            % original theta.
            thetaS(repmat(isnan(theta),1,1,size(thetaS,3))) = NaN;

            % also apply angle reflection about 90 deg (if applied to
            % original theta). Important to do this BEFORE conversion
            % applied below, since the reflection assumes angles are in
            % oringinal coordiates.
            if bitReflectTheta
                thetaS = thetaS - 2*max(thetaS-90,0);
            end

            % likewise, we have to apply scale and offset conversion to the
            % surrogate thetas
            thetaS = angleScale * ...
                (repmat(angleProfile(i).(s2u).aBound(:,2),1,size(theta,2),size(thetaS,3))- thetaS);

            % instantiate
            angleProfile(i).(s2u).step_thetaMean = NaN(nTimes,nStep);
            angleProfile(i).(s2u).step_thetaSD = NaN(nTimes,nStep);
            angleProfile(i).(s2u).step_thetaSDist = NaN(nTimes,size(thetaS,3),nStep);
            angleProfile(i).(s2u).step_thetaSDistFitErr = NaN(nTimes,nStep);
            angleProfile(i).(s2u).step_thetaP = NaN(nTimes,nStep);
            angleProfile(i).(s2u).step_thetaSDist95 = NaN(nTimes,nStep);
            angleProfile(i).(s2u).step_thetaSDist99 = NaN(nTimes,nStep);
            angleProfile(i).(s2u).step_thetaSDist999 = NaN(nTimes,nStep);
            if strcmpi(ct,'cis') % ~isempty(regexp(s2u,'^main_step(_\d+)*$','once'))
                % gamma fits
                angleProfile(i).(s2u).step_thetaSDistFit = NaN(nTimes,2,nStep);
            else
                % normal fits. Turns out both are 2 parameter distributions
                % but if one selects different distributions, the
                % parameterization may change.
                angleProfile(i).(s2u).step_thetaSDistFit = NaN(nTimes,2,nStep);
            end
            % loop through steps (NOTE: this accommodates case of >1 step)
            for st = 1:nStep
                % loop through tSet
                for t = 1:nTimes
                    % extract step start and stop times
                    if strcmp('paramSearch',angleProfile(i).(s2u).fitMethod)
                        % for the paramSearch method, start and stop times should be
                        % well behaved)
                        start = angleProfile(i).(s2u).b(t,2,st);
                        stop = sum(angleProfile(i).(s2u).b(t,2:3,st),2);
                    else
                        % if step start or stop exceeds trial dur, set to min or max time.
                        % NOTE: this also has effect of including the first/last time point
                        % for any step that starts/stops just before/after the first/last
                        % time point.
                        start = max(min(tSet),angleProfile(i).(s2u).b(t,2),st);
                        stop = min(max(tSet),sum(angleProfile(i).(s2u).b(t,2:3)),st);
                        %         stop = min(max(tSet),angleProfile(i).(s2u).stopParam(t));
                    end

                    % NOTE that for some rows, there will not be a step
                    % function fit to the real data. This happens when all
                    % theta are NaN, which can occur if all angles are
                    % out-of-range. This occurs specifically for the "trans"
                    % fits (90 to 180 deg). When this occurs, start and stop
                    % times are NaN and consequently the mean angles under the
                    % step for each of the surrogate datasets is NaN.

                    % determine time range
                    bitT = tSet >= start & tSet <= stop;

                    % compute data mean and SD of real data
                    angleProfile(i).(s2u).step_thetaMean(t,st) = ...
                        nanmean(theta(t,bitT));
                    angleProfile(i).(s2u).step_thetaSD(t,st) = ...
                        nanstd(theta(t,bitT));

                    % extract distribution of mean theta in time range from surrogates
                    angleProfile(i).(s2u).step_thetaSDist(t,:,st) = ...
                        nanmean(thetaS(t,bitT,:),2);

                    % PDF fits and compute P value
                    if sum(~isnan(angleProfile(i).(s2u).step_thetaSDist(t,:,st))) > 1

                        if strcmpi(ct,'cis') 
                            % we use gamma fit when fitting to the folded (CIS)
                            % angles, which are gamma-distributed 0 to 90deg.
                            
                            % Using custom distribution fitting tool. Since
                            % for folded distribution we are interested in
                            % if the veridical angles are greater than
                            % chance, we a 1-sided p-value (i.e., 1 - CDF),
                            % NOTE that stats are performed prior to scale
                            % and offset, which cause angle size to get
                            % smaller with increasing similarity (the
                            % intuitive case). Prior to scale and offset,
                            % angle size increases with increasing
                            % A (which is why we take the upper
                            % tail here).
                            [angleProfile(i).(s2u).step_thetaP(t,st),...
                                angleProfile(i).(s2u).step_thetaSDistFit(t,:,st),...
                                angleProfile(i).(s2u).step_thetaSDistFitErr(t,st)] = ...
                                sigTest(angleProfile(i).(s2u).step_thetaMean(t,st),...
                                angleProfile(i).(s2u).step_thetaSDist(t,:,st)',...
                                'upper','gamma','errThresh',errThreshSigTestGamma_temp);
                            
                            % store fit type
                            if st==1 && t==1
                                angleProfile(i).(s2u).step_thetaSDistFitType = 'gamma';
                            end                        
                            
                            % Find percentiles based on fit
                            angleProfile(i).(s2u).step_thetaSDist95(t,st) = ...
                                gaminv(0.95,angleProfile(i).(s2u).step_thetaSDistFit(t,1,st),...
                                angleProfile(i).(s2u).step_thetaSDistFit(t,2,st));
                            % Find percentiles based on fit
                            angleProfile(i).(s2u).step_thetaSDist99(t,st) = ...
                                gaminv(0.99,angleProfile(i).(s2u).step_thetaSDistFit(t,1,st),...
                                angleProfile(i).(s2u).step_thetaSDistFit(t,2,st));
                            % Find percentiles based on fit
                            angleProfile(i).(s2u).step_thetaSDist999(t,st) = ...
                                gaminv(0.999,angleProfile(i).(s2u).step_thetaSDistFit(t,1,st),...
                                angleProfile(i).(s2u).step_thetaSDistFit(t,2,st));
                        else
                            % We use normal (Gaussian) fits when fitting to the
                            % unfolded (TRANS) angles, which (prior to scale
                            % and shift transformations) are centered at 0 deg
                            % and range [-90 to 90].

                            % Assuming normal distribution, fit probability
                            % distribution object:
                            pd_obj = fitdist(angleProfile(i).(s2u).step_thetaSDist(t,:,st)','normal');
                            
                            % store fit type
                            if st==1 && t==1
                                angleProfile(i).(s2u).step_thetaSDistFitType = 'normal';
                            end                        
                            
                            % store parameters from fit object
                            angleProfile(i).(s2u).step_thetaSDistFit(t,:,st) = pd_obj.ParameterValues;

                            % This analysis for focuses on SCALED and
                            % OFFSET angles > 90, i.e., angles representing
                            % a change in the sign of encoding. As such,
                            % prior to transformation these angles are in range [0 90], while the null distribution is supported on [-90 90]). Since we are
                            % interested in if the observed UNtransformed angles are greater than
                            % the null angles (i.e., changed MORE than
                            % chance), we take the p-value as 1 - CDF(), or
                            % the null distribution's upper tail:  
                            angleProfile(i).(s2u).step_thetaP(t,st) = 1 - cdf(pd_obj,...
                                angleProfile(i).(s2u).step_thetaMean(t,st));

                            % Find percentiles based on fit
                            angleProfile(i).(s2u).step_thetaSDist95(t,st) = ...
                                norminv(0.95, pd_obj.mean, pd_obj.sigma);
                            % Find percentiles based on fit
                            angleProfile(i).(s2u).step_thetaSDist99(t,st) = ...
                                norminv(0.99, pd_obj.mean, pd_obj.sigma);
                            % Find percentiles based on fit
                            angleProfile(i).(s2u).step_thetaSDist999(t,st) = ...
                                norminv(0.999, pd_obj.mean, pd_obj.sigma);
                            
                            clear pd_obj

                        end
                    end
                end
            end
            % compute mean, SD, Nth pctile across surrogates
            angleProfile(i).(s2u).step_thetaSMean = squeeze(nanmean(angleProfile(i).(s2u).step_thetaSDist,2));
            angleProfile(i).(s2u).step_thetaSMedian = squeeze(nanmedian(angleProfile(i).(s2u).step_thetaSDist,2));
            angleProfile(i).(s2u).step_thetaSSD = squeeze(nanstd(angleProfile(i).(s2u).step_thetaSDist,[],2));
            angleProfile(i).(s2u).step_thetaS95 = squeeze(prctile(angleProfile(i).(s2u).step_thetaSDist,95,2));
            angleProfile(i).(s2u).step_thetaS99 = squeeze(prctile(angleProfile(i).(s2u).step_thetaSDist,99,2));
            angleProfile(i).(s2u).step_thetaS999 = squeeze(prctile(angleProfile(i).(s2u).step_thetaSDist,99.9,2));

            % compute empirical p-value across surrogates (1-tailed). Here we
            % need to both repeat the observed theta mean for the number of
            % surrogates, and permute it so the nStep dimension is third and
            % the surrogate-repeated dimension is 2nd
            angleProfile(i).(s2u).step_thetaPe = squeeze(1-(sum(permute(repmat(angleProfile(i).(s2u).step_thetaMean,1,1,size(thetaS,3)),[1 3 2]) ...
                >= angleProfile(i).(s2u).step_thetaSDist,2) ./ sum(isfinite(angleProfile(i).(s2u).step_thetaSDist),2)));
            
            %%%
            % CLEAN UP
            clear foo thetaS theta tSet nTimes nStep

        end

        toc

        clear tdrTemp
    end
end

%% save summary structs 

% list of structs to save
structs_to_save = {'oTDRSummary','oTDRSummary_serial','oTDRSummary_PP',...
    'TDRSummary','TDRSummary_PP','TDRSummary_presRespPrevCond','BSSummary',...
    'BSSummary_PP','angleProfile'};

% path for saving results
dataPath_results = fullfile(dataPath,'results');

% ask if user wants to save structs
reply_save = input('Do you want to save the output structs? Y/N [Y]:','s');
if isempty(reply_save)
  reply_save = 'Y';
end
if strcmpi(reply_save,'Y')
    % ask to what path to save output
    reply_path = input(sprintf('Enter absolute path for saving [%s]:',dataPath_results),'s');
    if isempty(reply_path)
      reply_path = dataPath_results;
    end
    
    % check path exists
    if ~exist(reply_path,'dir')
        % attempt to make path
        foo = mkdir(reply_path);
        if foo==1
            fprintf('Created new dir for saving: %s\n',reply_path)
        else
            error('Could not create dir for saving at %s',reply_path);
        end
        clear foo
    end
    
    % save structs
    for fn = structs_to_save
        if exist(fn{:},'var')
            % save
            save(fullfile(reply_path,sprintf('results_monkey%s_%s',monkeyName(1),fn{:})),fn{:});
        else
            fprintf('Struct named "%s" was not found and will not be saved.\n',fn{:})
        end
    end    
end