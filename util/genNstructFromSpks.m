function [N, NprevPres, NpresRespPrevCond, N_withLateReject, timeIntervals] = ...
    genNstructFromSpks(dataSpks, binWidth, samplesBefore, samplesAfter, ...
    alignment, samplesBeforePP, samplesAfterPP, lateRejectThresh)
% [N, NprevPres, NpresRespPrevCond, N_withLateReject, timeIntervals] = ...
%     genNstructFromSpks(dataSpks, binWidth, samplesBefore, samplesAfter, ...
%     alignment, samplesBeforePP, samplesAfterPP, lateRejectThresh)
%
% Function for generating "N structs" for cost-benefit task (see Kimmel et
% al., Nat Comms, 2020). The N structs contain the trial-level PSTHs for
% each neuron and serve as input to the task-agnostic getTrialAvgPSTHs(),
% which is called at the end of the present function and adds the
% trial-average responses to the N structs. Finally, the new N structs
% serve as input to getDataStruct(), which groups the data across neurons
% and must be called outside of the present function.
%
% Note that this function will compute trial-average PTSHs for ALL unique
% combinations of the task-relevant parameters (i.e., condition) found
% across ALL neurons, even if condition is not represented in a given
% neuron or if condition is not a "condition of interest". Conditions can
% be filtered in the subsequent call to getDataStruct().
%
% Accepts:
%   dataSpks -- Structure (1 x N), provided with public release of data,
%       where each element refers to one of N neurons and contains fields:
%           .spikeTime: Vector (S x 1) of absolute spike times (e.g.,
%               relative to beginning of experiment). Time units are
%               arbitrary but must be consistent across various
%               time-related inputs.
%           .trialInfo: struct of trial-level events, where each field
%               contains a (1 x R) vector of values for each of R trials.
%               Field names must include: bitAccept, nRwdOffered,
%               bitSingleton, fixationTime, offerTime,
%               rejectTime, rewardTime
%           .bitUnitPresent: logical vector (1 x R) of whether unit was
%               present (i.e., available for analysis) for each of R
%               trials.
%           .bitExcludeTrialFromPrevTrialAnal: logical vector (1 x R) of 
%               whether to exclude each of R trials from the previous-trial
%               epoch and analysis. 
%           .bitExcludeUnitFromPrevTrialAnal: logical scalar of whether to
%               include unit (i.e., ALL trials) from the previous-trial
%               epoch and analysis. 
%           .unitType: scalar indicating type of unit -- (1) single unit,
%               (2) multi-unit cluster, (3) background "hash". Note that
%               types (2) and (3) were combined and classified as
%               "multi-units" in Kimmel et al., 2020.
%
%   binWidth: scalar of duration of time bin to use for PTSHs.

%   samplesBefore: integer scalar of number of time bins to include in PTSH
%       BEFORE the event of interest (i.e., time=0) on the PRESENT trial.
%
%   samplesAfter: integer scalar of number of time bins to include in PTSH
%       AFTER the event of interest (i.e., time=0) on the PRESENT trial.
%
%   alignment: name of the task event to which to align PTSHs. Options
%   include:
%       "offer" -- align to the onset of the offer period on the present
%           trial.
%       "rejection" -- align to the time of the rejection on the present
%           trial ("accept" trials are ignored).
%       NOTE that alignment to fixation on the NEXT trial (used for the
%           previous-trial epoch and analysis) is done automatically and
%           in addition to the present-trial alignment when
%           alignment="offer". There is no previous-trial epoch when
%           aligning to the rejection.
%        
%   samplesBeforePP: integer scalar of number of time bins to include in
%       PTSH BEFORE fixation on the next trial (i.e., time=0) when
%       computing the previous-trial PSTHs.
%
%   samplesAfterPP: integer scalar of number of time bins to include in
%       PTSH AFTER fixation on the next trial (i.e., time=0) when
%       computing the previous-trial PSTHs.
%
%   lateRejectThresh: scalar specifying the time from the offer event after
%       which to classify rejections as "late". Generally use the median
%       rejection time across all trials. 
%
%
% Returns:
%   N: struct (1 x n) where each element refers to one of n neurons and
%       contains the following fields:
%           .deadNeuron: logical as to whether neuron was "dead" (i.e., not
%               consistently available for recording)
%           .psth: matrix (T * R) of average firing rate on each of R trials
%               and T time bins relative to the event of interest.
%           .all_times: vector (1 x T) of times of time bins in .psth (see
%               above)
%           .Benefit: vector (1 x R) of the offer size on each of R trials
%           .Choice: vector (1 x R) of the animals choice (accept = 1) on
%               each of R trials 
%           .OB: vector (1 x R) of whether the stimulus was the singleton
%               (=1), aka "oddball (OB)", or normal icons (=0) on each of R
%               trials.  
%           .LateReject: logical vector (1 x R) of whether the rejection was
%               classified as "late" (=1) based on input arg
%               lateRejectThresh (see above) for each of R trials. If offer
%               was rejected "early" or accepted, LateReject(r) = 0.
%           .unitType: scalar indicating type of unit -- (1) single unit,
%               (2) multi-unit cluster, (3) background "hash". Note that
%               types (2) and (3) were combined and classified as
%               "multi-units" in Kimmel et al., 2020.
%           .bitUnitPresent: logical vector (1 x R) as to whether the unit
%               was present for each of R trials.
%           .variableValues4Trial: matrix (R x P) of the values of each of
%               P task variables on each of R trials. A unique combination
%               of the variables defines a condition. Note that these task
%               variables need not be the same task variables that define
%               the sRAs in oTDR. In the case of the present task,
%               conditions are defined by the unique combination of
%               Benefit, Choice and OB (i.e., singleton stimulus) -- which
%               is the order of columns in variableValues4Trial -- whereas
%               (elsewhere) the sRAs are defined by Benefit, Choice, and
%               Expected Reward (i.e., Benefit * Choice).
%           .cond: struct (1 x C) where each element refers to each of C
%               experimental conditions (i.e., unique combination of task
%               variables across ALL units) with the following fields:
%               .psth: vector (T x 1) of trial-average firing rate for
%                   current condition in each of T time bins.
%               .numTrials: integer number of trials contributing to
%                   trial-average in .psth
%               .variableValues4Condition: vector (1 x P) of the values of
%                   the P task variables defining the present condition.
%
%   NprevPres: Similar to the N struct (above) except aligned to fixation
%       on the present trial and extending retrospectively into the
%       previous trial (according to input arg samplesBeforePP). Trials
%       are grouped according to the unique combination of Benefit, Choice
%       and OB (i.e., singleton stimulus) on the PREVIOUS trial. For
%       example, NprevPres(1).Benefit(1)=8 indicates that for Neuron 1, the
%       first row of the PSTH refers to a trial in which the previous trial
%       had an offer size of 1. Note that because certain trials are not
%       eligible for the previous trial analysis (e.g., because the
%       interval between trials was excessive), the number and indexing of
%       trials in NprevPres does not necessarily match those in N (i.e.,
%       present-trial analysis).
%   NpresRespPrevCond: Similar to N struct (above) except that trials are
%       grouped according to the unique combination of Benefit, Choice and
%       OB (i.e., singleton stimulus) on the PREVIOUS trial. (In this way,
%       the output is like NprevPres, but unlike NprevPres, the alignment
%       is to the offer on the present trial and the window extends
%       primarily forward in time, like the N struct.)  
%   N_withLateReject: Similar to N struct (above) except that the task
%       variables used to determine the conditions also includes
%       "LateReject", i.e., logical as to whether the rejection occured
%       after the time specified in the input arg, lateRejectThresh. As
%       such, the columns in .variableValues4Trial, refer to Benefit,
%       Choice, OB (i.e., singleton stimulus), and LateReject,
%       respectively.
%   timeIntervals: structure summarizing various time intervals between
%       task events. Serves as input to various downstream functions.

%% hard-coded params

% percentiles to use for time intervals
prctile2Use = [2.5 25 50 75 97.5];


%% compute single-trial PTSHs

numNeus = length(dataSpks);
if ~isempty(samplesBefore) && ~isempty(samplesAfter)
    bins = -binWidth*samplesBefore:binWidth:binWidth*samplesAfter;
else
    bins = [];
end
if ~isempty(samplesBeforePP) && ~isempty(samplesAfterPP)   
    binsPrevPres = -binWidth*samplesBeforePP:binWidth:binWidth*samplesAfterPP;
else
    binsPrevPres = [];
end


if isempty(lateRejectThresh)
    lateRejectThresh = Inf;
end

for n = 1:numNeus
     trialFlg = dataSpks(n).bitUnitPresent(:);
     if strcmpi(alignment,'offer')
        alignTimes = (dataSpks(n).trialInfo.offerTime);
     elseif strcmpi(alignment,'rejection')
        alignTimes = (dataSpks(n).trialInfo.rejectTime);
     else
        alignTimes = (dataSpks(n).trialInfo.fixationTime);
     end
     fixationTimes = (dataSpks(n).trialInfo.fixationTime);
     Benefit = dataSpks(n).trialInfo.nRwdOffered;
     Choice = dataSpks(n).trialInfo.bitAccept;
     OB = [dataSpks(n).trialInfo.bitSingleton];
 
     dataSpks(n).trialInfo.Benefit = Benefit(trialFlg);
     dataSpks(n).trialInfo.Choice = Choice(trialFlg);
     dataSpks(n).trialInfo.OB = OB(trialFlg);
     Times(n).alignTimes = alignTimes(trialFlg);
     
     % additional variable for whether rejection was late
     responseTime = dataSpks(n).trialInfo.rejectTime - dataSpks(n).trialInfo.offerTime;
     dataSpks(n).trialInfo.LateReject = responseTime(trialFlg) > lateRejectThresh; 
     
     % extract timing of events relative to aligntime
     ff = ['align2',upper(alignment(1)),lower(alignment(2:end))];
     timeIntervals.byNeuron(n).(ff).fix = dataSpks(n).trialInfo.fixationTime(trialFlg)...
         - Times(n).alignTimes;
     timeIntervals.byNeuron(n).(ff).offer = dataSpks(n).trialInfo.offerTime(trialFlg)...
         - Times(n).alignTimes;
     timeIntervals.byNeuron(n).(ff).reject = dataSpks(n).trialInfo.rejectTime(trialFlg)...
         - Times(n).alignTimes;
     timeIntervals.byNeuron(n).(ff).reward = dataSpks(n).trialInfo.rewardTime(trialFlg)...
         - Times(n).alignTimes;
     
     % Alignment of previous trial events. 
     % a) apply trialFlg to presPrevTrialFlg
     % b) apply presPrevTrialFlg to Benefit/choice/ER in same trial
     % reference 
     % c) reference previous benefit, choice, and OB for trial N to the r
     % espective variables for trial N-1
     presPrevTrialFlg = trialFlg' & ~(dataSpks(n).bitExcludeTrialFromPrevTrialAnal | dataSpks(n).bitExcludeUnitFromPrevTrialAnal);
     % also eliminate trials whose previous trial is invalid (assume the
     % first trial is always invalid since it does not have a previous
     % trial):
     presPrevTrialFlg = presPrevTrialFlg & [false, trialFlg(1:end-1)'];
     % store flags for later use:
     Flags(n).presPrevTrialFlg = presPrevTrialFlg;
     Flags(n).trialFlg = trialFlg';
     
     % get fixation times on current trial N:
     Times(n).fixationTimesNext = fixationTimes(presPrevTrialFlg);
     

     % use the absolute index of valid fixation trials, and simply
     % subtract 1 to find the previous trial:
     posPrev = find(presPrevTrialFlg) - 1;
     % check that all of posPrev is > 0. It should be, since
     % presPrevTrialFlg always eliminates the first trial:
     if any(posPrev <= 0)
         error('position of previous trial must always be > 0')
     end
     % Also check that all previous trials are valid. They should be, since
     % presPrevTrialFlg eliminates trials whose previous trial was invalid
     if ~all(trialFlg(posPrev))
         error('all previous trials must be valid. some were not')
     end
     
     %NOW WE CAN PROCEED EXTRACTING THE VARIABLES FROM THE PREVIOUS TRIALS: 
     dataSpks(n).prevPrestrialInfo.Benefit = Benefit(posPrev);
     dataSpks(n).prevPrestrialInfo.Choice = Choice(posPrev);
     dataSpks(n).prevPrestrialInfo.OB = OB(posPrev);
     
     % also extract timing of previous trial events relative to
     % fixationTimesNext: 
     timeIntervals.byNeuron(n).align2Fix.prevFix = dataSpks(n).trialInfo.fixationTime(posPrev)...
         - Times(n).fixationTimesNext;
     timeIntervals.byNeuron(n).align2Fix.prevOffer = dataSpks(n).trialInfo.offerTime(posPrev)...
         - Times(n).fixationTimesNext;
     timeIntervals.byNeuron(n).align2Fix.prevReject = dataSpks(n).trialInfo.rejectTime(posPrev)...
         - Times(n).fixationTimesNext;
     timeIntervals.byNeuron(n).align2Fix.prevRwd = dataSpks(n).trialInfo.rewardTime(posPrev)...
         - Times(n).fixationTimesNext;
     
     % eliminate any trial for whih an alignment time or predictor has a
     % NaN value 
     bitElim = isnan(Times(n).alignTimes) | isnan(dataSpks(n).trialInfo.Benefit) ...
         | isnan(dataSpks(n).trialInfo.Choice) | isnan(dataSpks(n).trialInfo.OB);
     if any(bitElim)
         
         Times(n).alignTimes(bitElim) = [];
         dataSpks(n).trialInfo.Benefit(bitElim) = [];
         dataSpks(n).trialInfo.Choice(bitElim) = [];
         dataSpks(n).trialInfo.OB(bitElim) = [];
         dataSpks(n).trialInfo.LateReject(bitElim) = [];
         
         % eliminate trials from interval measurements as well:
         for f = fieldnames(timeIntervals.byNeuron(n).(ff))'
            f = f{:};
            timeIntervals.byNeuron(n).(ff).(f)(bitElim) = [];
         end         
     end
     
     bitElimPP = isnan(Times(n).fixationTimesNext) | isnan(dataSpks(n).prevPrestrialInfo.Benefit) ...
         | isnan(dataSpks(n).prevPrestrialInfo.Choice) | isnan(dataSpks(n).prevPrestrialInfo.OB);
     if any(bitElimPP)
                  
         % we don't expect any trials to be eliminated here except when
         % aligning to rejection or subsampling
         if ~strcmpi(alignment,'rejection') 
             error('dud not expect trials to be eliminated here.')
         end
         Times(n).fixationTimesNext(bitElimPP) = [];
         dataSpks(n).prevPrestrialInfo.Benefit(bitElimPP) = [];
         dataSpks(n).prevPrestrialInfo.Choice(bitElimPP) = [];
         dataSpks(n).prevPrestrialInfo.OB(bitElimPP) = [];
         
         % eliminate trials from interval measurements as well:
         for f = fieldnames(timeIntervals.byNeuron(n).align2Fix)'
            f = f{:};
            timeIntervals.byNeuron(n).align2Fix.(f)(bitElimPP) = [];
         end
     end
     
     % store
     Flags(n).bitElim = bitElim;
     Flags(n).bitElimPP = bitElimPP;
     clear bitElim bitElimPP
end

%% compute some summary stats on the time intervals

% initialize summary matricies based on first neuron
n = 1;
for ff = fieldnames(timeIntervals.byNeuron(n))'
    ff = ff{:};
    for f = fieldnames(timeIntervals.byNeuron(n).(ff))'
        f = f{:};
        timeIntervals.(ff).(f).mean = NaN(numNeus,1);
        timeIntervals.(ff).(f).std = NaN(numNeus,1);
        timeIntervals.(ff).(f).std= NaN(numNeus,length(prctile2Use));
    end
end
    
% populate summary matrices
for n = 1:numNeus
    for ff = fieldnames(timeIntervals.byNeuron(n))'
         ff = ff{:};
         for f = fieldnames(timeIntervals.byNeuron(n).(ff))'
             f = f{:};
             timeIntervals.(ff).(f).mean(n,1) = nanmean(timeIntervals.byNeuron(n).(ff).(f));
             timeIntervals.(ff).(f).std(n,1) = nanstd(timeIntervals.byNeuron(n).(ff).(f));
             timeIntervals.(ff).(f).prctile(n,:) = prctile(timeIntervals.byNeuron(n).(ff).(f),prctile2Use);
         end
     end
end

% store percentiles used
timeIntervals.prctile2Use = prctile2Use;

%% creat N structs

N = [];
NprevPres = [];
NpresRespPrevCond = [];
uniqueConds = [];
uniqueConds_withLateReject = [];
for n = 1:numNeus

    if ~isempty(bins)
        [psth, N(n).deadNeuron] = getPSTHS(dataSpks(n).spikeTime, bins, Times(n).alignTimes);    
        N(n).psth = psth;
        N(n).all_times =bins(1:end-1)+binWidth./2;
        N(n).Benefit = dataSpks(n).trialInfo.Benefit;
        N(n).Choice = dataSpks(n).trialInfo.Choice;
        N(n).OB = dataSpks(n).trialInfo.OB;
        N(n).LateReject = dataSpks(n).trialInfo.LateReject;
        N(n).unitType = dataSpks(n).unitType;
        
        % because PSTH and condition flags have been re-indexed to valid
        % trials, we must limit the following field to valid trials only
        % (which makes it redundant, since by definition, the unit is
        % present on valid trials).
        % In addition, we must apply the "bitElim" flag, which removes NaN
        % trials that may result when aligning to rejection (thus
        % eliminating accept trials) or when subsampling (thus eliminating
        % random trials)
        N(n).bitUnitPresent = dataSpks(n).bitUnitPresent(Flags(n).trialFlg);
        N(n).bitUnitPresent = N(n).bitUnitPresent(~Flags(n).bitElim);
        
        % store matrix of variable values for each trial
        N(n).variableValues4Trial = [N(n).Benefit' N(n).Choice' N(n).OB'];
        
        % make copy of N struct...
        N_withLateReject(n) = N(n);
        % ... and add LateReject variable
        N_withLateReject(n).variableValues4Trial = [N_withLateReject(n).variableValues4Trial, N_withLateReject(n).LateReject'];
        
        % append conditions to master list 
        uniqueConds = [uniqueConds; N(n).variableValues4Trial];
        uniqueConds_withLateReject = [uniqueConds_withLateReject; N_withLateReject(n).variableValues4Trial];
                
        % the following is used to compute PSTHs of present trial data
        % sorted WRT previous offer conditions and are only used for
        % measuring certain projections and NOT used for discovering sRAs.
        % Therefore we can and should skip this part when not needed (e.g.,
        % aligning to rejections and subsampling). Including this part
        % would careful consideration of the bitElimPP flag.
        if ~strcmpi(alignment,'rejection') 
            
            % also generate PTSH of present trial responses with respect to
            % previous trial conditions.
            NpresRespPrevCond(n).deadNeuron = N(n).deadNeuron;
            % PSTH is same as above, but limited to valid previous trials.
            % Because Flags(n).presPrevTrialFlg is indexed to all trials, but
            % PSTH is index to valid present trials, we have to re-index
            % Flags(n).presPrevTrialFlg to valid present trials:
            NpresRespPrevCond(n).psth = psth(:,Flags(n).presPrevTrialFlg(Flags(n).trialFlg));
            NpresRespPrevCond(n).all_times = N(n).all_times;
            % following fields are based on previous trial conditions
            NpresRespPrevCond(n).Benefit = dataSpks(n).prevPrestrialInfo.Benefit;
            NpresRespPrevCond(n).Choice = dataSpks(n).prevPrestrialInfo.Choice;
            NpresRespPrevCond(n).OB = dataSpks(n).prevPrestrialInfo.OB;
            % following fields do not depend on previous trial
            NpresRespPrevCond(n).unitType = dataSpks(n).unitType;
            % limit to valid previous trials
            NpresRespPrevCond(n).bitUnitPresent = dataSpks(n).bitUnitPresent(Flags(n).presPrevTrialFlg);
            
            % store matrix of variable values for each trial. Note that for
            % the NpresRespPrevCond struct, the fields "Benefit", "Choice",
            % and "OB" refer to events in the PREVIOUS trial.
            NpresRespPrevCond(n).variableValues4Trial = ...
                [NpresRespPrevCond(n).Benefit' NpresRespPrevCond(n).Choice' NpresRespPrevCond(n).OB'];
        end
    end 
    
    % The following computes PSTHs aligned to fixation on the PRESENT trial
    % and looks retrospectively for responses on the PREVIOUS trial sorted
    % WRT the task variables on the previous trial.
    if ~isempty(binsPrevPres)
        [prevPresPsth, NprevPres(n).deadNeuron] = getPSTHS(dataSpks(n).spikeTime, binsPrevPres, Times(n).fixationTimesNext);
        NprevPres(n).psth = prevPresPsth;
        NprevPres(n).all_times = binsPrevPres(1:end-1)+binWidth./2;
        NprevPres(n).Benefit = dataSpks(n).prevPrestrialInfo.Benefit;
        NprevPres(n).Choice = dataSpks(n).prevPrestrialInfo.Choice;
        NprevPres(n).OB = dataSpks(n).prevPrestrialInfo.OB;
        NprevPres(n).unitType = dataSpks(n).unitType;
        % limit to valid previous trials
        NprevPres(n).bitUnitPresent = dataSpks(n).bitUnitPresent(Flags(n).presPrevTrialFlg);
        
        % store matrix of variable values for each trial. Note that for
        % the NprevPres struct, the fields "Benefit", "Choice",
        % and "OB" refer to events in the PREVIOUS trial.
        NprevPres(n).variableValues4Trial = ...
            [NprevPres(n).Benefit' NprevPres(n).Choice' NprevPres(n).OB'];
        
        % We can skip this step if present offer data was computed
        if isempty(bins)
            % append conditions to master list 
            uniqueConds = [uniqueConds; NprevPres(n).variableValues4Trial];
        end            
    end
    
end

%% Compute trial-averages for each condition within each neuron

% find unique combinations of variables from master list
uniqueConds = unique(uniqueConds,'rows');
% sort according to OB (normal offers first), Choice, enefit
uniqueConds = sortrows(uniqueConds,[3 2 1]);

% find unique combinations of variables including LateReject conditions 
uniqueConds_withLateReject = unique(uniqueConds_withLateReject,'rows');
% sort according to OB, Choice, LateReject, Benefit
uniqueConds_withLateReject = sortrows(uniqueConds_withLateReject,[3 2 4 1]);

if ~isempty(bins)
    % add trial-average responses to N struct
    N = getTrialAvgPSTHs(N, uniqueConds);
    % add trial-average responses to NpresRespPrevCond struct
    NpresRespPrevCond = getTrialAvgPSTHs(NpresRespPrevCond, uniqueConds);
    % add trial-average responses to N_withLateReject struct
    N_withLateReject = getTrialAvgPSTHs(N_withLateReject, uniqueConds_withLateReject);
end
if ~isempty(binsPrevPres)
    % add trial-average responses to NprevPres struct
    NprevPres = getTrialAvgPSTHs(NprevPres, uniqueConds);
end

