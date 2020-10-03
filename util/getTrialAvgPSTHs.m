function N = getTrialAvgPSTHs(N, uniqueConds, resampleFlg)
% N = getTrialAvgPSTHs(N, uniqueConds, [resampleFlg])
%
% Computes trial-average PSTHs from single-trial PTSHs.
% Args:
%   N: 1 x n struct where each element refers to individual neuron with fields:
%       .psth: matrix (t * r) of firing rate of neuron n on trial r
%           (columns) at time t (rows) relative to event-of-interest. 
%       .variableValues4Trial: matrix (r * p) of the values of the p
%           task variables (columns) for trial r (rows).
%       (Additional fields will be ignored)
%
%   uniqueConds: matrix (c * p) of the values of the p task variables
%       (columns) for the c condition (rows) for which you wish to
%       compute the trial-average. Number of columns p and the specific
%       values for each variable must match those provided in
%       N.variableValues4Trial (see above).
%
%   resampleFlg (Bool): (Optional) Whether to sample trial-level data with
%       replacement when computing trial-average for each condition.
%       Default = False.
%
% Returns:
%   N struct with additional substruct (1 * c) .cond, where each element
%       refers to one condition provided in uniqueConds (even if condition
%       was not present for current neuron), with the following fields:
%           .psth: vector (1 * t) of trial-average firing rate at time t
%               relative to event-of-interest.
%           .numTrials: int of number of trials contributing to
%               trial-average.
%           .variableValues4Condition: vector (1 * p) of the values of
%               the p task variables for the present condition.

% default to no resampling if not specified
if ~exist('resampleFlg','var')
   resampleFlg = false;
end

% loop through n neurons
for n = 1:length(N)    
    % extract for convenience
    conds = N(n).variableValues4Trial;
    
    % check that number of variables matches that provided in uniqueConds
    if size(conds,2) ~= size(uniqueConds,2)
        error('Number of variables (columns), %d, specified in N().variableValues4Trial for neuron %d does not equal number of variables provided in uniqueConds, %d. Consider removing irrelevant variables.',...
            size(conds,2), n, size(uniqueConds,2));
    end
    
    % loop through unique conditions for trial-average
    for j = 1:size(uniqueConds,1)
        
        % generate bitmask of those trial matching present condition
        condsFlg = ismember(conds,uniqueConds(j,:),'rows');
        
        % extract PTSH for trials matching present condition
        trialsPSTH = N(n).psth(:,condsFlg);
        
        if resampleFlg && ~isempty(trialsPSTH)
           % if requested, resample trials with replacement for present condition 
           trialsPSTH = trialsPSTH(:, randi(sum(condsFlg), sum(condsFlg),1)); 
        end
        
        % compute and store trial-average 
        N(n).cond(j).psth =  mean(trialsPSTH,2);
        % store number of trials contributing to average
        N(n).cond(j).numTrials =  sum(condsFlg);
        % store values of the task variables for present condition
        N(n).cond(j).variableValues4Condition = uniqueConds(j,:);
        
    end
    
end
