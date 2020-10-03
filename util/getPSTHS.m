function [psth, deadNeuronFlg] = getPSTHS(spkTimes, bins, eventTimes, binWidth)
% [psth, deadNeuronFlg] = getPSTHS(spkTimes, bins, eventTimes, [binWidth])
%
% Compute PSTH from spike times.
%
% Accepts:
% spkTimes: vector of absolute spike times relative to arbitrary time
%   point that must match eventTimes, with time units that must match
%   those in bins and eventTimes (see below).
% bins: vector (1 * t) of time bin left (i.e., leading) edges, within which
%   spikes will be pooled. Time units must match those in spkTimes and
%   eventTimes. Note that spikes occuring exactly at or after the last bin
%   edge are not counted.
% eventTimes: vector (1 * r) of times (relative to arbitrary time
%   point that must match spkTimes) of r events to which firing rate will
%   be aligned. 
% binWidth: (Optional) Integer or 1 * t vector specifying duration of bins.
%   If not provided, will compute duration from bins. If duration is not
%   uniform across bins, error will be thrown. By providing value of 1,
%   will cause function to return PSTH of spike counts instead of firing
%   rate.
%
% Returns:
% psth: matrix (t * (r-1)) of firing rate (but see binWidth, above, for
%   spike counts) in time bin t (rows) relative to event r (colums). Note
%   that last bin edge (i.e., bins(end)) sets upper bound on penultimate
%   bin (i.e., bins(end-1)), but spikes occuring exactly at or after the
%   last bin edge are not counted.
% deadNeuronFlg (Bool): Whether firing rate across all elements of psth has
%   zero variance (and thus assumed that neuron is "dead").

% optional binwidth input
if ~exist('binWidth','var')
    % compute unique durations of bins (allowing for some numerical
    % tolerance relative to system eps)
    binWidth = unique(round(diff(bins),round(-log10(eps))-3));
    if length(binWidth) > 1
        error('If binWidth not provided, bin durations must be uniform.')
    end
else
    % get binwidth in correct orientation (column vector)
    binWidth = binWidth(:);
end

% initialize
deadNeuronFlg = false;
numTrials = length(eventTimes);
psth = nan(length(bins)-1, numTrials);

% loop through events
parfor c = 1:numTrials
    % compute spike times relative to current event
    spkTimes_c = spkTimes-eventTimes(c);
    % compute firing rate (spike count / binWidth) for current event
    PSTH = histc(spkTimes_c,bins)./binWidth;
    % store PSTH for event, eliminating last bin
    psth(:,c) = PSTH(1:end-1);
end

% flag for zero-variance PTSH (i.e., "dead neuron")
if var(psth(:))==0
    deadNeuronFlg = true;
end
