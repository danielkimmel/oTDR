% [dataset_times] = getTimes4oTDR(params,[tMinMax],[magThresh])
%
% Function to specify temporal epoch(s) from which to derive dataset on
% which oTDR is run. 
%
% INPUTS
% params -- "Summary" struct output of TDR(). If struct contains optional
%   field .dataset_times, this field will be returned unalertered,
%   effectively overriding all other processing. See below for format of
%   .dataset_times
%
% Optional inputs:
% tMinMax -- E x 2 matrix of minimum start time (Col 1) and maximum end
%       time (Col 2) for inclusion of data in a given dataset, with each
%       dataset specified in a given row for a total of E datasets. When
%       magThresh is provided (below), it overrides tMinMax.
% magThresh -- scalar specifying threshold of dRA magnitude as returned
%       by TDR() and passed as input params.normdRAs. When dRA magnitude
%       exceeds magThresh, time bin will be included for time bins > 0
%       (time bins < 0 are excluded). Note that dRA magnitude is first
%       normalized by the peak magnitude for that dRA across all times > 0,
%       thus magThresh is supported on [0, 1]. When magThresh is provided,
%       a dataset will be defined for each dRA (i.e., columns of
%       params.normdRAs), including the dRA for the constant term in the
%       regression (which is otherwise generally excluded from most dRA
%       analyses). When magThresh is empty or not provided, times in
%       tMinMax are used.
%       
% RETURNS
% dataset_times -- 1 x E cell array, where each cell specifies the time bins
%   to include for each of the E datasets. Time bins are specified as a
%   vector and reference the time bin centers provided in
%   params.Times.all_times. 

function [dataset_times] = getTimes4oTDR(params,varargin)
    % Collection optional input params
    if nargin > 1
        tMinMax = varargin{1};
    end
    if nargin <= 1 || isempty(tMinMax)        
        tMinMax = [];                
    end     
    
    if nargin > 2
        magThresh = varargin{2};
    else
        magThresh = [];
    end 
    
    if isfield(params, 'dataset_times')
        dataset_times = params.dataset_times;
    elseif ~isempty(magThresh)
        
        % times need to be based on regression times
        regressTimes = params.Times.regressTimes;
        
        %  mask out negative times.
        mskTime = regressTimes > 0;        
        regressTimes = regressTimes(mskTime); % limit to time during trial

        % store some useful values
        all_times = params.Times.all_times;
        binwidth = mean(diff(all_times)); 
        regressBins = params.Times.regressBins; % number of original time 
                                                % bins grouped together for
                                                % purposes of regression
        
        % loop through each dRA
        for i = 1:size(params.normdRAs,2)
            
            % extract norm of dRA (i.e., beta vector magnitude) for
            % non-negative times
            mag = params.normdRAs(mskTime,i);

            % normalize by max value:
            mag = mag / max(mag);

            % find times where mag > magThresh
            t = regressTimes(mag > magThresh);

            % convert back from regression times to all times
            t = round(bsxfun(@plus, repmat(t(:), 1, regressBins), (-((regressBins-1)/2):1:((regressBins-1)/2))*binwidth),2);
            dataset_times{i} = sort(t(:));
        end
        
    elseif ~isempty(tMinMax)
        % use all_times
        all_times = params.Times.all_times;
        
        % loop through each epoch
        for i = 1:size(tMinMax,1)

            % find times meeting criteria. Round to 2 decimal places
            dataset_times{i} = round(all_times(all_times >= tMinMax(i,1) & all_times <= tMinMax(i,2)),2);
        end
        
    else
        error('If dataset_times is not provided in PARAMS input, then either tMinMax or magThresh must be provided')
    end
end