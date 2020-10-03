function [vect] = TDRAlignDataAndRegressTimes(mat,RAN,varargin);
% [vect] = TDRAlignDataAndRegressTimes(mat,RAN,[name,value]);
%
% Helper function to extract data for regression axis RAN from matrix MAT
% as contained in the TDR.varanalysis struct (from TDR.m), where MAT is of
% the form: data time bins x regression time bins x regression axis [x
% signal], where the 4th dimension is optional.
%
% MAT contains variance data based on data from time t as projected onto a
% regression axis RAN from regression time rt, and in some cases broken
% down by relevant signal variance relative to signal SN. In many cases, we
% wish to extract a single vector VECT of variance data where t = rT, RAN =
% RAN, and SN = RAN. Because the timebines in rT can be a multiple of the
% timebins in t, we cannot take a simple diagonal across dimensions 1 and 2
% in MAT. Instead, we have to align positions in dimensions 1 and 2
% according to their respective time vectors.
% 
% INPUTS as name value, pairs:
%
% Optionally provide either T AND RT (preferred), OR nRegressBin. When T
% and RT are provided, we check that the absolute times align. When
% nBinRegress is provided, we assume that T and RT start at the same time.
% Note that when T and RT are provided, we compute nBinRegress based on the
% ratio of bin sizes between the two vectors. You can instead provide
% nBinRegress in addition to T and RT, in which case the provided value
% will be used. If none of the above are provided, we compute nBinRegress
% based on the ratio of sizes of dimensions 2 to 1 in MAT.
%
% T -- 1 x t vector of data times (corresponds to dimension 1 in MAT)
%
% RT -- 1 x rt vector of regression times (corresponds to dimension 2 in
%       MAT)
%
% nBinRegress -- scalar specifying the ratio of data timebins to regression
%       timebins.
%
% SN -- scalar specifying the signal (corresponds to dimension 4 in MAT)
%       to extract. Defaults to RAN.
%
% RETURNS
% VECT -- t x 1 vector of variance data
%
% Daniel Kimmel, 10 Dec 2016

%% Default param values
T = [];
RT = [];
nBinRegress = [];
SN = RAN;

%% collect optionally provided parameters
warnopts(assignopts(who, varargin));

%% analyze time vectors

dT = mean(round(diff(T),3));
dRT = mean(round(diff(RT),3));

if isempty(nBinRegress)
    if ~isempty(T) && ~isempty(RT)
        % compute based on ratio of RT binwidth to T binwidth
        nBinRegress = dRT / dT;
    else
        nBinRegress = size(mat,1) ./ size(mat,2);
    end
end
if nBinRegress ~= round(nBinRegress)
    error('Not integer number of time bins per regress bins')
end

if ~isempty(T) && ~isempty(RT)
    % find T and RT alignment explicitly
    
    % for now, just check that T and RT start in the same place. If they
    % don't, throw an error. We'll worry about custom alignment when it's a
    % problem
    if (T(1) - dT/2) - (RT(1) - dRT/2) > 10*eps
        error('T and RT vector do not start at the same time')
    else
%         % determine bin start and stop
%         Tedge = [T(:)-dT/2, T(:)+dT/2];
%         RTedge = [RT(:)-dRT/2, RT(:)+dRT/2];
% 
%         % index of Tedge where starting edge matches starting edge in RT
%         ismember(T(:,1),RT(:,1));
        
    end
end
            
% initialize
vect = NaN(size(mat,1),1);

% specify position in 4th dimension. When there are 3 dims, set to 1
if ndims(mat) <= 3
    % for VE
    SN = 1;
elseif ndims(mat) > 4 || ndims(mat) < 2
    error('MAT must have 2, 3 or 4 dimensions');
end

% loop through regress bins to extract time bins
for i = 1:size(mat,2)
    vect(i*nBinRegress-(nBinRegress-1):i*nBinRegress) = ...
        mat(i*nBinRegress-(nBinRegress-1):i*nBinRegress,i,RAN,SN);
end
