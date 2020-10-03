% function [m,n] = idealSubplotDim(nAxes,aspectRatioIdeal,[bitAllowEmptyAxes])
%
% Returns the ideal m x n dimensions of a figure subplot based on the
% number of axes needed and the dimensions of the display device. 
% ARGUMENTS:
% 1. nAxes = number of axes required for subplot
% 2. aspectRatioIdeal = aspectRatio (height/width) of display device (e.g., a
% flat panel 19" LCD monitor is approx. 12/15.
%
% OPTIONAL:
% 3. bitAllowEmptyAxes = logical on whether to allow for empty axes
% positions (e.g., 3 axes arranged on 2 x 2 grid), which will allow axes to
% most closely match aspectRatioIdeal. [Default = TRUE]

% OUTPUT:
% 1. m and n are the ideal number of rows and columns, respectively.
%
% Daniel Kimmel, originally written ~2005, updated with bitAllowEmptyAxes
% 13 Nov 2017

function [m,n] = idealSubplotDim(nAxes,aspectRatioIdeal,varargin)

if nargin <= 2 || isempty(varargin{1})
    bitAllowEmptyAxes = true;
else
    bitAllowEmptyAxes = varargin{1};
end

if bitAllowEmptyAxes
    % compute number of rows such that aspect ratio of axes most closely
    % matches aspectRatioIdeal
    mFloor = floor(sqrt(nAxes * aspectRatioIdeal));
    mCeil = ceil(sqrt(nAxes * aspectRatioIdeal));
%     m = round(sqrt(nAxes * aspectRatioIdeal));

    % tentative number of columns based on larger number of rows
    nCeil = ceil(nAxes / mCeil); 
    % choose m to protect against empty rows
    if nCeil * (mCeil - 1) == nAxes 
        m = mFloor;
    else
        m = mCeil;
    end
    
%     % choose m to minimize empty axes, with bias toward fewer rows
%     if ceil(nAxes / mFloor) * mFloor <= ceil(nAxes / mCeil) * mCeil
%         m = mFloor;
%     else
%         m = mCeil;
%     end
    
    n = ceil(nAxes / m);
else
    n = 1;
    m=1;
    % introduce extra axes for prime number of files
    if isprime(nAxes) & nAxes~=2
        nAxes = nAxes+1;
    end
    
    for i=1:nAxes/2
        if round(nAxes/i) == nAxes/i
            n(i) = nAxes/i;
            m(i) = nAxes/n(i);
        end
    end
    warning off MATLAB:divideByZero;
    aspectRatioDiff = abs(aspectRatioIdeal - m./n);
    warning on MATLAB:divideByZero;
    [foo,i] = min(aspectRatioDiff);
    n = n(i);
    m=m(i);
end