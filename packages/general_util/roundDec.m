% function newnum = roundDec(number,ndecimal)
% rounds "number" to nearest decimal place as secified by "ndecimal".
% works only in base 10.
% ndecimal must be an integer. If ndecimal is not provided, it is relaced
% by floor(abs(log10(eps)))


function newnum = roundDec(number,ndecimal)

if nargin < 2
    ndecimal = floor(abs(log10(eps)));
end

if round(ndecimal) ~= ndecimal
    error('ndecimal must be an integer');
end

if ndecimal > 0
    newnum = round(number.*10.^ndecimal) ./ 10^ndecimal;
elseif ndecimal == 0
    newnum = round(number);
else
    error('ndecimal cannot be a negative number');
end

