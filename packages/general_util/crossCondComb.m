function [condComb] = crossCondComb(condSet)
%
% [condComb] = crossCondComb(condSet);
%
% Generates condComb with all combinations of the condition values passed
% in condSet.
%
% ACCEPTS
% condSet -- 1xN cell array where each cell corresponds to one of N
% conditions and contains a vector of condition values 
%
% RETURNS
% condComb -- MxN matrix of M unique combinations of the values across each
% of the N conditions.
%
% The code was suggested by John Cunningham, 23 Aug 2007. Written here by
% Daniel Kimmel, 08 Feb 2017

nParamToSplit = length(condSet);
nRow = 1;
for i = 1:nParamToSplit
    nRow = nRow * length(condSet{i});
end
sublength(nParamToSplit) = 1; % last column changes value at every row
for i = nParamToSplit-1:-1:1
    sublength(i) = sublength(i+1) * length(condSet{i+1});
end
for i = 1:nParamToSplit
    nCopies = nRow / sublength(i) / length(condSet{i});
    foo = repmat(condSet{i},sublength(i),1);
    foo = foo(:);
    condComb(:,i) = repmat(foo,nCopies,1);
end
condComb = unique(condComb,'rows');
