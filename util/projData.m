function proj = projData(dataTensor,RA,varargin)
% proj = projData(dataTensor,RA,[bitMeanSubtract])
% 
% function for computing projections
% takes dataTensor -- times by neurons by conditions
% RA -- regression axes -- neurons by axes
% [bitMeanSubtraction] -- optional flag to subtract mean from data prior to
%       projection (default = true)
%
% returns PROJ -- times by conditions by axes
%
% Daniel Kimmel, 2017 January 15

if length(varargin) > 0
    bitMeanSubtraction = varargin{1};
else
    bitMeanSubtraction = false;
end

% data tensor is times by neurons by conditions
[T,N,C] = size(dataTensor);

% bigA is times*conditions by neurons
bigA = reshape(permute(dataTensor, [1 3 2]), T*C, N);

% mean center bigA:
if bitMeanSubtraction
    bigA = bsxfun(@minus,bigA,mean(bigA,1));
end

% normalize axes. RA is neurons x axes
[~,RA] = normVects(RA);

% project (PROJ is TC x axes)
proj = bigA * RA;

% reshape project into T x C x axes
proj = reshape(proj,T,C,size(RA,2));

end