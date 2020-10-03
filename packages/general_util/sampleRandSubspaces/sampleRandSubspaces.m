%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <subspaces>
% Copyright (C) 2016 Gamaleldin F. Elsayed and John P. Cunningham (see full notice in README)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Ix] = sampleRandSubspaces(dim, useCov, IndexFn, numSamples)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function evaluates the alignment index between two datasets
% occupying two subspaces.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%       - dim: dimensionality of the subspace.
%       - useCov: sample subspaces based on some data covariance. This
%       ensures that the subspaces are sampled from the same data space.
%       This data alignment option can be disabled by setting useCove to
%       identity matrix.
%       - IndexFn: contains function name that calculates an index between
%       two samples spaces.
%       - numSamples: number of samples. 
%       - varargin: other optional inputs. These can be used for other
%       inputs for the IndexFn.
% Outputs:
%       - Ix: this is the result distribution of indices based on random
%       subspaces.
%       - Qout: 1 x S cell array for S samples where each cell contains the
%       random vector(s) generated in the subspace of the dataset that are
%       passed to IndexFn. Added by Daniel Kimmel, 23 Jan 2017.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ix,Qout] = sampleRandSubspaces(dim, useCov, IndexFn, numSamples, varargin)
%% Project Neural Activity on random directions for significance analysis 
%%% generate random unit basis (random unit vectors in neural space)

rng shuffle
N = size(useCov, 1);                                                       % space size
[U, S, ~] = svd(useCov);                                                   
% NOTE: this step may be needlessly complex. Rather than linearize S then
% re-diagonalize it, why not just take the sqrt of the diagonal matrix S,
% i.e., diag(sqrt(diag(S))) = sqrt(S)
biasMtx = U*diag(sqrt(diag(S)));                                           % evaulates the biasing matrix that bias the random vectors to data space
numSubspaces = length(dim);
parfor s = 1:numSamples
    %% sample random subspace
    Q = cell({});
    for j = 1:numSubspaces
        Q{j} = randn(N, dim(j));                                                    % initially sample uniformly
        [~, Q{j}] = normVects(Q{j}); 
        Q{j} = orth(biasMtx*Q{j});                                                 % bias the sampled vectors to data space
    end
    % Store random vector/matrix
    Qout{s} = Q;
    %% evaluate statistic measure between space 1 and 2
    if ~isempty(IndexFn)
%         if isempty(varargin)
%             Ix{s} = feval(IndexFn, Q);
%         else
            Ix{s} = feval(IndexFn, Q, varargin{:});
%         end
    end
end

% create for output
if isempty(IndexFn)
    Ix = [];
end

end

