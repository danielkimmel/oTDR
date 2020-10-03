%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <subspaces>
% Copyright (C) 2016 Gamaleldin F. Elsayed and John P. Cunningham (see full notice in README)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ix = align_ix(w1, w2, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function evaluates the alignment index between two datasets
% occupying two subspaces.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%       - w: are the orthonormal basis of size (neurons x dimensionality2) 
%         that defines the subspace that we want to test if response 1 is
%         aligned to. 
%       - D: is the matrix of size (neurons x observations) that contains 
%         the responses. 
% Outputs:
%       - Ix: is the resulting alignment index of response1 to subspace 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ix = align_ix(w, D)
    if iscell(w)
        w = w{1};
    end
    if iscell(D)
        D = D{1};
    end
    w = orth(w);                                                           % make sure the basis are orthogonal
    dim = size(w, 2);                                                      % dimensionality of the subspace
    C = (D*D');                                                            % calculate covariance (note, scale of covariance is not important as it cancels out in the alignment index formula)
    s = svd(C);                                                            % singular values of the covariance
    Ix = trace(w'*C*w)./sum(s(1:dim));                                     % calculate the alignment index
end