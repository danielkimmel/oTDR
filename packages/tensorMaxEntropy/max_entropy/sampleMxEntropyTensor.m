%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <tensorMaxEntropy>
% Copyright (C) 2016 Gamaleldin F. Elsayed and John P. Cunningham 
%       (see full notice in README)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x, Summary] = sampleMxEntropyTensor(covConst, numSamples, bigCovEigValues)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function draws tensor samples from the maximum entropy 
% distribution with marginal covariance constraints 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%       - covConst: is a cell where each element contains a specified
%       marginal covariance matrix.
%       - numSamples: number of draws to return
%       - bigCovEigValues: (optional) this is an optional input. If the
%       solution to the maximum entropy distripution is known; provide the
%       eigenvalues solution in the vector bigCovEigValues.
% Outputs:
%       - x: the resulting tensor samples from the maximum entropy
%       distribution
%       - Summary:
%               .optSummary: optimization results:
%               .bigCovEigValues: eigenvalues solution to the maximum
%                entropy problem.
%               .eigValues: cell with each element containing the vector of
%                the eigenvalues of the specified marginal covariances.
%               .eigVectors: cell with each element containing the matrix 
%                of the eigenvectors of the specified marginal covariances.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, Summary] = sampleMxEntropyTensor(covConst, numSamples, bigCovEigValues)
x = nan;
tensorSize = length(covConst);                                             % tensor size; ie the number of different dimensions of the tensor
dim = nan(tensorSize, 1);                                                  % tensor dimensions
eigVectors = cell(tensorSize, 1);                                          % eigenVectors of each of the specified marginal covariances
eigValues = cell(tensorSize, 1);                                           % eigenValues of each of the specified marginal covariances
trSigma = nan(tensorSize, 1);                                              % sum of each of the eigenValues of each of the specified marginal covariances
for i = 1:tensorSize                                                       % load all the inputs
    Sigma = covConst{i};
    dim(i) = size(Sigma,1);
    [Q, S] = svd(Sigma);
    [S, ix] = sort(diag(S), 'descend');
    Q = Q(:, ix);
    eigVectors{i} = Q;
    eigValues{i} = S;
    trSigma(i) = trace(Sigma);
end
Summary.eigValues = eigValues;
Summary.eigVectors = eigVectors;
%% the marginal covariances should all have the same trace (i.e. the sum of their eigenvalues should be equal)
if ~(sum((trSigma-mean(trSigma))>=-(sum(dim)*sqrt(eps)) & (trSigma-mean(trSigma))<=(sum(dim)*sqrt(eps)))==length(trSigma))
    display('Error: the covariance matrices should have exactly the same trace')
    return                                                                 % abort if the marginal covariance specified are not correct
end

%% find the solution to the maximum entropy distribution optimization problem
if ~exist('bigCovEigValues','var')  
    [bigCovEigValues, optSummary] = optMaxEntropyTensor(eigValues);        % This function performs the optimization for the lagrangians of the max entropy problem and provides the corresponding eigenvalues of the max entropy sampling distribution 
    Summary.optSummary = optSummary;
    Summary.bigCovEigValues = bigCovEigValues;
end
%% sample from maximum entropy distribution
x = randn(prod(dim),numSamples);                                           % draw random samples from a normal distribution
x = bsxfun(@times, bigCovEigValues.^0.5, x);                               % multiply the samples by the eigenvalues of the covariance matrix the maximum entropy distribution
Qs = cell(tensorSize, 1);                                                  % load the eigenvectors of the covariance matrix of the maximum entropy distribution           
for i = 1:tensorSize
    Qs{i} = eigVectors{i};
end
x = real(kron_mvprod(Qs, x));                                              % efficiently multiply the samples by the eigenvectors of the covariance matrix of the maximum entropy distribution           
end