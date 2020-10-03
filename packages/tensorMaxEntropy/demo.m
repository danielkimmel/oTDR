%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <tensorMaxEntropy>
% Copyright (C) 2016 Gamaleldin F. Elsayed and John P. Cunningham 
%       (see full notice in README)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a demonstration of how to use this code package 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% startup
tensorSize = 3; % 3D-tensor size example
% tensorSize = 5; % 5D-tensor size example
dim = randi(20, 1, tensorSize); % tensor dimensions
%% genarate marginal covariance matrices of size specified in dim
covConst = cell(tensorSize, 1);
for i = 1:tensorSize
    eigVectors = orth(randn(dim(i)));
    eigValues = abs(randn(dim(i),1));
    % the covariances should be positive semidef with the same trace
    Sigma = eigVectors*diag(eigValues)*eigVectors';
    covConst{i} = Sigma./trace(Sigma);
end
%% sample from the maximum entropy distribution
numSamples = 100;
[x] = sampleMxEntropyTensor(covConst, numSamples);

%% check if covariance of sampled data match the covariance constraints
estSigma = cell(tensorSize, 1);
for s = 1:numSamples
   Xs = reshape(x(:,s), dim);
   
   for i = 1:tensorSize
       if s == 1
           estSigma{i} = zeros(dim(i));
       end
       jSet = setdiff(1:tensorSize, i);
       Xi = reshape(permute(Xs, [i jSet]), dim(i), prod(dim(jSet)));
       estSigma{i} = estSigma{i}+(Xi*Xi')./numSamples;
   end
end

Er = nan(tensorSize,1);
for i = 1:tensorSize
    Er(i) = norm(mean(estSigma{i},3)-covConst{i},'fro')^2./norm(covConst{i},'fro')^2*100;
    fprintf('\n Error in Sigma_%d: %.2f %% \n', i, Er(i))
end