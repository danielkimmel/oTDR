%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <tensorMaxEntropy>
% Copyright (C) 2016 Gamaleldin F. Elsayed and John P. Cunningham 
%       (see full notice in README)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [D, Summary] = optMaxEntropyTensor(eigValues)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves for the eigenvalues of the covariance matrix of the
% maximum entropy distribution with specified eigenvalues of the marginal 
% covariance matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%       - eigValues: cell with each element containing the vector of
%         eigenvalues of the specified marginal covariance matrices.
% Outputs:
%       - D: is the vector of the eigenvalues of the covariance of the
%       maximum entropy distribution
%       - Summary:
%               .Lagrangians: are the eigenvalues of the largrangian
%                multipliers of the optimization program
%               .objCost: is the final original objective cost.
%               .logObjperIter: is the final log objective cost.
%               .logObjperIter: is the log objective cost at each
%                iteration. Note, optimization is performed only on the
%                log objective because the original objective can be
%                problematic and slowly converges when the specified
%                marginal covariance matrices are low rank. The optimal
%                solution of the log objective and the original objective 
%                is the same and both the log objective and original objective
%                values at the global optimum should be 0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D, Summary] = optMaxEntropyTensor(eigValues)
 % if the marginal covariances are low rank then the number of variables 
 % that we solve for are less. If full rank the number of variables that we 
 % solve for are equal to the sum of the tensor dimensions.
figFlg = true;                                                             % display summary figure flag
tensorSize = length(eigValues);                                            % tensor size; ie the number of different dimensions of the tensor
dim = nan(tensorSize,1);                                                   % tensor dimensions
tensorIxs = 1:tensorSize;                                                  
threshold = -10;                                                           % if an eigenvalue is below this threshold it is considered 0. 
for x = tensorIxs
   dim(x) = length(eigValues{x});
end
preScale = (sum(eigValues{1})./mean(dim));                                 % prescale the eigenvalues for numerical stability
logEigValues = cell(tensorSize,1);                                         % the log of the eigenvalues
optDim = nan(tensorSize,1);                                                % true number of variables that we solve for, which is equal to the sum of the ranks of the marginal covariances                                                
for x = tensorIxs
    logEigValues{x} = log(eigValues{x}./preScale);
    logEigValues{x} = logEigValues{x}(logEigValues{x}>threshold); % eigenvalues should be order apriori
    optDim(x) = length(logEigValues{x});
end

%% instead of solving for the largrangians directly we optimize latent variables that is equal to the log of the lagrangians
%% initialization of the optimization variables                                         
logL0 = cell(tensorSize,1);                                                
for x = tensorIxs
    nxSet = tensorIxs(tensorIxs~=x);
    logL0{x} = log(sum(optDim(nxSet)))-logEigValues{x};
end    

%% this is the optimization step
maxiter = 10000;                                                           % maximum allowed iterations
[logL, logObjperIter] = minimize(vertcat(logL0{:}) ,...
    'logObjectiveMaxEntropyTensor' , maxiter, logEigValues);               % this function performs all the optimzation heavy lifting
L = exp(logL);                                                             % calculates the optimal Largrangians from the latent by taking the exponential
Lagrangians = cell(tensorSize,1);                                          % save the lagrangians to the output 
for x = tensorIxs
    Lagrangians{x} = [L(sum(optDim(1:x-1))+(1:optDim(x)));...
        Inf*ones(dim(x)-optDim(x),1)]./preScale;                           % add the lagrangians known solution (Infinity) of the zero marginal covariance eigenvalues (if low rank)
end
%% caculates the eigenvalues of maximum entropy covariance matrix from the lagrangians
D = 1./diagKronSum(Lagrangians);                                           
%% save and display summary
objCost = objectiveMaxEntropyTensor(vertcat(Lagrangians{:}), eigValues);
logObjCost = logObjectiveMaxEntropyTensor(logL, logEigValues);
Summary.Lagrangians = Lagrangians;
Summary.logObjperIter = logObjperIter;
Summary.objCost = objCost;
fprintf('\n log objective: \n')
fprintf(' - gradient inconsistency with numerical gradient = %.5f \n',...
    checkgrad('logObjectiveMaxEntropyTensor', randn(sum(optDim),1), 1e-5, logEigValues))
fprintf(' - function value at convergence = %.5f \n', logObjCost)


fprintf('\n original objective: \n')
fprintf(' - gradient inconsistency with numerical gradient = %.5f \n',...
    checkgrad('objectiveMaxEntropyTensor', abs(randn(sum(dim),1)), 1e-5, eigValues));
fprintf(' - function value at convergence = %.5f \n', objCost)

if figFlg
figure;
plot(logObjperIter, 'r.-');
xlabel('Iteration #');
ylabel('objective function value');
ylim([0 1])
set(gca, 'FontSize', 16)
end
end
