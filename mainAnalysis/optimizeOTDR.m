function [Q, D, f, flag] = optimizeOTDR(R, codedParams, trCountmtx, spcConstraint, orthFlg)
%{
Optimizes OTDR.

Args:
  R: P x 1 cell array specifying firing rates tensors each of 
    size (T x N x C) (P: coded parameters, T:times, N: neurons, C: conditions).
  codedParams: P x 1 cell array specifying encoded parameter for each condition.
  trCountmtx: C x P matrix where each row is a unique condition and each
       column is a task variable (i.e. regression predictor) for which an
       sRA will be computed.
  spcConstraint: matrix (N x V) defining V orthogonal, N-dimensional
       vectors, the space spanned by which will constrain the sRAs. Leave
       empty to apply no constraint. 
  orthFlg: P x 1 cell array specifying the coded parameters for which the
       corresponding static regression axes (sRA) should be orthogonalized. 

Returns:
  Q: matrix of size N x P representing sRAs encoding P task parameters 
    (normalized to 1 norm). 
  D: scale of sRA.
  f: loss function.
  flag: Boolean whether algorithm converged or not.
%}
if ~exist('orthFlg','var')
    orthFlg = true;
end

% set tolerance for convergence
tol = 10^-5;

% initialize
flag = 0; % flag for success of optimization

N = size(R{1},2);
C = size(codedParams{1},1);
if isempty(spcConstraint)
    projSpc = eye(N);
else
    spcConstraint = orth(spcConstraint); % make sure contraint space is orthonormal;
    projSpc = spcConstraint*spcConstraint'; % projection matrix on the space
end

if isempty(trCountmtx)
    trCountmtx = ones(N, C);
end
numSignals = [];
% orthSignals = [];
for i = 1:length(R)
    numSignals(i) = size(codedParams{i},2);
    % DK -- this appears to be an effort to linearize orthFlg, which is a
    % cell array of potentially different sized vectors. Since vectors
    % of different sizes cannot be combined into a matrix, we perform this
    % linearization below.
%     orthSignals(:, i) = orthFlg{i};
end
Q = projSpc(:,1:sum(numSignals));
orthSignals = horzcat(orthFlg{:}); % DK -- see above
% orthSignals = orthSignals(:);
%% random initialization within data space

if sum(orthSignals)>1
   maxIter = 30000;
%    D = ones(sum(numSignals),1);
   Q0 = zeros(N, length(numSignals));
   stp = 0.1;
   for k = 1:maxIter
       for i = 1:length(numSignals)
           if i > 1
              mskQ = sum(numSignals(1:i-1))+1:sum(numSignals(1:i));
           else
               mskQ = 1:numSignals(i);
           end
           [f_i(i), gradBetas] = TDRobjectiveFn(permute(R{i}, [3, 2, 1]), [codedParams{i} ones(C,1)], permute([Q(:,mskQ), Q0(:,i)], [3 1 2]), trCountmtx);
           gradBetas = squeeze(gradBetas);
           Q(:,mskQ) = Q(:,mskQ) - stp * gradBetas(:,1:numSignals(i));
           Q0(:,i) = Q0(:,i) - stp * gradBetas(:,numSignals(i)+1);
           Q = projSpc*Q; % Q projected into constrained space
           
           Q = orthFn(Q, orthSignals);
           
%            % apply orthogonization constraint
% %            if orthFlg
%            orthFn(Q, orthFlg)
%                [u,~,v] = svd(Q, 0);
%                D = diag((Q')*(u*v')); % note that we used to use unconstrained Q here,
%                % but now are constraining Q by the space
%                % defined in projSpc. This should not change
%                % things because u*v' applies the space
%                % constraint to Q regardless of whether Q
%                % was already space-contrained.
%                Q = (u*v')*diag(D);
% %            end
       end
       % DK -- I'm guessing this is supposed to be a reassignment of "f",
       % not "f_t"
       f(k) = mean(f_i);
%         f(k) = mean(f);

       if k>200 && abs(mean(f(k-200:k-100)) - mean(f(k-99:k))) < tol
           break
       end
   end
   if k >= maxIter
       warning('Maximum number of iterations (%d) reached. Optimization may not have converged.',maxIter);
       flag = abs(mean(f(k-200:k-100)) - mean(f(k-99:k)));
   end
[D, Q] = normVects(Q);
  % convert sRAs (Q) to unit vectors 
%   if orthFlg
%       Q = Q * diag(1./D);
%   else
%       [D, Q] = normVects(Q);
%   end
%%
  %{ 
manifold = stiefelfactory(N, numSignals, 1);
problem.M = manifold;

problem.cost = @(Q) get_output( 'oTDRobjFn' , 1 , Q, R, codedParams, trCountmtx); 
problem.egrad = @(Q) get_output( 'oTDRobjFn' , 2 ,Q, R, codedParams, trCountmtx);
% options.verbosity = 0;
options.tolgradnorm = 1e-5;
options.maxiter = maxIter; 
warning('off', 'manopt:getHessian:approx');

% figure;
% checkgradient(problem);

[Q, fQ , info, options] = trustregions(problem , Q , options );
[f, gradQ, Rscale1, Cscale1, D] = oTDRobjFn(Q, R, codedParams, trCountmtx);
f_t = [info.cost];
figure;
hold on
plot(f_t)
xlabel('iter')
ylabel('cost fn.')

Q = Q*diag(sign(D));
D = abs(D);
  %}


else
    for i = 1:length(numSignals)
        if i > 1
           mskQ = sum(numSignals(1:i-1))+1:sum(numSignals(1:i));
        else
           mskQ = 1:numSignals(i);
        end
        [dRAs, ~, Summary] = runTDR(permute(R{i}, [3, 2, 1]), N, [codedParams{i} ones(C,1)], trCountmtx, false);
        Q(:,mskQ) = squeeze(dRAs(:, :, 1:numSignals(i)));
        f_i(i) = Summary.f_t;
    end
    f = mean(f_i);
    [D, Q] = normVects(Q);
end

end

function Q = orthFn(Q, orthSignals)
    orthSignals = orthSignals>0;
    Qo = Q(:, orthSignals);
    [u,~,v] = svd(Qo, 0);
    D = diag((Qo')*(u*v')); 
    Qo = (u*v')*diag(D);
    Q(:, orthSignals) = Qo;
end
