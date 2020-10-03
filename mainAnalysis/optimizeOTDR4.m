function [Q] = optimizeOTDR4(dataTensor, codedParams, trCountmtx, spcConstraint)
[T, N, C] = size(dataTensor);
if isempty(spcConstraint)
projSpc = eye(N);
else
    spcConstraint = orth(spcConstraint); % make sure contraint space is orthonormal;
    projSpc = spcConstraint*spcConstraint'; % projection matrix on the space
end

if isempty(trCountmtx)
    trCountmtx = ones(N, C);
end

numSignals = size(codedParams, 2);
Q = orth(randn(N, numSignals));
Beta0 = randn(N, T);
a = ones(T, numSignals);
%% random initialization within data space
maxIter = 2000;
stp = 1;
IT = eye(T);
IK = eye(numSignals);
onesVect = ones(C,1);
As{1} = IT;
As{3} = IK;
for i = 1:maxIter
 i
% %  for t = 1:T
% %     Betas(t, :, :) = permute([Q*diag(a(t,:)) Beta0(:,t)], [3 1 2]);
% %  end
% % [f, gradf, ~,~] = TDRobjectiveFn(dataTensor, [codedParams onesVect], Betas, trCountmtx);
% % Betas = Betas - stp * gradf;
% % Q = squeeze(mean(Betas(:, :, 1:numSignals)));
% % [u,~,v] = svd(projSpc*Q, 0);
% % Q = u*v';
% % Beta0 = squeeze(Betas(:, :, numSignals+1)).';
% % 
% % B = Betas(:, :, 1:numSignals);
% % 
% % for t = 1:T
% %     a(t,:) = diag(squeeze(B(t, :, :))'*Q);
% % end
% % 
% % 
% % 
% % a(t,:) = diag((Q')*(u*v'));
for t = randperm(T)
        Betas = [Q*diag(a(t,:)) Beta0(:,t)];
        [f(i,t), gradf, ~,~] = TDRobjectiveFn(dataTensor(t, :, :), [codedParams onesVect], permute(Betas, [3 1 2]), trCountmtx);
        gradf = squeeze(gradf);
%         Betas(:,j) = Betas(:,j)-stp*gradf(:,j); % update gradient with respect to data R{j}
        Betas(:,1:numSignals) = Betas(:,1:numSignals)-stp*gradf(:,1:numSignals); % update gradient with respect to all data
        Betas(:,numSignals+1) = Betas(:,numSignals+1)-stp*gradf(:,numSignals+1);
        Q = Betas(:, 1:numSignals);
        Beta0(:,t) = Betas(:, numSignals+1);
        [u,~,v] = svd(projSpc*Q, 0);
        a(t,:) = diag((Q')*(u*v'));
        Q = (u*v');
end
%%   
    
end
figure;
plot(mean(f,2))
end











