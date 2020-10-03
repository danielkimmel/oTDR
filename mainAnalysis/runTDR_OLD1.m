% This is an older version of runTDR in which lambdas were scaled by the
% predictors separately for each cross-validation, which is incorrect. It
% also has alternative code to use a consistent, non-scaled set of lambdas.
% The code has been consolidated into using a consistent, SCALED set of
% lambdas in the newer runTDR.m

function [dRAs, normdRAs, Summary] = runTDR(dataTensor, numPCs, codedParams, trCountmtx, loocvFlg)
[T, N, C] = size(dataTensor);
%% denoise step
if numPCs < N
    XN = reshape(permute(dataTensor,[1 3 2]), [], N);
    PCs = pca(XN);
    XN = bsxfun(@minus, XN, mean(XN))*(PCs(:, 1:numPCs)*PCs(:, 1:numPCs)');
    dataTensor = permute(reshape(XN', N, T, C), [2 1 3]);
end

%%
K = size(codedParams, 2);
% normalize trial count by total number per neuron, but only if not already
% normalized (since we don't wish to introduce extra machine precision
% error) 
% if abs(sum(trCountmtx,2) - 1) > eps
%     trCountmtx = bsxfun(@rdivide,trCountmtx,sum(trCountmtx,2)); 
% end
sqrtTrCountmtx = sqrt(trCountmtx);
Betas = nan(T, N, K); % beta coefficients
lambda = nan(T, N); %regularization term
err = nan(T, N); % fit error
for t = 1:T
   parfor n = 1:N
       scaledCodedParams = repmat(sqrtTrCountmtx(n, :)', 1, K).*codedParams; % scaled coded parameters by the square root of trial counts
       scaledRbar_tn = sqrtTrCountmtx(n, :)'.*reshape(dataTensor(t, n, :), C, 1); % scaled coded parameters by the square root of trial counts
       if loocvFlg
           [Betas(t, n, :), lambda(t, n), err(t,n)] = crossValRegress(scaledRbar_tn, scaledCodedParams, C);
       else
           Betas(t, n, :) = regress(scaledRbar_tn, scaledCodedParams);
       end
   end
end
[Summary.f_t, ~, Summary.R2_t, Summary.R2_tk] = TDRobjectiveFn(dataTensor, codedParams, Betas, trCountmtx);
Summary.lambda = lambda;
Summary.testError = err;
%% divide beta by magnitude and direction  
[normdRAs, dRAs] = normVects(reshape(permute(Betas, [2, 1, 3]), N, T*K));
dRAs(:, normdRAs<=eps) = 0; %set the direction to zero if the magnitude of the vector is zero
dRAs = permute(reshape(dRAs, N, T, K), [2 1 3]);
normdRAs = reshape(normdRAs, T, K);

%% compute coefficients used for projecting neural data, these equal to transpose(Beta) when Beta is orthogonal. Otherwise, they should be the pseudoinverse.

%%% FOR NOW, we are not confident in the pseudoinverse method, so we will
%%% not apply it here, but may wish to in the future. Note that it only
%%% influences the projection of data onto dRAs. We can can continue to
%%% compare dRAs.
warning('TDR:ProjectOntoNonorthogAxes','Incorrect to project data onto non-orthogonalized, non-transformed dRAs. Still deciding on correct transformation')
Betas_proj = Betas;

% Betas_proj = nan(T, N, K); % Coefficients used for proj
% parfor t = 1:T
%    B = squeeze(Betas(t, :, :));
%    Betas_proj(t, :, :) = ((B'*B)\B')';
% end

%% divide beta for proj by magnitude and direction  
[normdRAs_proj, dRAs_proj] = normVects(reshape(permute(Betas_proj, [2, 1, 3]), N, T*K));
dRAs_proj(:, normdRAs_proj<=eps) = 0; %set the direction to zero if the magnitude of the vector is zero
dRAs_proj = permute(reshape(dRAs_proj, N, T, K), [2 1 3]);
normdRAs_proj = reshape(normdRAs_proj, T, K);
Summary.dRAs_proj = dRAs_proj;
Summary.normdRAs_proj = normdRAs_proj;
end


function [B, lambda, Er] = crossValRegress(R,P, numConds)
% Gamal's lambdas
% lambdas =  [0 2.^[0:15]*0.01 inf];
% DK's lambdas (need to scale given different range
lambdas =  [0 2.^[0:19]*0.01 inf];

T = size(R,1)/numConds;
K = size(P, 2);
R = reshape(R, T, numConds);
P = reshape(P, T, numConds, K);
erTrain = nan(length(lambdas), numConds);
erTest = nan(length(lambdas), numConds);

% DK -- for the set of lambdas when using an "RMS" based error term (see
% below), here was the original set: L = [0 2.^[0:12]*0.01 inf];
% By converting to "norm-2" based error term (see below), the results
% should be same if we let the new set of lambda = L^2 / N, where N is
% number of predictor values used in training set (predictors * (conditions-1).
lambdas_trad = lambdas.^2 / ((numConds-1)*K);
% lambdas = lambdas.^2 / ((numConds-1)*K);

for c = 1:numConds
    msk = false(numConds, 1);
    msk(c) = true;
    RTest = R(:, msk);
    PTest = reshape(P(:, msk, :), T, K);
    
    RTrain =  R(:, ~msk).';
    PTrain = reshape(P(:, ~msk, :), T*(numConds-1), K);
    
    % use gamal's ridge regression
%     [erTrain(:, c), erTest(:, c)] =  ridgeRegress(PTrain,RTrain(:), PTest, RTest, lambdas,lambdas_trad);
    
    % use DK's ridge regression
    [erTrain(:, c), erTest(:, c)] =  ridgeRegress_DK(PTrain,RTrain(:), PTest, RTest, lambdas);
    
%     % compute error for each level of lambda
%     for i = 1:size(betaTemp,2)
%         erTrain(i,c)=norm(PTrain*betaTemp(i)-RTrain(:),'fro')^2;
%         erTest(i,c)=norm(PTest*betaTemp(i)-RTest(:),'fro')^2;
%     end   
    
end
% Here we select the minimum test error as the "error" to report for
% regularization. This emphasizes the generalizability of the betas over
% the goodness-of-fit, which we had previously emphasized by returning the
% error from the final fit that includes all conditions; see
%        [~,~, B] =  ridgeRegress_DK(P, R(:), P, R(:), lambda);
[Er,ix] = min(nanmean(erTest,2));
% [Er,ix] = min(bsxfun(@times, sum(erTest,2), 1./norm(R(:)).^2)); % problem of nan if R =0, but good for plotting

% check that lambda_RMS and lambda_traditional match
if ~isequal(lambdas(ix).^2/((numConds-1)*K),lambdas_trad(ix))
    error('lambda via RMS and lambda_traiditional do not match')
end
lambda = lambdas(ix);
% reconvert lambda_trad, because when using all conditions (not C-1, as for
% cross-validation), the conversion factor changes
lambda_trad = lambda.^2 / (numConds*K); % using traditional lambda instead

R =  R(:, :).';
P = reshape(P(:, :, :), T*(numConds), K);

% use gamal's ridge regression
% [~,~, B] =  ridgeRegress(P, R(:), P, R(:), lambda,lambda_trad);
% DK's regression
[~,~, B] =  ridgeRegress_DK(P, R(:), P, R(:), lambda);

end


function [erTrain,erTest,Xlast]=  ridgeRegress(Atrain,Btrain, Atest, Btest, lambdas,lambdas_trad)
[~,n]=size(Atrain);


erTrain = nan(length(lambdas),1);
erTest = nan(length(lambdas),1);
for i=1:length(lambdas)
    % DK -- here, Gamal is using RMS Beta and squaring lambda.
    lambda = (lambdas(i)*sqrt(mean(Atrain(:).^2)))^2;
    % add a little bit of error to test idea that tiny error in lambda will
    % cause bigger error in X
%     lambda = lambda + rand(1)/1e5;

    % DK -- here, Daniel is implementing the more traditional norm-2 of
    % beta and not squaring lambda. If one uses norm_fro instead of norm-2, 
    % they should be related by:
    % lambda_traditional = lambda_RMS^2 / N, where N is the number of
    % elements in Atrain. 
    % However norm-2 < norm-fro, and technically L2 regularization uses
    % norm-2
%     lambda = lambdas(i)*norm(Atrain(:),2)^2;

    if ~isempty(lambdas_trad)
        lambda_trad = lambdas_trad(i)*norm(Atrain(:),2)^2;

        if abs(lambda - lambda_trad) > 1e-5
            error('Gamal''s Lambda (%0.2f) and DK''s lambda_trad (%0.2f) differ',...
                lambda,lambda_trad);
        end
        % now that we've tested lambda traditional, change lambda to it
        lambda = lambda_trad;
        
    else
        % when lambdas_trad is empty, this means that the incoming lambdas
        % is under the lambdas_trad transformation
        lambda = lambdas(i)*norm(Atrain(:),2)^2;
        
        % confirm that lambda_RMS and lambda_trad are matching 
        
        
        % just confirm
%         if i==1
%             fprintf('lambda is traditional\n')
%         end

    end
    
    if isinf(lambda)
        X = zeros(size(Atrain,2), size(Btrain,2));
    else
    X=(Atrain'*Atrain+diag(lambda*ones(n,1)))\(Atrain'*Btrain);
%     X=(Atrain'*Atrain+diag([lambda*ones(n-1,1); 0]))\(Atrain'*Btrain);
    end
    erTrain(i)=norm(Atrain*X-Btrain,'fro')^2;%./norm(Btrain,'fro')^2; 
    erTest(i)=norm(Atest*X-Btest,'fro')^2;%./norm(Btest,'fro')^2; 
end
Xlast = X;

end

%% DK's own ridge regression based on matlab's ridge but without rescaling

function [erTrain,erTest,Xlast]=  ridgeRegress_DK(Atrain,Btrain, Atest, Btest, lambdas)

[~,n]=size(Atrain);

erTrain = nan(length(lambdas),1);
erTest = nan(length(lambdas),1);

for i=1:length(lambdas)
    if isinf(lambdas(i))
        X = zeros(size(Atrain,2), size(Btrain,2));
    else
        % note that here we may wish to square lambda, since technically
        % the solution takes the form B = (X'*X + L'*L)\(X'*y), 
        % where L = lambda * eye(n), and thus L'*L = lambda^2 * eye(n)
        % See https://en.wikipedia.org/wiki/Tikhonov_regularization
        X=(Atrain'*Atrain+lambdas(i)*eye(n))\(Atrain'*Btrain);
    end
    erTrain(i)=norm(Atrain*X-Btrain,'fro')^2;%./norm(Btrain,'fro')^2; 
    erTest(i)=norm(Atest*X-Btest,'fro')^2;%./norm(Btest,'fro')^2; 
end
Xlast = X;


end % end DK's ridge regression
