function [dRAs, normdRAs, Summary] = runTDR(dataTensor, numPCs, codedParams, trCountmtx, loocvFlg)
% TODO: add docstring
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
% normalize trial count by total trial count per condition per neuron, then
% normalize by max normalized trial count across neurons. This ensures that
% all neurons sum to the same number, but that the range of relative trial
% counts is [0, 1] for numerical purposes.
nrm = sum(trCountmtx,2); % save normalizing term
trCountmtx = bsxfun(@rdivide,trCountmtx,nrm);
nrm = nrm .* max(trCountmtx(:)); % update normalizing term.
trCountmtx = trCountmtx ./ max(trCountmtx(:));

% % normalize by max trial count over all neurons
% trCountmtx = trCountmtx ./ max(trCountmtx(:));

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
           [Betas(t, n, :),~,errTemp] = regress(scaledRbar_tn, scaledCodedParams);
           % store error as the squared norm of the residuals (same as used
           % when performing cross validation). But note that this error is
           % final fit error, not "test" error as with cross validation
           err(t,n) = norm(errTemp).^2
       end
   end
end

%%%%
% The following is useful if one wishes to compare error between fits using
% unnormalized vs. normalized trial counts. To do so, the error from
% normalized trial counts should be "unscaled". Do NOT unscale the error if
% comparing between normalized trial counts.
% Unscale error by normalizing factor used to scale trial count
% matrix. Let err and err_bar be the squared norm of residuals
% of the unnormalized and normalized trial-count problems,
% respectively. Let nrm be the trial count scaling factor such
% that
% trCountmtx_normalized = bsxfun(@rdivide,trCountmtx_original,nrm)
% Thus err = err_bar * nrm;
% err = bsxfun(@times,err,nrm');

[Summary.f_t, ~, Summary.R2_t, Summary.R2_tk] = TDRobjectiveFn(dataTensor, codedParams, Betas, trCountmtx);
Summary.lambda = lambda;
if loocvFlg
    % store test error for cross-validated regularization
    Summary.testError = err;
else
    % store final error
    Summary.finalError = err;
end   
% store normalization term
Summary.trialCountNormalization = nrm;
% store codedParams
Summary.codedParams = codedParams;

%% divide beta by magnitude and direction  
[normdRAs, dRAs] = normVects(reshape(permute(Betas, [2, 1, 3]), N, T*K));
dRAs(:, normdRAs<=eps) = 0; %set the direction to zero if the magnitude of the vector is zero
dRAs = permute(reshape(dRAs, N, T, K), [2 1 3]);
normdRAs = reshape(normdRAs, T, K);

end

%% SUBFUNCTION
function [B, lambda, Er] = crossValRegress(R,P, numConds)

T = size(R,1)/numConds;
K = size(P, 2);
R = reshape(R, T, numConds);
P = reshape(P, T, numConds, K);

% Here we define the set of lambda values to sample. 
% lambdas =  [0 2.^[0:19]*0.01 inf]; % (A) set useful when scaling or not scaling lambda and when not squaring lambda below
% lambdas =  [0 2.^[0:13]*0.02 inf]; % (B) close approximation to above (A)
% lambdas = [0 2.^[0:0.5:9.5]*0.1 inf]; % (C) equivalent to sqrt(A). useful if squaring lambda below
% lambdas = [0 2.^[0:1:10]*0.1 inf]; % (D) Comparable to (C), but with fewer steps
lambdas =  [0 2.^[-5:14]*0.01 inf]; % (E) attempting to be simliar to (A) but while scaling trial counts

% take square root of lambda, given that it will be squared below
lambdas = sqrt(lambdas);

% We then scale this set of lambda based on the predictor matrix P, which
% conditions the values of lambda to a range appropriate for the
% predictors. By performing the scaling at this stage, we customize the set
% of lambda for the particular neuron, since different neurons could have
% different predictors (unlikely) and certainly have different scaled
% predictors (as predictors are scaled by trial count). As such, the final
% set of lambdas will different across neurons. Lastly, by performing the
% lambda scaling here, we ensure that the same values of lambda are used
% across cross-validation sets.
% In terms of the particular scaling, we are scaling by the 
% root-mean-squared (i.e., the mean-squared) of all predictors, i.e.,
% across all times, conditions, and predictors. 
lambdas = lambdas*sqrt(mean(P(:).^2));

% instantiate
erTrain = nan(length(lambdas), numConds);
erTest = nan(length(lambdas), numConds);

for c = 1:numConds
    msk = false(numConds, 1);
    msk(c) = true;
    RTest = R(:, msk);
    PTest = reshape(P(:, msk, :), T, K);
    
    RTrain =  R(:, ~msk).';
    PTrain = reshape(P(:, ~msk, :), T*(numConds-1), K);
        
    % use DK's ridge regression
    [erTrain(:, c), erTest(:, c)] =  ridgeRegress(PTrain,RTrain(:), PTest, RTest, lambdas);
        
end

% Here we select the minimum test error as the "error" to report for
% regularization. This emphasizes the generalizability of the betas over
% the goodness-of-fit, which we had previously emphasized by returning the
% error from the final fit that includes all conditions; see
[Er,ix] = min(nanmean(erTest,2));

% select the lambda with the lowest test error
lambda = lambdas(ix);

% reshape response and predictor matrices
R =  R(:, :).';
P = reshape(P(:, :, :), T*(numConds), K);

% repeat regression with selected verion of lambda and return final set of
% beta coefficients.
[~,~, B] =  ridgeRegress(P, R(:), P, R(:), lambda);

end


%% ridge regression based on matlab's ridge 

function [erTrain,erTest,Xlast]=  ridgeRegress(Atrain,Btrain, Atest, Btest, lambdas)

[~,n]=size(Atrain);

erTrain = nan(length(lambdas),1);
erTest = nan(length(lambdas),1);

% loop through values of lambda
for i=1:length(lambdas)
    if isinf(lambdas(i))
        X = zeros(size(Atrain,2), size(Btrain,2));
    else
        % here we square lambda, since the solution takes the form 
        %       B = (X'*X + L'*L)\(X'*y),
        % where L = lambda * eye(n), and thus 
        %       L'*L = lambda^2 * eye(n)
        % See https://en.wikipedia.org/wiki/Tikhonov_regularization
        X=(Atrain'*Atrain+lambdas(i).^2*eye(n))\(Atrain'*Btrain);
    end
    erTrain(i)=norm(Atrain*X-Btrain,'fro')^2;%./norm(Btrain,'fro')^2; 
    erTest(i)=norm(Atest*X-Btest,'fro')^2;%./norm(Btest,'fro')^2; 
end
Xlast = X;


end % end ridge regression
