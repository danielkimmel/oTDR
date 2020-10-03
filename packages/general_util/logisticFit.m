function [b, ci, bp, bcov, err, X, pSuc, coiSet, errAll, figH, gof, posSat, poiSet] = logisticFit(nSuc,nTrial,X,varargin)
% [b ci bp bcov err X pSuc coiSet errAll figH] = logisticFit(nPos,nTrial,X,varargin)
% 
% Performs maximum likelihood estimation (MLE) of the parameters B that
% maximize the following log likelihood:
% L = P * log(nSuc) + (1-P) * log(nTrial-nSuc)
% summed across all conditions that gave rise to nSuc and nTrial.
% P is the probability of success on a single Bernoilli trial for a given
% set of predictor conditions.  P is given by:
% P = 1 / (1 + exp(-B*X))
% where B are the set of parameters we're attempting to maximize, and X is
% a set of predictor values across all conditions, and therefore B*X is the
% linear combination of predictor values weighted by B.
%
% ACCEPTS:
% nSuc = M x 1 vector of the number of successful trials in a given
%        condition from 1 to M.
% nTrial = M x 1 vector of the number of trials in a given condition from 1
%          to M. 
% X = M x p matrix of predictor values, with each condition given rowwise
%     from 1 to M, and each predictor given column wise from 1 to p.  The
%     first column must be all ones to serve as a constant.  In addition to
%     the main predictors, additional columns can be added as interaction
%     terms between predictor values by taking the product of two or more
%     predictors.
% [Optional name,value pairs -- see below].
%
% RETURNS:
% B = 1 x p vector of parameters estimates. Note that for multiple
%       predictors of interest (see posPOI), the parameter associated with
%       each predictor are returned in the column order the predictors were
%       pass in X.
% CI = 2 x p matrix of confidence intervals.
% BP = 1 x p vector of p values for parameter estimates.
% BCOV = p x p matrix of covariance values between the parameters.
% ERR = negative log likelihood function evaluated at the MLE parameters.
% X = predictor matrix (see above) after transformation, if any, as used
%     during fitting 
% PSUC = c x M matrix of the proportion of successes as a function
%        of c values of the covariate given by posCOI (see below) at M
%        unique combinations of the predictors of interest (POI) as
%        specified in POISET.
% COISET = c x 1 vector of unique values of the covariate given in posCOI
%          (see below). Useful because COISET contains the transformed
%          values used in the logistic fit.
% ERRALL = vector of negative log likelihoods at each iteration of
%          minimization routine. 
% FIGH = handle to figure if plotted.
% GOF = structure with goodness of fit info, including
%       .ssr = sum of squared residuals
%       .chi2 = sum of normalized squared residuals, which is chi2
%       distributed
%       .df = degrees of freedom of chi2 distribution, given by nCond -
%       nPredictor (less the constant) - 1
%       .pModel = p-value of model.
% POSSAT = position in B vector that contains saturation parameter. Empty
%       when form does not contain a saturation term.
% POISET = M x N matrix of M sets of unique predictor values for the N
%       predictors specified in POSPOI. These M sets are used to define
%       the M columns in PSUC.
%
% OPTIONAL NAME,VALUE PAIRS (DEFAULT):
% alpha = level needed to achieve signifiance of a parameter estimate
%         differing from zero.  Also sets width of confidence intervals as
%         CI = alpha * 100%.  (0.95)
% bitEstStart = logical whether to estimate the starting conditions of the
%               minimization routine empirically based on the data (=true)
%               or based on the null hypotheses (=false). Currently only
%               applies to the saturation parameter and its offset.
%               Defaults to false.
% bitPlot = logical whether to plot logistic fit.  Requires posCOI.
%           (false).
% form = string specifying form of the logistic link function. ('basic').
%        'basic': P = 1 / (1 + exp(-B*X))
%        'sat': P = B(end) / (1 + exp(-B(1:end-1)*X)) -- allows function to
%               saturate at level less than 1.  Note that the linear
%               weights are represented in B(1:end-1).  The saturation
%               level is set by the last element in B.  B(end) is bounded
%               on [0,1]
%        'sat2': P = (B(end-1)+B(end)*X(:,posPOI)) / (1 + exp(-B(1:end-2)*X)) 
%                allows function to saturate at level less than 1 and for
%                this level to vary depending on the value of a BINARY
%                predictor given in the column of X specified by posPOI
%                (see below).  The saturation level is set by the
%                penultimate element in B, or B(end-1). The potential
%                additive effect of predictor X(:,posPOI) is given by
%                B(end). B(end-1) is bounded on [0,1], B(end) is bounded on
%                [-1 1], and B(end-1) + B(end) is bounded on [0,1].  Note
%                that the linear weights are represented in B(1:end-2).
%        'sat2cov': P = (B(end-1)+B(end)*X(:,posPOI)+B(end-2)*X(:,posCov4Sat)) / (1 + exp(-B(1:end-3)*X))
%                Same as 'sat2' except that Sat parameter B(end-2) can be
%                modified by an additional covariate (see posCov4Sat). The
%                influence of the covariate is additive, as is the
%                influence of the predictor (see posPOI). The influence of
%                the covariate is bounded by [-1 1], with the sum of
%                B(end-2) + B(end-1) + B(end) bounded on [-1 1].
% bitSymmetricLapse = logical on whether to treat lapse rate as symmetric.
%          Specifically, the logistic is assumed to saturate at some level
%          delta (see FORM), where delta < 1, the difference (1 - delta) is
%          considered the lapse rate, L. When bitSymmetricLapse == False
%          (default), these lapses are assumed to only apply at the top
%          end, i.e., p = logistic() spans [0 delta]. When
%          bitSymmetricLapse == True, the lapses are assumed to occur at
%          the bottom and top ends with equal frequency, i.e., p =
%          logistic() spans [L/2 delta+L/2]. Note that in either case, the
%          returned parameter delta refers to the range of the logistic,
%          which is independent of the minimum value. bitSymmetricLapse
%          applies only when form is not 'basic'.
% nIter = number of iterations with random start points that should be
%         attempted before returning the parameter set with the lowest
%         negative log likelihood.  When nIter = 1, the start points are
%         the null hypothesis (not random). (1). 
% posCOI = scalar specifying the column position in X of the predictor that
%          serves as a covariate of nSuc/nTrial. This predictor is used to
%          more accurately estimate the saturation parameter when present,
%          and is used when bitPlot==true as the x-axis. Defaults to [].
% posPOI = vector specifying the column position(s) in X of the predictor(s) of
%          interest.  This predictor must be a binary variable.  It is
%          required with certain forms of the logistic (e.g., 'sat2' -- see
%          above) that estimate an additive effect of this predictor.  It
%          may also be exploited when plotting to emphasize a comparison
%          between the presence or absence of that predictor. Defaults to
%          empty [].
% posCov4Sat = scalar specifying the column position in X of the covariate
%          that will influence the height of the curve saturation as a
%          function of the covariate (e.g., trial number).  Applies only in
%          saturating forms of the logistic function for which special
%          allowance has been made, i.e., 'sat2cov'.  Defaults to empty =
%          [];
% posCov4Gamma = scalar specifying the column position in X of the
%          covariate that you wish to transform by raising it to the gamma
%          -- a new parameter.  When this field points to a column, the
%          gamma parameter is automatically included and is the
%          last parameter in the output vectors B, etc.  Also, expected to
%          be the last parameter in any optional input vectors. When the
%          field is empty, the gamma parameter is not included (i.e., set
%          to 1). 
% fixValue = M x 1 vector specifying values of parameters to fix and not
%          allowed to vary in the MLE optimization. This is achieved by
%          including the fixed value as the upper and lower bounds. The
%          vector must be as long as the number of parameters, or that
%          returned as B. The total number of parameters is not obvious
%          from the input params, and must be computed based on the form of
%          the logistic or by running the function once without fixed
%          params and obtaining the size of B. The fixed values should be
%          in the position of the parameter in B. Use NaN if you do not
%          wish to fix a given parameter value. Leave empty [] to not fix
%          any parameter value.
% randSD = standard deviation of normal distribution used to generate
%          random start points for minimization (see nIter).
% transform = string specifying what transformation, if any, should be
%             applied to the predictors in X. Note that all transformations
%             spare the first predictor, which is for computing the
%             constant. ('scale') 
%             'none' or [] = no transfomration applied.
%             'scale' = predictors are scaled to live on range [0,1].
%                       Specifically, X_scale = (X - min(X)) ./ range(X),
%                       columnwise.
%             'z' = predictors are z-transformed: 
%                   X_z = (X - mean(X)) ./ sd(X), columnwise
%             'center' = predictors are centered about zero: 
%                        X_center = X - mean(X), columnwise.
% bitUnwindTransform = logical on whether to reverse the effect of
%           transformation of X on b, ci, 
%           and bcov (note that I'm only really sure that the reversal works
%           correctly on b, not necessarily the other vars). This is useful if we
%           want the absolute value of the parameters to map onto the oringinal
%           predictors X. However, if we want to compare the relative strength of the
%           parameters independent of the values in X, then we do not want to unwind
%           the transformation. Default = false.
%
% algorithm = string specifying method of optimization by FMINCON.
% maxFunEvals = scalar of number of function evaluations allowed by
%               minimization routine. (2000)
% maxIter = scalar of number of iterations allowed by minimization routine.
%           (2000)
% tolFun = tolerance parameter of minimization routine (1e-4).
% tolX = tolerance parameter of minimization routine (1e-4).
%
% regularization = string on type of regularization to perform:
%           'none' (Default), 'lasso' (L1 regularization), 'ridge' (L2
%           regularization -- not currently supported) 
% lambda = scalar value of coefficent on regularization term. Default = 0.
% bIn = Array of parameter values to use in place of finding fit parameter
%           values This is of use when we wish to compute the error of a
%           given set of input params. Default = [].
%
% Daniel Kimmel, January 21, 2009.
% Update June 5, 2011 -- Includes multiple predictors of interest. 
% Updated 2015 Dec 23 to fix parameter values -- D Kimmel
% Updated 2016 Jul 21 to include regularization
% Updated 2016 Oct 17 to accommodate unique combinations of predictors of
%       interest and to output POISET.
% Update April 9, 2020 -- includes bitSymmetricLapse for modeling symmetric
%       lapse rate
%% default parameter values

alpha = 0.95;
bitEstStart = false;
bitPlot = false;
form = 'sat2';
bitSymmetricLapse = false;
nIter = 1;
posCOI = [];
posPOI = [];
poiSet = [];
posCov4Sat = [];
posCov4Gamma = [];
fixValue = [];
randSD = 0.1;
transform = 'scale';
% logical on whether to reverse the effect of transformation of X on b, ci,
% and bcov (note that I'm only really sure that the reversal works
% correctly on b, not necessarily the other vars). This is useful if we
% want the absolute value of the parameters to map onto the oringinal
% predictors X. However, if we want to compare the relative strength of the
% parameters independent of the values in X, then we do not want to unwind
% the transformation.
bitUnwindTransform = false; 

algorithm = 'interior-point';
maxFunEvals = 2000;
maxIter = 2000;
tolFun = 1e-4;
tolX = 1e-4;

% %%% REGULARIZATION PARAMS
regularization = 'none'; % string on type of regularization to perform:
    % 'none', 'lasso' (L1 regularization), 'ridge' (L2 regularization)
lambda = 0; % scalar value of coefficent on regularization term

bIn = []; % Array of parameter values to use in place of finding fit parameter values
          % This is of use when we wish to compute the error of a given set
          % of input params.

%% collect optional input values
warnopts(assignopts(who, varargin));

% instantiate:
figH = [];
coiSet = [];

% MODIFY INPUTS:
% We sort the posPOI index so that the order of the saturation offset
% parameter estimates remains constant, and the order of conditions in pSuc
% remains constant ,regardless of the order of the posPOI index supplied.
posPOI = sort(posPOI);

% do the same for posCOI, though this shouldn't really matter
% since we only expect a single entry in posCOI.
posCOI = sort(posCOI);

%% instantiate

% these vars may not be set and will cause an error if not defined.
pSuc = [];

%% checks

% certain forms require a predictor of interest be specified:
if strcmpi(form,'sat2')
    if isempty(posPOI) 
        error('Must specify a predictor of interest (posPOI) when using the %s form of the logistic',form);
    end
end

% checks on posPOI
if ~isempty(posPOI)
    % posPOI must point to a predictor with binary data
    for i = 1:length(posPOI)
        if ~isequal([0 1]',unique(X(:,posPOI))) && ~isequal([-1 1]',unique(X(:,posPOI)));
            error('posPOI must point to a BINARY predictor with both 1s and 0s')
        end
    end
    
    % posPOI cannot point to the constant column. This check should be
    % superceded by the previous check for BINARY predictor.
    if posPOI == 1
        error('posPOI must point to a predictor other than the constant term')
    end
    
    % WE HAVE ELIMINATED THIS RESTRICTION. I'M NOT SURE WHY WE HAD IT. MAY
    % BE IMPORTANT FOR MICROSTIM, OR FOR SAT2 FORMS.
    % predictors of interest are mutually exclusive
%     if any(sum(X(:,posPOI),2) > 1)
%         error('Predictors of interest must be mutually exclusive')
%     end
end

% predictor matrix must have constant term:
if ~all(X(:,1) == 1)
    error('First column of predictor matrix X must be all 1s to serve as a constant term')
end

% nSuc and nTrial should be m x 1 vectors:
if all(size(nSuc) > 1)
    error('nSuc must be an M x 1 vector');
elseif size(nSuc,1) < size(nSuc,2)
    % transpose from row to column vector:
    nSuc = nSuc';
end
if all(size(nTrial) > 1)
    error('nTrial must be an M x 1 vector');
elseif size(nTrial,1) < size(nTrial,2)
    % transpose from row to column vector:
    nTrial = nTrial';
end

% length of nSuc and nTrial should match X
if size(nSuc,1) ~= size(X,1)
    error('nSuc and X must have same number of rows');
end
if size(nTrial,1) ~= size(X,1)
    error('nTrial and X must have same number of rows');
end

% if using regularization, must provide lambda
if ~strcmp('none',regularization) && isempty(lambda)
    error('LAMBDA input cannot be empty when using regularization')
end
    

%% clean up inputs

% remove conditions with no trials 
foo = nTrial == 0;
X(foo,:) = [];
nSuc(foo) = [];
nTrial(foo) = [];
clear foo

%% transform predictors if requested:

if ~isempty(transform) && ~strcmp(transform,'none')
    % CONSIDER: remove special parameter of interest. Since it's binary, it shouldn't
    % be scaled.
    
    % make temporary matrix
    XTemp = X(:,2:end);
    
    switch transform
        case 'scale'
%             XTemp = scaleParam(XTemp); % use this function if you want a
                                         % range other than [0 1]
            % if range is 0, then dividing by range results in NaN
            XSigma = range(XTemp);
            XSigma(XSigma==0) = 1;
            XMu = min(XTemp);
            XTemp = (XTemp - repmat(XMu,size(XTemp,1),1)) ./ repmat(XSigma,size(XTemp,1),1);
        case 'center'
            XMu = mean(XTemp);
            XSigma = 1;
            XTemp = XTemp - repmat(mean(XTemp),size(XTemp,1),1);
        case 'z'
            [XTemp,XMu,XSigma] = zscore(XTemp);
        otherwise
            error('Transformation %s not recognized.',transform);
    end
    
    X(:,2:end) = XTemp;
    clear XTemp
else
    XSigma = 1;
    XMu = zeros(1,size(X,2)-1);
end

%% Set search contraints

% example constraints for a sat2 form with size(X,2) == 3 (3 predictor columns):
% 0 <= b4 <=1 (saturation) 
% -1 <= b5 <= 1 (stim influence on saturation)
% 0 <= b4 + b5 <= 1

% Only need contraints on linear combinations of parameters for sat2 or
% greater form: 
if strncmpi(form,'sat2',4)
    % build weight matrix, with two rows for every predictor of interest.
    bWeight = zeros(2*length(posPOI),size(X,2)+1+length(posPOI));
    
    % weights on parameters: no influence of standard beta parameters
    % (leave weights as zero), but put equal positive influence of
    % b(end-1) as negative influence of b(end).
    % first set weight for saturation parameter
    bWeight(:,size(X,2)+1) = repmat([1; -1],length(posPOI),1);
    % next, iteratively set weight for each predictor of interest
    for i = 1:length(posPOI)
        bWeight(1+(i-1)*2:2+(i-1)*2, size(X,2)+1+i) = [1; -1];
    end
    
    % bound on linear combination of coefficients, so that bWeight * b <=
    % bBound  
    bWeightBound = repmat([1; 0],length(posPOI),1);
    
    % if form also includes a covariate that can influence total saturation
    % point, include it in the linear combination as a negative influence:
    if strcmpi(form,'sat2cov')
        bWeight(:,end+1) = repmat([1; -1],length(posPOI),1);
    end
else
    bWeight = [];
    bWeightBound = [];
end

% upper and lower bound on stardard beta parameters:
bUBound = Inf(1,size(X,2));
bLBound = -bUBound;

% for special forms, we have additional bounds on saturation parameters
if strncmpi(form,'sat',3)
    bUBound = [bUBound 1];
    bLBound = [bLBound 0];
end
% ... and on additive parameters
if strncmpi(form,'sat2',4)
    bUBound = [bUBound ones(1,length(posPOI))];
    bLBound = [bLBound -ones(1,length(posPOI))];
end
% ... and on additional additive parameters
if strcmpi(form,'sat2cov')
    bUBound = [bUBound 1];
    bLBound = [bLBound -1];
end

% add weights and bounds for gamma param
if ~isempty(posCov4Gamma)
    % only if being used already
    if ~isempty(bWeight)
%         bWeight(:,end+1) = [0;0];
        % NEW 14 May 2013: it appears this was written assuming there would
        % be only 1  posPOI and thus 2 rows in bWeight.  But with multiple
        % posPOIs (as in Kirby with multiple stim parameters), we need
        % a whole column of bWeight.
        bWeight(:,end+1) = zeros(size(bWeight,1),1);
    end
    bUBound(end+1) = Inf;
    bLBound(end+1) = -Inf;
end

%% compute prob of success
% only can compute if posCOI is given
if ~isempty(posCOI)
    
    % get covariate of interest (x-axis).
    [coiSet,~,coiInd] = unique(X(:,sort(posCOI)));
    
    % prepare "Single" pSuc in case of common saturation point (form ==
    % 'sat')
%     pSucSingle = NaN(length(coiSet),1);
    pSucSingle = NaN(length(coiSet),1);
    % compute pSucSingle for each value in coiSet.
    for i = 1:length(coiSet)
        bitFoo = X(:,posCOI)==coiSet(i);
        pSucSingle(i,1) = sum(nSuc(bitFoo)) ./ sum(nTrial(bitFoo));
    end
    clear bitFoo
    
    % compute pSuc seperately depending on the binary predictor of interest
    % given in posPOI. 
    if ~isempty(posPOI)
        
        % determine set of combinations of the predictor(s) of interest. Do
        % not use the UNIQUE() function, because we want the order of each
        % predictor to be the same as the column order passed in X.
        poiSet = zeros(length(posPOI)+1,length(posPOI));
        for i = 1:length(posPOI)
            poiSet(i+1,i) = 1;
        end
        % I actually think we can use unique. It won't change the column
        % order in poiSet. Of course, the columns in pSuc will need to
        % reference the rows in poiSet to understand each condition.
        [poiSet,~,poiInd] = unique(X(:,posPOI),'rows');
 
        pSuc = NaN(length(coiSet),size(poiSet,1));

        % compute pSuc for each value in coiSet.
        for i = 1:size(coiSet,1)
            % compute pSuc for each combination of values across all
            % predictors of interest
            for j = 1:size(poiSet,1)
                % find conditions that match particular value of covariate
                % and predictor(s) of interest.
%                 bitFoo = coiInd==i & poiInd==j;
                bitFoo = coiInd==i & ismember(X(:,posPOI),poiSet(j,:),'rows');
%                 bitFoo = ismember(X(:,posCOI),coiSet(i,:),'rows') & ismember(X(:,sort(posPOI)),poiSet(j,:),'rows');
                pSuc(i,j) = sum(nSuc(bitFoo)) ./ sum(nTrial(bitFoo));
            end
        end        
    else
        pSuc = pSucSingle;
    end
    clear bitFoo;
end

%% log likelihood maximization

switch regularization
    case 'none'
        % define negative Log likelihood (summed across all conditions)
        nll = @(b) -(nansum(nSuc .* log(logistic(b,X,form,'posPOI',posPOI,'posCov4Sat',posCov4Sat,'posCov4Gamma',posCov4Gamma,'bitSymmetricLapse',bitSymmetricLapse)) ...
            + (nTrial-nSuc) .* log(1 -  logistic(b,X,form,'posPOI',posPOI,'posCov4Sat',posCov4Sat,'posCov4Gamma',posCov4Gamma,'bitSymmetricLapse',bitSymmetricLapse))));
    case 'lasso'
        % define negative Log likelihood (summed across all conditions)
        nll = @(b) -(nansum(nSuc .* log(logistic(b,X,form,'posPOI',posPOI,'posCov4Sat',posCov4Sat,'posCov4Gamma',posCov4Gamma,'bitSymmetricLapse',bitSymmetricLapse)) ...
            + (nTrial-nSuc) .* log(1 -  logistic(b,X,form,'posPOI',posPOI,'posCov4Sat',posCov4Sat,'posCov4Gamma',posCov4Gamma,'bitSymmetricLapse',bitSymmetricLapse)))) ...
            + lambda*norm(b,1);
    case 'ridge'
        error('Ridge regularization not yet supported')
    otherwise
        error('Regularization type %s not recognized.',regularization)
end

% starting conditions
% null hypothesis is that standard beta parameters are 0:
startCond = zeros(1,size(X,2));

if strncmpi(form,'sat',3)
    % saturation parameter is 1:
    startCond(1,end+1) = 1;
    
    % get better estimate if posCOI given and requested -- saturate at
    % maximum nSuc/nTrial, not 1:
    if bitEstStart && exist('pSuc','var')
        % if a single saturation point (form == 'sat'), use average across
        % saturation points:
        if strcmpi(form,'sat')
            startCond(1,end) = max(pSucSingle);
        else
            % isolate case when POI == 0, which we assume is the first
            % column in pSuc (safe assumption given how pSuc is constructed
            % -- see above)
            startCond(1,end) = max(pSuc(:,1));
        end
    elseif bitEstStart && ~exist('pSuc','var')
        error('Position of Covariate not given and therefore better starting parameter estimates could not be computed as requested.');
    end
end
if strncmpi(form,'sat2',4)
    
    % get better estimate if posCOI given (i.e., pSuc computed) and
    % requested -- modify saturation point by difference between max when
    % posCOI == 0 and posCOI == 1
    if bitEstStart && exist('pSuc','var')
        % isolate cases when POI == 1 for each predictor of interest. Since
        % pSuc is in order of the predictors as passed in X, we can go as
        % in order of pSuc:
        for i = 1:length(posPOI)
            startCond(1,end+1) = max(pSuc(:,1+i)) - max(pSuc(:,1));
        end
    else 
        % additive effect on saturation point is zero
        startCond(1,end+1:end+length(posPOI)) = 0;
    end
end
if strcmpi(form,'sat2cov')
    % covariate-dependent effect on saturation point is zero
    startCond(1,end+1) = 0;
    
%     % get better estimate if posCOI given (i.e., pSuc computed) and
%     % the estimate is requested -- modify saturation point by difference in
%     % pSuc between lowest and highest value of covariate.
%     if bitEstStart && exist('pSuc','var')
%         % isolate case when POI == 1
%         [~,foo] = max(X(:,posCov4Sat));
%         [~,goo] = min(X(:,posCov4Sat));
%         startCond(1,end) = mean(pSuc(foo,:),2) - mean(pSuc(goo,:),2);
%         clear foo goo
%     end
end

% add starting point for gamma param
if ~isempty(posCov4Gamma)
    startCond(end+1) = 1;
end

% predimension vars
b = NaN(nIter,size(startCond,2));
errAll = NaN(nIter,1);
exitflag = NaN(nIter,1);
% output = struct([]);
% lambda = struct([]);
grad = NaN(size(startCond,2),nIter);
hessian = cell(nIter,1);


% Do the following outside the parfor loop:
% fix specific parameter values
if ~isempty(fixValue)
    if length(fixValue) ~= size(b,2)
        error('To fix parameter values, the vector fixValue must be same length as B, the number of parameters')
    end
    foo = ~isnan(fixValue);
    bLBound(foo) = fixValue(foo);
    bUBound(foo) = fixValue(foo);
    if ~isempty(bWeight) || ~isempty(bWeightBound)
        error('When fixing parameter values, support for fixaing bWeight and bWeightBound has not yet been considered or implemented.')
    end
else
    foo = [];
end

%%% DECIDE WHETHER TO FIND PARAM VALUES, OR COMPUTE ERROR ON INPUT VALUES
if ~isempty(bIn)
    % compute nll and gof for input parameters. Set other fields to NaN;
    err = nll(bIn);
    [ci bp bcov errAll] = deal([]);
    b = bIn;

else
    % Iterate minimization with random start points
    parfor i = 1:nIter
        % randomize start points only on subsequent interations
        if i > 1
            startCondTemp = startCond + randn(size(startCond)) * randSD; 
        else
            startCondTemp = startCond;
        end

        % fix specific values of startCondTemp
        if ~isempty(fixValue)
            startCondTemp(foo) = fixValue(foo);
        end

        % FMINCON minimization
        [b(i,:), errAll(i), ~, ~, ~, ~, hessian{i}] =  ...
            fmincon(nll,startCondTemp,bWeight,bWeightBound,[],[],...
            bLBound,bUBound,[],...
            optimset('MaxIter',2000,'MaxFunEvals',2000,'Algorithm',algorithm,...
            'Display','off'));
    end
    clear foo

    % use parameters from minimum nll (maximum log likelihood)
    [~, i] = min(errAll);
    b = b(i,:);
    err = errAll(i);
    % exitflag = exitflag(i);
    % output = output(i);
    % lambda = lambda(i);
    % grad = grad(:,i);
    hessian = hessian{i};

    %%% ALT METHOD OF COMPUTING COVARIANCE MATRIX
    % % compute covariance, sd and CI using FMINCON derived parameters, but the
    % % Matlab MLE-based covariance function
    % [bcov] = mlecov(b,nSuc,'nloglf',nll)
    % %     'options',statset('DerivStep',eps^(1/4)));

    % compute covariance, sd, CI using FMINCON hessian:
    bcov = inv(hessian);
    % compute SD:
    sd = sqrt(diag(bcov))';
    % compute CI
    ci = [b - sd * norminv(0.975); ...
        b + sd * norminv(0.975)];

    % compute Probability of getting each param by chance
    % null hypothesis is that standard beta parameters are 0:
    nullH = zeros(1,size(X,2));
    if strncmpi(form,'sat',3)
        % saturation parameter is 1:
        nullH(1,end+1) = 1;
    end
    if strncmpi(form,'sat2',4)
        % additive effect on saturation point is zero
        nullH(1,end+1:end+length(posPOI)) = 0;
    end
    if strcmpi(form,'sat2cov')
        % additive effect of covariate on on saturation point is zero
        nullH(1,end+1) = 0;
    end
    % add null hypothesis for gamma param
    if ~isempty(posCov4Gamma)
        nullH(end+1) = 1;
    end

    % flip null hypothesis and distribution for cases where nullH is to the
    % left of the mean:
    signFlip = sign(nullH -b);
    % probability
    bp = (1-normcdf(nullH .* signFlip,b .* signFlip,sd)) * 2;
    
    %%% Unwind transformation
    if bitUnwindTransform
        % first excluding the constant term
        b(2:end) = bsxfun(@rdivide,b(2:end),XSigma);
        ci(:,2:end) = bsxfun(@rdivide,ci(:,2:end),XSigma);

        % and intercept
        b(1) = b(1) - XMu * b(2:end)';
        ci(:,1) = ci(:,1) - XMu * b(2:end)';
        
        % I have no idea how to transform the covariance matrix:
        bcov = NaN(size(bcov));
        
        % also transform X, skipping constant
        X(:,2:end) = bsxfun(@plus,bsxfun(@times,X(:,2:end),XSigma),XMu);
    end    
end

%% compute sum of squared error
% see "Hallet - Goodness of fit in logistic regression.pdf"
% see http://krex.k-state.edu/dspace/bitstream/handle/2097/530/YingLiu2007.pdf?sequence=1

% derive expected P from fit values
pExp = logistic(b,X,form,'posPOI',posPOI,'posCov4Sat',posCov4Sat,'posCov4Gamma',posCov4Gamma,'bitSymmetricLapse',bitSymmetricLapse);

% To compute the difference between the choice and the probability of
% accept, we we must decide how to code choice. Given that P(accept) is
% limited to the range [0 delta], we encode choice as reject=0,
% accept=delta 

% find position of saturation param
posSat = size(b,2);
if ~isempty(posCov4Gamma)
    posSat = posSat - 1;
end
if strcmpi(form,'sat2cov')
    posSat = posSat - 1;
end
if strncmpi(form,'sat2',4)
    posSat = posSat - 1;
end
if ~strncmpi(form,'sat',3)
    posSat = [];
end

% make sure we didn't screw up in finding saturation param
if ~isempty(posSat) && (b(posSat) > 1 || b(posSat) < 0)
    error('Value for saturation parameter (%g) is unexpected, likely due to mis-positioning parameter with vector of terms',b(posSat));
end

% for each condition, compute difference between nSuc expected and
% observed. Square this difference and sum them across conditions.  
residNorm = NaN(size(pExp));
d = NaN(size(pExp));
nSucExp = nTrial.*pExp;
for c = 1:size(pExp,1)

    
    %%% for pearsons residual
    % residual normalized (standardize) by std deviation of binomial distribution,
    % aka Pearson's Residual.
    residNorm(c,1) = (nSuc(c)-nSucExp(c)) / sqrt(nSucExp(c)*(1-pExp(c)));
    
    %%% for Deviance residual, a la, Hosmer and Lemeshow, 1989
    if nSuc(c) == 0
        d(c) = -sqrt(2*nTrial(c)*abs(log(nTrial(c)*(1-pExp(c)))));
    elseif nSuc(c) == nTrial(c)
        d(c) = sqrt(2*nTrial(c)*abs(log(nSucExp(c))));
    else
        d(c) = sign(nSuc(c)-nSucExp(c))*...
            sqrt(2*(nSuc(c)*log(nSuc(c)/nSucExp(c))+...
            (nTrial(c)-nSuc(c))*log((nTrial(c)-nSuc(c))/(nTrial(c)*(1-pExp(c))))));
    end
end

gof.pExpected = pExp; % model-based probability of success
gof.pCorrByModel = 1 - 1/sum(nTrial) * sum(abs(nSuc-round(nSucExp)));

%%% Pearson's residual
% sum of squared normalized residual = chi sq stat
gof.pearResid.ssr = nansum(residNorm.^2);
if any(isnan(residNorm))
    warning('SOME VALUES IN PEARSON RESIDUALS ARE NAN -- this happens when nTrial=0 for some conditions.')
end
% degrees of freedom = nCond - nPredictors (excluding the intercept) - 1
gof.pearResid.df = size(pExp,1) - (size(X,2) - 1) - 1;
% compute p Value on model
gof.pearResid.pModel = 1-chi2cdf(gof.pearResid.ssr,gof.pearResid.df);

%%% Deviance Residual
% sum of squared deviance is chi-square distributed
gof.deviance.ssd = nansum(d.^2);
gof.deviance.df = size(pExp,1) - (size(X,2) - 1) - 1;
gof.deviance.pModel = 1-chi2cdf(gof.deviance.ssd,gof.deviance.df);

%%% Error

%% plot fits
if bitPlot && ~isempty(posCOI)
    
    % make figure
    figH = figure;
    whitebg(gcf,'k');
    set(gcf,'Color','k','Name',['MLE - form ',form]);
    
    % plot
    [hData,hFit] = CB_psychCurveFitPlotLogistic(coiSet,pSuc,b,X,form,posPOI,posCOI,ci,[],[],bp,err,[],posCov4Sat,[],[],posCov4Gamma);
    
    findfigs
elseif bitPlot && isempty(posCOI)
    warning('Plot requires position of covariate of interest (posCOI). Plot aborted.');
end
