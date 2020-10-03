function p = logistic(b,X,form,varargin)
% p = logistic(b,X,form,[name,value pairs])
% 
% Logistic function that relates predictors given rowwise in X to a
% probability p according to parameters b.  Used as a link function in
% logistic regression by logisticFit().  
%
% B = 1 x p+n vector of parameters estimates for p beta coefficients and
%     n additional, optional parameters depending on the logistic form that
%     are concatenated to B in the following order: 1 saturation parameter
%     (for 'sat' at minimum), q offset parameters for each predictor of
%     interest (for 'sat2' at minimum), and 1 parameter that scales the
%     saturation point with the covariate of interest (for 'sat2cov' at
%     minimum).
% X = M x p matrix of predictor values, with each condition given rowwise
%     from 1 to M, and each predictor given column wise from 1 to p.  The
%     first column must be all ones to serve as a constant.  In addition to
%     the main predictors, additional columns can be added as interaction
%     terms between predictor values by taking the product of two or more
%     predictors.
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
%       'sat2cov': P = (B(end-1)+B(end)*X(:,posPOI)+B(end-2)*X(:,posCov4Sat)) / (1 + exp(-B(1:end-3)*X))
%                Same as 'sat2' except that Sat parameter B(end-2) can be
%                modified by an additional covariate (see posCov4Sat). The
%                influence of the covariate is additive, as is the
%                influence of the predictor (see posPOI). The influence of
%                the covariate is bounded by [-1 1], with the sum of
%                B(end-2) + B(end-1) + B(end) bounded by [-1 1].
%
% The following are input as optional name, value pairs:
% posPOI = vector specifying the column position(s) in X of the
%          predictor(s) of interest.  This predictor(s) must be a binary
%          variable.  It is required with certain forms of the logistic
%          (e.g., 'sat2' -- see above) that estimate an additive effect of
%          this predictor(s). Defaults to empty [].
% posCov4Sat = scalar specifying the column position in X of the covariate
%          that will influence the height of the curve saturation as a
%          function of the covariate (e.g., trial number).  Applies only in
%          saturating forms of the logistic function for which special
%          allowance has been made, i.e., 'sat2cov' 
% posCov4Gamma = scalar specifying the column position in X of the
%          covariate that you wish to transform by raising it to the gamma
%          -- a new parameter.  When this field points to a column, the
%          gamma parameter is automatically included and expected to be the
%          last parameter in the input vector B.  When the field is empty,
%          the gamma parameter is not included (i.e., set to 1).
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
%
% Daniel Kimmel, January 21, 2010.
% Update March 31, 2010 -- includes sat2cov
% Update June 5, 2011 -- Includes multiple predictors of interest. For
%       clarity, breaks input parameter vector into separate vectors of
%       beta coefficients and other parameters.
% Update April 9, 2020 -- includes bitSymmetricLapse for modeling symmetric
%       lapse rate
%% default values for optional input vars
posPOI = [];
posCov4Sat = [];
posCov4Gamma = [];
bitSymmetricLapse = false;

%% collect optional input values
warnopts(assignopts(who, varargin));

%% modify inputs

% for clarity, break up single vector of parameters into beta coefficients
% and other parameters

% the first p elements of b are the beta coefficients.
nBeta = size(X,2);

% if including a gamma parameter
if ~isempty(posCov4Gamma)
    nParam4Gamma = 1;
else
    nParam4Gamma = 0;
end

% check that length of b is as expected
switch form
    case 'basic'
        if length(b) ~= nBeta + nParam4Gamma
            error('Number of parameters provided in B (%d) is more or fewer than expected (%d) for %s logistic',length(b),nBeta + nParam4Gamma,form);
        end            
    case 'sat'
        if length(b) ~= nBeta + nParam4Gamma + 1
            error('Number of parameters provided in B (%d) is more or fewer than expected (%d) for %s logistic',length(b),nBeta + nParam4Gamma+1,form);
        end
    case 'sat2'
        if length(b) ~= nBeta + nParam4Gamma +1+length(posPOI)
            error('Number of parameters provided in B (%d) is more or fewer than expected (%d) for %s logistic',length(b),nBeta + nParam4Gamma+1+length(posPOI),form);
        end
    case 'sat2cov' 
        if length(b) ~= nBeta + nParam4Gamma +2+length(posPOI)
            error('Number of parameters provided in B (%d) is more or fewer than expected (%d) for %s logistic',length(b),nBeta + nParam4Gamma+2+length(posPOI),form);
        end        
end

% collect the additional parameters, 
% without any checks since these should
% of all been done above.

% the next element after beta coefficients is the saturation parameter
if length(b) >= nBeta+1
    satVal = b(nBeta+1);
end

% the next element(s) is the offset parameter(s) for the predictors of
% interest. 
% NEW (14 May 2013): The satOffset should only apply if form = 'sat2'. This
% used to not matter because the only time length(b) >= nBeta+2 was when
% form was sat2.  But we since added gamma term, which makes length(b) ==
% nBeta+2 for the 'sat' form.  This too wasn't a problem because the
% command below could still run against length(b) == nBeta+2 when there was
% only one posPOI.  But now that we have a gamma term AND are running
% against Kirby (who has multiple stim parameters), the call below exceeds
% length(b).  So now we limit the call below to when form = sat2.  This
% shouldn't be a problem, but it's a change in an otherwise arcane file, so
% I don't want to screw it up!
if strcmp(form,'sat2') && length(b) >= nBeta+2
    satOffset = b(nBeta+2:nBeta+1+length(posPOI));
end

% the next element is the parameter that scales the saturation point with
% the covariate of interest
if length(b) >= nBeta+2+length(posPOI)
    cov4SatVal = b(nBeta+2+length(posPOI));
end

% last element is the gamma parameter, if included
if ~isempty(posCov4Gamma)
    gamma = b(end);
else 
    % default -- no effect of gamma
    gamma = 1;
end

% reduce b to only the beta coefficients
b = b(1:nBeta);

% % the next element after beta coefficients is the saturation parameter
% if length(b) >= nBeta+1
%     if ismember(form,{'sat','sat2','sat2cov'})
%         satVal = b(nBeta+1);
%     else 
%         error('More parameters provided in B than expected')
%     end
% elseif ismember(form,{'sat','sat2','sat2cov'})
%     error('For saturating forms of logistic, a single saturation value must be passed in B');
% end
% 
% % the next element(s) is the offset parameter(s) for the predictors of
% % interest
% if length(b) >= nBeta+2
%     if ismember(form,{'sat2','sat2cov'})
%         satOffset = b(nBeta+2:nBeta+1+length(posPOI));
%     else 
%         error('More parameters provided in B than expected')
%     end        
% elseif ismember(form,{'sat2','sat2cov'})
%     error('For ''sat2'' and ''sat2cov'' forms of logisitic, the same number of saturation offset parameters must be provided as there are predictors of interest in posPOI');
% end
% 
% % the next element is the parameter that scales the saturation point with
% % the covariate of interest 
% if length(b) >= nBeta+2+length(posPOI)
%     if strcmp(form,'sat2cov')
%         cov4SatVal = b(nBeta+2+length(posPOI));
%     else
%         error('More parameters provided in B than expected')
%     end 
% elseif strcmp(form,'sat2cov')
%     error('For ''sat2cov'' form of logisitic, a coefficient for the covariate that scales the saturation point must be provided');
% end


%% checks -- take too long

% % check that satVal is the right length
% if ismember(form,{'sat','sat2','sat2cov'})
%     if length(satVal) ~= 1
%         error('For saturating forms of logistic, a single saturation value must be passed in SATVAL');
%     end
% end

% % compute number of putative predictors of interest
% nBinary = 0;
% for i = 1:size(X,2)
%     if isequal([0 1]',unique(X(:,i)))
%         nBinary = nBinary + 1;
%     end
% end

% % predictors of interest are mutually exclusive
% if any(sum(X(:,posPOI),2) > 1)
%     error('Predictors of interest must be mutually exclusive')
% end

% if ismember(form,{'sat2','sat2cov'})
%     % predictors of interest are specified.
%     if isempty(posPOI)
%         error('For ''sat2'' and ''sat2cov'' forms of logisitic, at least one predictor of interest must be specified.')
%     end
%     
%     % predictors of interest match number of binary variables
%     if length(posPOI) ~= nBinary
%         warning('Number of predictors of interest does not match the number of putative binary predictors provided in X')
%     end
% 
% %     % saturation offset params match predictors of interest
% %     if length(satOffset) ~= length(posPOI)
% %         error('For ''sat2'' and ''sat2cov'' forms of logisitic, the same number of saturation offset parameters must be provided as there are predictors of interest in posPOI');
% %     end
%         
% end

% if strcmp(form,'sat2cov')
% %     % value for covariate term is provided
% %     if length(cov4SatVal) ~= 1
% %         error('For ''sat2cov'' form of logisitic, a coefficient for the covariate that scales the saturation point must be provided')
% %     end
%     
%     % position of covariate
%     if length(posCov4Sat) ~= 1
%         error('For ''sat2cov'' form of logisitic, the position of the covariate that scales the saturation point must be provided.')
%     end
% end


%% compute logistic

% first raise the covariate of interest to the gamma
X(:,posCov4Gamma) = X(:,posCov4Gamma) .^ gamma;

switch form
    case 'basic'
        % compute numerator (moot in case of 'basic' form)
        numerator = 1;
        
        % alt
        p = (numerator ./ (1+exp(-b*X')))';
        % old way
        % p = (exp(b*X')./(1+exp(b*X')))';
    case 'sat'
        % compute numerator
        numerator = satVal;
        
        % alt
        p = (numerator.* 1./(1+exp(-b*X')))'; 
        % old way
        % p = (b(end).* exp(b(1:end-1)*X')./(1+exp(b(1:end-1)*X')))'; 
        
    case 'sat2'
        % compute numerator
        numerator = satVal + satOffset*X(:,posPOI)';

        % alt
        p = (numerator .* 1./(1+exp(-b*X')))'; 
        % old way
        % p = ((b(end-1) + b(end)*X(:,posPOI)').* exp(b(1:end-2)*X')./(1+exp(b(1:end-2)*X')))'; 

    case 'sat2cov'
        % compute numerator
        numerator = satVal + satOffset*X(:,posPOI)' + cov4SatVal*X(:,posCov4Sat)';
        
        % alt
        p = (numerator .* 1./(1+exp(-b*X')))'; 
        % old way
        % p = ((b(end-2) + b(end-1)*X(:,posPOI)' + b(end)*X(:,posCov4Sat)').* exp(b(1:end-3)*X')./(1+exp(b(1:end-3)*X')))'; 
end

% If modeling lapse rate (1 - numerator) as symmetric, shift the
% logistic upwards by half the lapse rate.
if bitSymmetricLapse
    p = p + (1-numerator)/2;
end
