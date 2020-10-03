
function [pVal,prms,err] = sigTest(dataStatistic, nullStatistic, tail, option,varargin)
% [pVal,prms,err] = sigTest(dataStatistic, nullStatistic, tail, option,[name,value])
%
% Compares dataStatistic to null distribution of values given in
% nullStatistic and returns probability pVal of observing dataStatistic.
%
% ACCEPTS
% dataStatistic -- 1 x N vector of data values for N samples
% nullStatistic -- S x N matrix of null (or surrogate) values for N
%       datasets corresponding to the N samples in dataStatistic, where
%       each dataset contains S surrogate samples. When
%       nullStatistic is S x 1, it is assumed that the single set of S
%       surrogate samples is to apply to all N data samples.
% tail -- string specifying whether pVal corresponds to the integral of the
%       null distribution from -Inf to dataStatistic ('lower') or
%       dataStatistic to Inf ('upper'), or a 2-tailed test ('both'), in
%       which case the smaller p-value from 'upper' and 'lower' is doubled
%       and returned.
% option -- string specifying the type of distribution to use. 
%       'empirical' -- dataStatistic is compared directly to nullStatistic,
%           e.g., pVal = sum(nullStatistic < dataStatistic) ./
%           length(nullStatistic)
%       'gamma' -- values in nullStatistic are fit to gamma distribution,
%           with gamma parameters shape and scale returned in PRMS(1) and
%           PRMS(2), respectively. Values in nullStatistic must therefore
%           be greater than or equal to 0. When values = 0, the MLE fit is
%           not possible and an approximation is used.
%       'gauss' OR 'normal' -- values in nullStatistic are fit to a
%           Gaussian or normal distribution, with parameters mu and sigma
%           returned in PRMS(1) and PRMS(2), respectively.
%       'halfguass' OR 'half normal' OR 'hn' -- values in nullStatistic are
%           fit to a half-Gaussian or half-normal, where if X is normally
%           distributed, then Y = abs(X) is distributed half-normal. 
%           NOTE:
%           Parameters mu and sigma returned in PRMS(1) and PRMS(2),
%           respectively, refer to distribution Z, which is normally
%           distributed and constructed by concatenating Y and -Y. Thus Y
%           is given by 2 * normpdf([0:Inf],mu,sigma) and the integral of Y
%           is 2 * normcdf(mu,sigma) evaluated from 0 to Inf.
%           NOTE 2:
%           In Matlab versions 9.0 (2016a) and later, half normal
%           distributions are fit using the built-in 'half normal'
%           distribution (see FITDIFF, PDF, CDF, etc) and the resulting
%           'upper' tail p-values are about 0.1% smaller than using a
%           home-spun alternative. The alternative reflects the input data
%           DATASTATISTIC about zero (data = [-data,data]) and fits a
%           normal distribution. 
% 
% OPTIONALLY ACCEPTS as name,value pairs
% 'pOORMax' -- scalar specifying the proportion of values in any column of
%       nullStatistic that are permitted to be out of range. For instance,
%       for the gamma fit, values must be non-negative. If negative values
%       are present but the proportion of illicit values is < pOORMax,
%       these values will be removed. Otherwise, error will be throw.
%       Default pOORMax = 0.001.
% 'prms' -- N x P matrix of P model-fit parameter values for each of N
%       samples. When provided, these values will be used in place of
%       fitting the model to the data de novo. When empty (default), model
%       is fit to the data and parameters returned in PRMS.
% 'errThresh' -- scalar specifying the threshold for sum of squared error
%       between model-based PDF and data frequency (normalized by sum of
%       squared empirical PDF) over which function will return error for
%       poor fit and plot distribution with fit. Generally, errThresh =
%       0.12 is a reasonable value (default). Set to Inf for no
%       threshold. To apply a threshold (and return NaNs when error exceeds
%       threshold) but avoid throwing an error, set bitIgnorePoorFit.
% 'bitIgnore' -- 1 x N logical vector, which when TRUE, will skip the
%       specified observation.
% 'bitIgnorePoorFit' = logical specifying whether to ignore (==TRUE)
%       observations for which the null distribution model was a poor fit
%       (based on errThresh), which case PVAL and PRMS associated with the
%       observation are set to NaN, or to throw an error (==FALSE,
%       default).
%
% RETUNRS
% pVal -- 1 x N vector of p-values corresponding to the values in
%       dataStatistic
% prms -- N x P matrix of P parameter values generating the PDF fit to each
%       column N in nullStatistic. 
% err -- 1 x N vector of mean absolute error between model-based PDF and
%       data frequency.
%
% By Daniel Kimmel and Gamal Elsayed, Last updated 2017 Jan 24

%% default param values

% proportion of nullStatistic value that can be out of range
pOORMax = 0.001;
% fit parameters
prms = [];
% threshold for sum of squared error between model-based PDF and data
% frequency (normalized by sum of squared empirical PDF) over which
% function will return error for poor fit and plot distribution with fit.
% Generally, errThresh = 0.115 is a reasonable value (default). Set to [] or
% Inf for no threshold.
errThresh = 0.12;
bitIgnore = false(size(dataStatistic));
bitIgnorePoorFit = false;

%% collect optionally provided parameters
warnopts(assignopts(who, varargin));

% if input value is empty, replace -- this is to cover older cases in which
% the function allowed for an empty threshold. To have no threshold, set to
% Inf. To apply a threshold (and return NaNs for poor fits) but avoid
% throwing errors, set bitIgnorePoorFit.
if isempty(errThresh)
    errThresh = 0.12;
end

%%
if sum(size(dataStatistic) > 1) > 1
    error('Can only provide a vector of data values')
end
if size(nullStatistic,2)~= length(dataStatistic)
    nullStatistic = nullStatistic.';
end
% determine datasets of nullStatistic 
if sum(size(nullStatistic) > 1) > 1
    nNull = size(nullStatistic,2);
else
    nNull = 1;
end

if size(nullStatistic,2) ~= length(dataStatistic) && nNull > 1
    error('Must provide set of nullStatistic samples for each value of data')
end

% convert option to recognized name
switch option
    case {'hn' , 'half normal'}
        option = 'halfgauss';
    case 'normal'
        option = 'gauss';
end     

% instatiate
pVal = nan(size(dataStatistic));
err = nan(size(dataStatistic));

switch option
    case 'gamma'
        P = 2; % number of parameters
    case 'gauss'
        P = 2; % number of parameters
    case 'halfguass'
        P = 2; % number of parameters
    otherwise
        P = 2; % number of parameters
end

if ~isempty(prms) && (size(prms,1) ~= size(dataStatistic,2) || ...
        size(prms,2) ~= P)
    error('Provided parameter matrix PRMS does not match expected size based on input samples and number of model parameters');
end

% instantiate PRMS matrix if not provided
if isempty(prms)
    prms = NaN(nNull,2);
    bitFitDist = true;
else 
    bitFitDist = false;
end


% numSamples = size(nullStatistic, 1);
for i = 1:length(dataStatistic)
    clear pValU pValL
    
    % skip ignored observations
    if bitIgnore(i)
        continue
    end
    
    % extract null statistic
    if nNull > 1
        temp = nullStatistic(:,i);
        prmI = i;
    else
        % common null statistic for all values of data
        temp = nullStatistic;
        prmI = 1;
    end
    
    if strcmpi(option, 'empirical')
        % empirical distripution based pVal value. 
        % constrain min p-value to 1/length(nullStatistic)
        if strcmpi(tail, 'upper') || strcmpi(tail, 'both')
            pValU = max(1,sum(temp>dataStatistic(i)))/sum(isfinite(temp));
        end
        if strcmpi(tail, 'lower') || strcmpi(tail, 'both')
            pValL = max(1,sum(temp<=dataStatistic(i)))/sum(isfinite(temp));
        end

    
    elseif strcmpi(option, 'gamma')
        % Handle out of range values
        bitFoo = temp < 0;
        foo = sum(bitFoo);
        if foo/sum(isfinite(temp)) > pOORMax 
            error('Proportion of negative null statistic values (%g) is greater than allowed (%g)',foo/sum(isfinite(temp)),pOORMax);
        elseif foo > 0
%             warning('Converting %d/%d = %g proportion of null statistic values to NaN because out of range for %s fit',foo,sum(isfinite(temp)),...
%                 foo/sum(isfinite(temp)),option);
            temp(bitFoo) = [];
        end
        
%         % additionally exclude zero values
%         bitFoo = temp == 0;
%         foo = sum(bitFoo);
%         if foo/sum(isfinite(temp)) > pOORMax * 2
%             error('Proportion of zero-valued null statistic values (%g) is greater than allowed (%g)',foo/sum(isfinite(temp)),pOORMax);
%         elseif foo > 0
% %             warning('Converting %d/%d = %g proportion of null statistic values to NaN because out of range for %s fit',foo,sum(isfinite(temp)),...
% %                 foo/sum(isfinite(temp)),option);
%             temp(bitFoo) = [];
%         end
% 

        % fit gamma to empirical null samples (if positive only)
        % suppress warning message regarding having zero values and thus
        % unable to do MLE 
        if bitFitDist
            warning('OFF', 'stats:gamfit:ZerosInData');
            [prms(prmI,:)] = gamfit(temp);
            warning('ON', 'stats:gamfit:ZerosInData');
        end
        
        clear pValU pValL
        if strcmpi(tail, 'upper') || strcmpi(tail, 'both')
            pValU = 1-gamcdf(dataStatistic(i),prms(prmI,1), prms(prmI,2));
        end
        if strcmpi(tail, 'lower') || strcmpi(tail, 'both')
            pValL = gamcdf(dataStatistic(i),prms(prmI,1), prms(prmI,2));
        end
        clear  foo bitFoo
    elseif strcmpi(option, 'gauss') 
        % fit gaussian to empirical null samples (if positive and negative only)
%         mu = mean(temp);
%         sigma = std(temp);
        if bitFitDist
            [prms(prmI,1),prms(prmI,2)] = normfit(temp);
        end
        if strcmpi(tail, 'upper') || strcmpi(tail, 'both')
            pValU = 1- normcdf(dataStatistic(i),prms(prmI,1),prms(prmI,2));
        end
        if strcmpi(tail, 'lower') || strcmpi(tail, 'both')
            pValL = normcdf(dataStatistic(i),prms(prmI,1),prms(prmI,2));
        end
    elseif strcmpi(option, 'halfgauss') 
        % faster to run as dedicated distribution if available
        if ~verLessThan('matlab','9.0')
            pd = fitdist(temp,'hn');
            [prms(prmI,1:2)] = pd.ParameterValues;
            if strcmpi(tail, 'upper') || strcmpi(tail, 'both')
                pValU = cdf(pd,dataStatistic(i),'upper');
            end
            if strcmpi(tail, 'lower') || strcmpi(tail, 'both')
                pValL = cdf(pd,dataStatistic(i));
            end
            clear pd
        else
            % fit half-gaussian to empirical null samples
            % "unfold" distribution
            temp = [temp;-temp];
            if bitFitDist
                [prms(prmI,1),prms(prmI,2)] = normfit(temp);
            end
            % first remove 0.5, or integral of distribution from -inf to
            % mu. Then normalize remaining integral mu to +inf by 0.5.
            if strcmpi(tail, 'upper') || strcmpi(tail, 'both')
                pValU = 1 - (normcdf(dataStatistic(i),prms(prmI,1),prms(prmI,2)) - 0.5)/0.5;
            end
            if strcmpi(tail, 'lower') || strcmpi(tail, 'both')
                pValL = (normcdf(dataStatistic(i),prms(prmI,1),prms(prmI,2))-0.5)/0.5;
            end
        end
    else
        error('Null distribution name %s not recognized',option);
        
    end
    
    % store p-value and compute 2-tailed p-value if necessary
    if strcmpi(tail, 'upper')
        pVal(i) = pValU;
    elseif strcmpi(tail, 'lower')
        pVal(i) = pValL;
    elseif strcmpi(tail, 'both')
        pVal(i) = 2*min(pValU,pValL);
    else
        error('Tail parameter ''%s'' not recognized',tail);
    end
    clear pValU pValL
    
    
    %%% QUALITY CONTROL
    if ~strcmpi(option, 'empirical') % && ~isempty(errThresh) && isfinite(errThresh)
        % histogram
        [freq,bin] = hist(temp,30);
        % convert model name
        switch option
            case 'gauss'
                optTemp = 'normal';
            case 'halfgauss'
                if ~verLessThan('matlab','9.0')
                    optTemp = 'hn';
                else
                    optTemp = 'normal';
                end
            otherwise
                optTemp = option;
        end
                
        % pdf
        foo = num2cell(prms(prmI,:));
        y = pdf(optTemp,bin,foo{:});
        clear foo
        
        % error
%         err(i) = nanmean(abs((freq / sum(freq) / diff(bin(1:2))) - y)) / ...
%             nanmean(freq / sum(freq) / diff(bin(1:2)));
        err(i) = sum(((freq / sum(freq) / diff(bin(1:2))) - y).^2) / ...
            sum((freq / sum(freq) / diff(bin(1:2))).^2);
        
        if ~isempty(errThresh) && isfinite(errThresh) && err(i) >= errThresh
            
            % choose to suppress error and skip value, or throw error 
            if bitIgnorePoorFit
                pVal(i) = NaN;
                prms(i,:) = NaN;
            else
                % plot histogram
                figure;
                plot(bin,(freq / sum(freq) / diff(bin(1:2))),'Color',[0.5 0.5 0.5]);
                hold on
                plot(bin,y,'k');
                legend({'data','model'});
                error('Error for sample number %d = %0.3g, which was above error threshold %0.3g',...
                    i,err(i),errThresh);
            end            
        end
    end
    clear temp freq bin y 
end



end