function [hData, hFit, xFit, yFit]  = CB_psychCurveFitPlotLogistic(xData,yData,varargin);
% [hData,hFit, xFit, yFit] = CB_psychCurveFitPlotLogistic(xData,yData,[b],[X], ...
%       [form],[posPOI],[posCOI],[ci],[nSample],[bitPlot],[p],[nll],[alpha],...
%       ,[posCov4Sat],[POIName],[nullName],[posCov4Gamma]) 
%
% Function to evaluate logistic fit curve for a given set of parameters,
% plot curve (optional), and return fit curve data (optional).
%
% xData = C x 1 vector of x-axes values at which to plot the data in yData.
%         Generally, the values in xData can be thought of as covariates in
%         the logistic regression, as such C is the number of discrete
%         covariates.  If xData is not provided and fit curve is desired,
%         must provide posCOI.  Leave empty to avoid plotting data.
% yData = C x 1 or C x 2 vector of y-values on range [0,1] representing the
%         proportion of successes.  If a binary predictor variable exists
%         (see posPOI) and you wish to plot data separately for the the
%         presence and absence of this predictor, then pass the y-values in
%         columns 1 and 2, respectively, of yData.  Leave empty to avoid
%         plotting data.
%
% %%% The following arguments are required for plotting fit curve.  If the
% %%% fit curve is not required, you may omit all subsequent arguments.
% b = 1 x p vector of parameters estimates from logisticFit().
% X = M x p matrix of predictor values, with each condition given rowwise
%     from 1 to M, and each predictor given column wise from 1 to p.  The
%     first column must be all ones to serve as a constant.  The matrix
%     must be the same as that used by logisticFit() to find b, including
%     any interaction columns.
% form = string specifying form of the logistic link function. Must be the
%        same as that used by logisticFit().  
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
%                B(end). Note that the linear weights are represented in
%                B(1:end-2). 
% posPOI = vector specifying the column position(2) in X of the predictor(s) of
%          interest.  This predictor must be a binary variable.  It is
%          required with certain forms of the logistic (e.g., 'sat2' -- see
%          above) that estimate an additive effect of this predictor.  It
%          is required when plotting the effect of the binary predictor as
%          implied when yData is a C x 2 matrix.
% posCOI = scalar specifying the column position in X of the predictor that
%          serves as a covariate of nSuc/nTrial. This predictor is used as
%          the x-axis when xData is not provided, and thus only required
%          when xData is not provided.
% ci = confidence intervals on parameter estimates b.  Optional for
%      plotting error bars.
% nSample = number of points at which to sample logistic curve. (100).
% bitPlot = logical specifying whether (=true) or not (=false) to actually
%           make plot.  This is useful because the function can also be
%           used to compute the fit curves without actually plotting them.
%           (true).
% p = 1 x p vector of p-values associated with the parameters in b. 
% nll = scalar with negative log likelihood value, which will be printed on
%       plot if provided.
% alpha = scalar of p value below which parameter value should be
%         considered significant. (0.05).
% posCov4Sat = scalar specifying the column position in X of the covariate
%          that will influence the height of the curve saturation as a
%          function of the covariate (e.g., trial number).  Applies only in
%          saturating forms of the logistic function for which special
%          allowance has been made, i.e., 'sat2cov' 
% POIName = cell array with each entry containing the name (description) of
%          a predictor of interest listed in posPOI.  The order of cells
%          should match the order of indecies in posPOI.
% nullName = string providing description of condition in which no
%          predictors of interest are present.
% posCov4Gamma = scalar specifying the column position in X of the
%          covariate that you wish to transform by raising it to the gamma
%          -- a new parameter.  When this field points to a column, the
%          gamma parameter is automatically included and expected to be the
%          last parameter in the optional input vector B.  When the field
%          is empty, the gamma parameter is not included (i.e., set to 1).
%   
% RETURNS:
% hData = handle to data plot. May contain 1 or two handles, depending on
%         whether a binary predictor variable was NOT or was specified,
%         respectively.
% hFit = handle to fit plot. May contain 1 or two handles, depending on
%         whether a binary predictor variable was NOT or was specified,
%         respectively.
% xFit, yFit = N x 1 or N x nStimCondition matrices of the x and y data used to plot the
%              fits, where N is the number of samples specfied in nSample.
%              nStimCond columns are returned when posPOI is given
%
% Daniel Kimmel, January 25, 2010
% Updated 17 Oct 2016 to accommodate multiple unique combinations of the
% predictors of interest (POI).



%% gather incoming vars:

if length(varargin) > 0
    b = varargin{1};
else
    b = [];
end
if length(varargin) > 1
    X = varargin{2};
else
    X = [];
end
if length(varargin) > 2
    form = varargin{3};
else
    form = [];
end
if length(varargin) > 3
    posPOI = varargin{4};
else
    posPOI = [];
end
if length(varargin) > 4
    posCOI = varargin{5};
else
    posCOI = [];
end
if length(varargin) > 5
    ci = varargin{6};
else
    ci = [];
end
if length(varargin) > 6
    nSample = varargin{7};
else
    nSample = 100; 
end
if isempty(nSample)
    nSample = 100;
end
if length(varargin) > 7
    bitPlot = varargin{8};
else
    bitPlot = true; 
end
if isempty(bitPlot)
    bitPlot = true;
end
if length(varargin) > 8
    p = varargin{9};
else
    p = []; 
end
if length(varargin) > 9
    nll = varargin{10};
else
    nll = []; 
end
if length(varargin) > 10
    alpha = varargin{11};
else 
    alpha = [];
end
if isempty(alpha)
    alpha = 0.05; 
end
if length(varargin) > 11
    posCov4Sat = varargin{12};
else
    posCov4Sat = [];
end
if length(varargin) > 12
    POIName = varargin{13};
else
    POIName = {};
end
if isempty(POIName)
    POIName = {};
end
if length(varargin) > 13
    nullName = varargin{14};
else
    nullName = {};
end
if isempty(nullName)
    nullName = {};
end
if length(varargin) > 14
    posCov4Gamma = varargin{15};
else
    posCov4Gamma = [];
end
%% checks

% if xData or yData is provided, then we'll plot the data points:
if ~isempty(xData) || ~isempty(yData)
    bitPlotData = true;
else
    bitPlotData = false;
end

% if data are to be plotted, xData and yData must be properly specified:
if bitPlotData
    if size(xData,1) ~= size(yData,1)
        error('For plotting data, xData and yData must be specified and have the same number of rows.')
    end
end

% if b is provided, then the fit curve is to be plotted:
if ~isempty(b)
    bitPlotFit = true;
else
    bitPlotFit = false;
end

% make sure at least one type of plot is specified:
if ~bitPlotData && ~bitPlotFit
    error('At least one type of plot -- data or fit -- must be specified');
end

% make sure data necessary for fit plot are present:
if bitPlotFit
    if isempty(X) 
        error('When plotting fits, X must be specified.');
    end
    
    if isempty(form)
        error('When plotting fits, form must be specified');
    end
    
    % covariate is required if not implied by xData
    if isempty(xData) && isempty(posCOI)
        error('Covariate specified in posCOI is required when xData is not provided');
    end
    
    % Binary predictor required for certain forms
    if strcmpi(form,'sat2') && isempty(posPOI)
        error('Binary predictor specified in posPOI is required when plotting form %s',form);
    end
    
    % Binary predictor required when yData imply that the presence and
    % absence of the predictor are to be compared.
    if size(yData,2) > 1 && isempty(posPOI)
        error('Binary predictor specified in posPOI required when size(yData,2) > 1 (implying that the presence and absence of the predictor are to be compared.');
    end
        
end

%% default values

% compute number of conditions as defined by the precense or absense of
% the predictors of interest
if ~isempty(posPOI)
%     % assume each entry in posPOI is a new condition, plus 1 for when all
%     % POI = 0.
%     nCond = length(posPOI) + 1;
    % alternatively, assume each unique combination of POI is a new
    % condition
    nCond = size(unique(X(:,posPOI),'rows'),1);
elseif ~isempty(yData)
    nCond = size(yData,2);
% elseif ~isempty(posPOI)
%     foo = unique(X(:,posPOI),'rows');
%     nCond = size(foo,1);
%     % make sure first row is all zeros. if not add it
%     if ~all(foo(1,:) == 0)
%         warning('No condition exists in which all predictors of interest are zero. We attempt to correct for this, but results may not be reliable.')
%         nCond = nCond + 1;
%     end
else
    nCond = 1;
end
    

% various colors
% sigColorStrFull = 'magenta';
textColor = 'w';
POIColorSetName = 'lines';

% we can optionally relabel the x-axes from the normalized units of nRwd
% ([0 1]) to the actual number of rewards.  Do this by multiplying existing
% values by some scalar.  Set the scalar to 1 to leave unchanged.
% Assume that if the range of xData is [0 1], then we need to upscale the
% axis to [0 8], otherwise leave alone
if isequal(minmax(xData'),[0 1])
    xAxisLabelScalar = 8;
else
    xAxisLabelScalar = 1;
end

xFit=[];
yFit=[];

%% prepare vars

hData = [];
hFit = [];

%% remove unspecified coefficients 

% some coefficients are passed that are not relevant for the present data.
% These coeffients are NaN'd out, as are the predictor values to which they
% would correspond.
posNaNCoeff = find(isnan(b));

% also find NaN'd columns in X matrix of predictor values
posNaNX = find(all(isnan(X)));

% expect match between NaN's coeffs and predictors up to the size of X
if (~isempty(posNaNX) || ~isempty(posNaNCoeff))...
        && ~isequal(posNaNX,posNaNCoeff(posNaNCoeff <= size(X,2)))
    error('NaN''d out predictor values do not match NaN''d out coefficient values')
end

% remove these predictors from all relevant variables:
X(:,posNaNX) = [];
b(posNaNCoeff) = [];
ci(:,posNaNCoeff) = [];
p(posNaNCoeff) = [];

% Generally speaking, a continuous predictor will screw up plotting unless
% the predictor is the covariate of interest (COI). Therefore, screen each
% predictor being being continuous and exclude it if it's not specified in
% posCOI.
% The first criterior requires that there be at least 60 values (often the
% min number of trials), which we test first
posCont2Elim = [];
if size(X,1) >= 60
    for i = size(X,2):-1:1
        % determine if continuous var, which we define as having at least 60
        % values (often the min number of trials) and 90% unique values
        if length(unique(X(isfinite(X(:,i)),i)))/sum(isfinite(X(:,i))) >= 0.9
            % check if col was named covariate of interest
            if posCOI ~= i
                % store column
                posCont2Elim(end+1) = i;
%                 warning('Continuous var was passed to logistic plot without being named as covariate of interest and therefore was REMOVED.')
            end
        end
    end
end

%%%% OVERRIDE ABOVE, AND DO NOT ELIMINATE ANY COLUMNS BASED ON CONTINUOUS
%%%% VARS
% posCont2Elim = [];
% 
% if ~isempty(posCont2Elim)
%     % remove columns and all references that depend on them
%     X(:,posCont2Elim) = [];
%     b(:,posCont2Elim) = [];
%     if ~isempty(ci)
%         ci(:,posCont2Elim) = [];
%     end
%     if ~isempty(p)
%         p(posCont2Elim) = [];
%     end
%     
%     % adjust posCOI if necessary
%     if length(posCOI) > 1
%         error('We assume that position pointer has only single element')
%     elseif ismember(posCOI,posCont2Elim)
%         error('Should not eliminate a predictor that was specified as a predictor of interest. Something is not right...');
%     elseif ~isempty(posCOI) && any(posCOI > posCont2Elim)
%         posCOI =  posCOI - sum(posCOI > posCont2Elim);
%     end
%     % adjust posCov4Gamma if necessary
%     if length(posCov4Gamma) > 1
%         error('We assume that position pointer has only single element')
%     elseif ismember(posCov4Gamma,posCont2Elim)
%         error('Should not eliminate a predictor that was specified as a predictor of interest. Something is not right...');
%     elseif ~isempty(posCov4Gamma) && any(posCov4Gamma > posCont2Elim)
%         posCov4Gamma =  posCov4Gamma - sum(posCov4Gamma > posCont2Elim);
%     end
%     % adjust posCov4Sat if necessary
%     if length(posCov4Sat) > 1
%         error('We assume that position pointer has only single element')
%     elseif ismember(posCov4Sat,posCont2Elim)
%         error('Should not eliminate a predictor that was specified as a predictor of interest. Something is not right...');
%     elseif ~isempty(posCov4Sat) && any(posCov4Sat > posCont2Elim)
%         posCov4Sat =  posCov4Sat - sum(posCov4Sat > posCont2Elim);
%     end
%     % adjust posPOI if necessary
%     if length(posPOI) > 1
%         error('We assume that position pointer has only single element')
%     elseif ismember(posPOI,posCont2Elim)
%         error('Should not eliminate a predictor that was specified as a predictor of interest. Something is not right...');
%     elseif ~isempty(posPOI) && any(posPOI > posCont2Elim)
%         posPOI =  posPOI - sum(posPOI > posCont2Elim);
%     end
% end
% 
% also remove predictors from posPOI, but it's tricky because we also need
% to readjust the values in posPOI to account for the now missing columns
pos2Del = [];
posPOINew = posPOI;
for i = 1:length(posNaNX)
    if ~ismember(posNaNX(i),posPOI)
        continue;
    end
    % find position to remove
    pos2Del(end+1) = find(ismember(posPOI,posNaNX(i)));
    
    % reduce all subsequent values by one
    posPOINew(pos2Del(end):end) = posPOINew(pos2Del(end):end) - 1;
end
% remove NaN'd entries
posPOINew(pos2Del) = [];

% % adjust nCond
% nCond = length(posPOINew) + 1;

% update
posPOI = posPOINew;

% % also remove colors that are no longer used. Add 1 to skip the "All POI =
% % 0 row".
% POIColorSet(pos2Del+1,:) = [];

% remove descriptions that are no longer used
if ~isempty(POIName)
    POIName(pos2Del) = [];
end

% compute number of conditions as defined by the precense or absense of
% the predictors of interest
if ~isempty(posPOI)
%     % assume each entry in posPOI is a new condition, plus 1 for when all
%     % POI = 0.
%     nCond = length(posPOI) + 1;
    % alternatively, assume each unique combination of POI is a new
    % condition
    nCond = size(unique(X(:,posPOI),'rows'),1);
elseif ~isempty(yData)
    nCond = size(yData,2);
% elseif ~isempty(posPOI)
%     foo = unique(X(:,posPOI),'rows');
%     nCond = size(foo,1);
%     % make sure first row is all zeros. if not add it
%     if ~all(foo(1,:) == 0)
%         warning('No condition exists in which all predictors of interest are zero. We attempt to correct for this, but results may not be reliable.')
%         nCond = nCond + 1;
%     end
else
    nCond = 1;
end
    
POIColorSet = eval(sprintf('%s(%d)',POIColorSetName, nCond));
POIColorSet(1,:) = [1 1 1]; % make no POI condtion white

% Define set of unique combinations of the predictor(s) of interest
[POISet,~,POIInd] = unique(X(:,posPOI),'rows');

% If names of these predictors were provided, flesh these out for each
% unique set
if ~isempty(POIName)
    POIDesc = cell(size(POISet,1),1);
    for i = 1:size(POISet,1)
        for j = 1:size(POISet,2)
            if j > 1
                foo = ', ';
            else
                foo = '';
            end
               
            POIDesc{i} = sprintf('%s%s%s = %d',POIDesc{i},foo,POIName{j},POISet(i,j));
            clear foo
        end
    end
end

%% plot data points

if bitPlotData && bitPlot

    % plot true pSuc vs COI
    hData = plot(xData,yData,'o','MarkerSize',8);
    for i = 1:size(yData,2)
        set(hData(i),'Color',POIColorSet(i,:),'MarkerFaceColor',POIColorSet(i,:));
    end
end


%% plot fits
if bitPlotFit 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % build predictor matrix that has entries for each parameter. However,
    % since only the covariate of interest (COI) varies continuously in the
    % plot, we have to compute the fitted value given the expected value
    % for all non-COI parameters. As a first approximation, we do this by
    % taking the average of all non-COI parameters for a given level of the
    % COI parameter. This is essentially the CONDITIONAL MARGINALIZED value
    % of the other non-COI predictors.
    % Note that we do the above procedure separately for each level of the
    % predictor of interest (POI).
    
    % initialize
    yFit = NaN(nSample,size(POISet,1));
    xFit = yFit;
    
    % loop through levels of the POI
    for i = 1:size(POISet,1)
        % temp predictors
%         XTemp = X;
        
        % find unique values of COI
        % Note that some values of COI may only be present with a given
        % value of POI, so COI may not be present in the data for all
        % levels of POI. Just tolerate this for now.

        if ~isempty(posCOI)
            COISet = unique(X(POIInd==i,posCOI));
        else
            % infer from XData
            COISet = unique(xData);
        end
        
        % remove POI from predictors.         
%         XTemp(:,posPOI) = [];

        % initialize holding var for predictors used for fit. Note that
        % this includes POI column
        XFitData = NaN(length(COISet),size(X,2));
        
        % build mask as to where to assign non-POI vars
        bitNonPOI = true(1,size(X,2));
        bitNonPOI(posPOI) = false;
        
        % loop through unique values of COI
        for j = 1:length(COISet)
            % take average of non-POI vars for each level of COI
            XFitData(j,bitNonPOI) = nanmean(X(POIInd==i & X(:,posCOI)==COISet(j),bitNonPOI),1);
        end
        
        % add back the POI
        XFitData(:,~bitNonPOI) = repmat(POISet(i,:),size(XFitData,1),1);
        
        % compute y for just the data points:
        yFitData = logistic(b,XFitData,form,'posPOI',posPOI,'posCov4Sat',posCov4Sat,'posCov4Gamma',posCov4Gamma);
        
        % plot fitted data points
        if bitPlot
            hold on;
            hFit = plot(COISet,yFitData,'x','MarkerSize',10);
            set(hFit,'Color',POIColorSet(i,:));
        end
        
        clear XFitData yFitData
        
        %%%%%%%%%%%%%%%%
        % Now we plot a smooth line. Since non-POI and non-COI predictors
        % are not defined for the in-between points on the COI (x-axis), we
        % use the MARGINALIZED value of the non-COI predictors (averaged
        % across all values of COI) for a given POI.
        
        % establish x-axis.
        XFit = NaN(nSample,size(X,2));
        xFit(:,i) = linspace(min(COISet),max(COISet),nSample)';
        XFit(:,posCOI) = xFit(:,i);
            
        % build mask for nonCOI columns
        bitNonCOI = true(1,size(X,2));
        bitNonCOI(posCOI) = false;
        
        XFit(:,bitNonPOI & bitNonCOI) = repmat(nanmean(X(POIInd==i,bitNonPOI & bitNonCOI),1),...
            nSample,1);

        % add back the POI
        XFit(:,~bitNonPOI) = repmat(POISet(i,:),nSample,1);
        
        % compute y for a smooth line:
        yFit(:,i) = logistic(b,XFit,form,'posPOI',posPOI,'posCov4Sat',posCov4Sat,'posCov4Gamma',posCov4Gamma);
        
        % plot fitted data points
        if bitPlot
            hold on;
            hFit = plot(xFit(:,i),yFit(:,i),'-','LineWidth',2);
            set(hFit,'Color',POIColorSet(i,:));
        end
        
        clear XFit yFitData 


    end
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%
%     % OLD WAY of handing non-plotted predictors. It did not marginalize
%     % and used some combination of eliminating continuous predictors and
%     % using continuous values for non-plotted predictors.
%     %%%%%%%%%%%%%%%%%%%%%%%%%
%     % generate fit curves
%     % make matrix of predictor values spanning range of input predictors:
%     XFit = NaN(nSample,size(X,2));
%     % if POI given, isolate that those conditions when all predictors of
%     % interest == 0
%     if ~isempty(posPOI)
%         bitFoo = all(X(:,posPOI) == 0,2);
%     else
%         bitFoo = true(size(X,1),1);
%     end
%     
%     % loop through all predictors and build vectors of the independent
%     % variable spanning the range of each predictor
%     for i = 1:size(X,2)
%         XFit(:,i) = linspace(min(X(bitFoo,i)),max(X(bitFoo,i)),nSample)';
%     end
%     clear bitFoo
%     % establish x-axis. 
%     if ~isempty(posCOI)
%         % use the covariate indicated in posCOI
%         xFit = XFit(:,posCOI);
%     else
%         % infer range from xData
%         xFit = linspace(min(xData),max(xData),nSample)';
%     end
%     % compute y:
%     yFit = logistic(b,XFit,form,'posPOI',posPOI,'posCov4Sat',posCov4Sat,'posCov4Gamma',posCov4Gamma);
%     
%     % plot first fit:
%     if bitPlot
%         hold on;
%         hFit = plot(xFit,yFit,'-','LineWidth',2);
%         set(hFit,'Color',POIColorSet(1,:));
%     end
%     clear XFit
%     
%     % if POI given, isolate each condition where a predictor of interest == 1 
%     if ~isempty(posPOI)
%         
%         for j = 1:length(posPOI)
% 
%             bitFoo = X(:,posPOI(j)) == 1;
%             % make matrix of predictor values spanning range of input predictors:
%             XFit = NaN(nSample,size(X,2));
%             % loop through predictors
%             for i = 1:size(X,2)
%                 XFit(:,i) = linspace(min(X(bitFoo,i)),max(X(bitFoo,i)),nSample)';
%             end
%             % establish x-axis. 
%             if ~isempty(posCOI)
%                 % use the covariate indicated in posCOI
%                 xFit(:,1+j) = XFit(:,posCOI);
%             else
%                 % infer range from xData
%                 xFit(:,1+j) = linspace(min(xData),max(xData),nSample)';
%             end    
% 
%             % compute y
%             yFit(:,1+j) = logistic(b,XFit,form,'posPOI',posPOI,'posCov4Sat',posCov4Sat,'posCov4Gamma',posCov4Gamma);
% 
%             % plot additional fit:
%             if bitPlot
%                 hFit(1+j) = plot(xFit(:,1+j),yFit(:,1+j),'-','LineWidth',2,'Color',POIColorSet(1+j,:));
%             end
%             clear bitFoo
%         end
%         clear XFit
%     end
end

%%%%% RETURN if not plotting
if ~bitPlot
    return
end

%%%%% OTHERWISE continue...

% set limits
if ~isempty(xData)
    % based on data
    xInfoTemp = xData;
    set(gca,'XTick',unique(xData'));
else
    % based on fits
    xInfoTemp = X(:,posCOI);
    xlim([min(X(:,posCOI))-0.125*range(X(:,posCOI)) max(X(:,posCOI))]);
end
xlim([min(xInfoTemp)-0.125*range(xInfoTemp) max(xInfoTemp)]);
ylim([0 1]);
% compute x position for plotting parameter value
xOffsetPlot = 0.0250*range(xInfoTemp);
xParam = min(xInfoTemp) - 0*xOffsetPlot;
clear xInfoTemp

%%%%%%%%%%%%%%%%%%
% plot parameters

% plot saturation point
if strncmpi(form,'sat',3)
    posParam = size(X,2)+1;
    % prepare errorbars:
    if ~isempty(ci)
        lb = b(posParam) - ci(1,posParam);
        ub = ci(2,posParam) - b(posParam);
        h = errorbar(xParam,b(posParam),lb,ub,'o','Color',POIColorSet(1,:));
    else
        h = plot(xParam,b(posParam),'o','Color',POIColorSet(1,:));
    end
    set(h,'MarkerSize',5)
    % fill in if significant
    if p(posParam) <= alpha
        set(h,'MarkerFaceColor',POIColorSet(1,:));
    end

    % plot POI offset on saturation point
    if strcmpi(form,'sat2')
        hold on;
        % loop through POI conditions
        for j = 2:nCond
            % prepare parameter
            bOffset = b(posParam + j-1);
            bAlt = b(posParam) + bOffset;
            if ~isempty(ci)
                lb = bOffset - ci(1,posParam + j-1);
                ub = ci(2,posParam + j-1) - bOffset;
                h = errorbar(xParam - xOffsetPlot*(j-1),bAlt,lb,ub,'o','Color',POIColorSet(j,:));
            else
                h = plot(xParam - xOffsetPlot*(j-1),bAlt,'o','Color',POIColorSet(j,:));
            end
            set(h,'MarkerSize',5);

            % fill in if significant
            if p(posParam + j-1) <= alpha
                set(h,'MarkerFaceColor',POIColorSet(j,:));
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WRITE TEXT FOR OTHER PARAMS

% % ignore saturation parameter, if present
% posParam = 1:length(b);
% if strncmpi(form,'sat',3)
%     posParam(size(X,2)+1) = [];
%     if strcmpi(form,'sat2')
%         % offset param is now just one away from predictor list, since the
%         % saturation param was just removed.
%         posParam(size(X,2)+1) = [];
%     end
% end

% printed parameters should just be the standard beta coefficients, which
% are as numerous as columns in X:
posParam = 1:size(X,2);

% prepare param text. 
paramText = [''];
for i = posParam
    paramText = [paramText,'{'];
    
    if ~isempty(p) && p(i) <= alpha
        % italicize significant values
        paramText = [paramText,'\it\fontsize{10}'];
        
%         % color significant values
%         paramText = [paramText,'\color{',sigColorStrFull,'}'];
    end
    
    % if parameter is a predictor of interest, color it by that predictor.
    if ismember(i,posPOI) %ismember(i,posPOI)
        % param matches a predictor of interest either if it's in the same
        % position as the predictor (in posPOI) or if it's in that same
        % position but offset by the number of predictors of interest
        % (and are therefore interactions with those predictors).
        if i <= max(posPOI)
            foo = find(posPOI==i) + 1; % add 1 to skip the no-POI color.
        else
            foo = find(posPOI+length(posPOI)==i) + 1;
        end
        % better to print without a color than request an [] color causing
        % an error
        if ~isempty(foo)
            paramText = [paramText,'\color[rgb]{',num2str(POIColorSet(foo,:)),'}'];
        end
    end
    
    paramText = [paramText,num2str(b(i),2),'}'];
    if i ~= posParam(end)
        paramText = [paramText,', '];
    end
end

% add text for gamma parameter
if ~isempty(posCov4Gamma)
    if ~isempty(paramText)
        paramText = [paramText,', '];
    end
    if ~isempty(p) && p(end) <= alpha
        % italicize significant values
        paramText = [paramText,'\it\fontsize{10}'];
    end
    paramText = [paramText, '\gamma = ',num2str(b(end),2)];
end

% write param text
xLim = get(gca,'XLim');
yLim = get(gca,'YLim');
hText = text(xLim(2),yLim(1),paramText);
set(hText,'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',9,'Color',textColor);

% write nll text
if ~isempty(nll)
    hText = text(0,yLim(2),['NLL=',sprintf('%0.1f',nll)]);
    set(hText,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',9,'Color',textColor);
end    

% relabel xaxis with actual (not normalized) number of rewards.
if xAxisLabelScalar ~= 1
    foo = get(gca,'XTick');
    set(gca,'XTickLabel',foo*xAxisLabelScalar);
end

% title figure only if a predictor of interest
if ~isempty(posPOI)
    if ~isempty(POIName) && length(POIName) ~= length(posPOI)
        error('Descriptions for predictors of interest must match the number of predictors of interst')
    end
    
    % loop through levels of POISet
    titleText = cell(size(POISet,1),1);
    for i = 1:size(POISet,1)
        titleText{i} = ['\color[rgb]{',num2str(POIColorSet(i,:)),'}',POIDesc{i}];
%         if i==1 && ~isempty(nullName)
%             titleText{i} = [titleText{i},nullName];
%         else
%             titleText{i} = [titleText{i},'No Manipulation'];
%         end
        
    end
    
%     titleText{1,1} = ['\color[rgb]{',num2str(POIColorSet(1,:)),'}'];
%     if ~isempty(nullName)
%         titleText{1,1} = [titleText{1,1},nullName];
%     else
%         titleText{1,1} = [titleText{1,1},'No Manipulation'];
%     end
%     
%     for i = 1:length(posPOI)
%         titleText{1+i,1} = ['{\color[rgb]{',num2str(POIColorSet(i+1,:)),'}'];
%         if ~isempty(POIName)
%             titleText{1+i,1} = [titleText{1+i,1},POIName{i}];
%         else
%             titleText{1+i,1} = [titleText{1+i,1},'Predictor ',num2str(i)];
%         end
%         titleText{1+i,1} = [titleText{1+i,1},'}'];
%     end
    title(titleText);
end

% turn off box
set(gca,'box','off');

return


