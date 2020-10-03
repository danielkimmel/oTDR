function [h,hL,bitObsPlotted,P,hT] = plotScatterN(x,y,varargin)
% [h,hL,bitObsPlotted,P,hT] = plotScatterN(x,y, ... [name,value pairs]) 
% 
% Plots scatter of X and Y in current axes. Sets limits on both axes to be
% equal. Draws unity line, colored magenta if X-Y is significantly
% different from 0 at P < ALPHA by Wilcoxon signed rank test SIGNRANK()
% (alpha must be provided or else no test is done).
%
% Optionally sets custom datatips with strings passed in cell array
% DATATIPLABEL, which must have same number of rows as X and Y. If multiple
% columns are provided, then datatips are customized for each point,
% otherwise if one column is provided, datatips are specific for all points
% in a row of X and Y.
%
% plotScatterN supercedes plotScatter and allows for plotting muliple
% vectors with different markers.  X and Y must have the same number of
% columns.  
%
% ONCLICKFUNCTION optionally contains a function handle to a callback
% function to be called when the cluser mouse-clicks on a point in the
% scatter plot.  If the callback function requires data to passed to it,
% make ONCLICKFUNCTION a cell array, the first cell contains the function
% handle and each subsequent cell contains 1 variable to be passed to the
% callback function.  The cell array must have the same number of rows and
% columns as X and Y.
%
% COLOR (full name or rgb values).  If mutliple N columns are provided in X
% and Y, color must be either a 1 x N or 1x1 cell array of color names, or
% an N x 3 matrix or a 1x3 vector of RGB values.  If only a single color is
% provided, it will be used for all columns in X and Y.
%
% BITSIG -- logical same size and X and Y.
%
% bitEqualAxes -- axes limits and unity line
%
% alpha = significance level at which to highlight statistical test.
% REQUIRED for any test (comparison of medians or correlation) to be
% performed.
% 
% outlierZThresh = threshold for outliers in z-scores, beyond which data 
% points will be excluded from the computation of bin edges. (default =
% Inf).
% 
% statTest = string to specify the type of statistical test to be
% performed. Note that alpha must be set for any test to be performed.
% Possible values: 
% 'signrank' -- for Wilcoxon sign rank test to compare medians (Default)
% 'corrPear' -- for Pearson's linear correlation
% 'corrSpear' -- for Spearman's rank-based regression
% 'ttest' -- for ttest to compare means
% 'totregress' -- for linear "total" regression calling totregress()
%
% bitStatOnSig = logical to perform statistical test only on the data
% marked significant by the bitSig vector. Defaults to FALSE.
% 
% pTextFontSize = font size for probability text
%
% Returns H and HL, handles to the scatter and unity line, respectively.
% Returns bitObsPlotted, a logical vector with the same number of rows as X
% and Y that specifies whether that observation was plotted.
% Also return P, the probability that the median difference is zero.
% Also return hT, handle to the probability text.
%
% Daniel Kimmel, 9 Dec 2009

%% default values:
alpha = [];
datatipLabel = {};
onClickFunction = [];
onClickData = {};
color = [];
bitSig = [];
bitEqualAxes = true;
outlierZThresh = Inf;
statTest = 'signrank';
bitStatOnSig = false;
pTextFontSize = 7; % font size for probability text
%% collect optional name,value pairs
warnopts(assignopts(who, varargin));

%% instantiate:
hL = [];
hT = [];

%% prepare and plot

% if X and Y (and other inputs) are vectors, make sure they are column
% vectors 
if min(size(x)) == 1 && size(x,1) < size(x,2)
    x = x';
end
if min(size(y)) == 1 && size(y,1) < size(y,2)
    y = y';
end
if min(size(datatipLabel)) == 1 && size(datatipLabel,1) < size(datatipLabel,2)
    datatipLabel = datatipLabel';
end
if min(size(color)) == 1 && ~ischar(color) && size(color,1) < size(color,2)
    color = color';
end
if min(size(bitSig)) == 1 && size(bitSig,1) < size(bitSig,2)
    bitSig = bitSig';
end
if min(size(onClickData)) == 1 && size(onClickData,1) < size(onClickData,2)
    onClickData = onClickData';
end

% save axes overwrite property
nextPlot = get(gca,'NextPlot');
hold on;

% handle for attaching custom datatips:
hdt = datacursormode;
set(hdt,'UpdateFcn',{@scatterDatatipCallback});

% instatiate:
bitObsPlotted = false(size(x));
obsind = cell(1,size(x,2));
h = cell(1,size(x,2));
hL = cell(1,size(x,2));
P = NaN(1,size(x,2));
rho = NaN(1,size(x,2));
lineParam = NaN(size(x,2),2);
pTxt = cell(size(x,2),1);
xyTempLim = NaN(size(x,2),2);
xTempLim = NaN(size(x,2),2);
yTempLim = NaN(size(x,2),2);

% prepare color
if iscell(color) && ~isempty(color)
    % repeat color for all columns
    if length(color) < size(x,2)
        foo = color;
        clear color;
        color(1:size(x,2)) = foo(1);
        clear foo
    end
elseif ischar(color)
    foo = color;
    clear color
    % repeat color for all columns
    color(1:size(x,2)) = {foo};
    clear foo
elseif isempty(color)
    goo = get(gca,'ColorOrder');
    color = goo(1:size(x,2),:);
    clear goo
end
    

% X and Y can be matrices with multiple columns to be plotted against one
% another columnwise.
for i = 1:size(x,2)
        
    % prepare x and y data, eliminate nans
    foo = isnan(x(:,i)) | isnan(y(:,i));
    xTemp = x(~foo,i);
    yTemp = y(~foo,i);
    if ~isempty(datatipLabel)
        % reuse datatip labels across columns unless column-specific labels
        % provided:
        if size(datatipLabel,2) > 1
            dtLabelTemp = datatipLabel(~foo,i);
        else
            dtLabelTemp = datatipLabel(~foo);
        end
    else
        dtLabelTemp = {};
    end
    if ~isempty(onClickData)
        onClickDataTemp = onClickData(~foo,i);
    else
        onClickDataTemp = {};
    end
    obsind{i} = 1:length(xTemp);
    bitObsPlotted(:,i) = ~foo;

    if ~isempty(bitSig)
        bitSigTemp = bitSig(~foo,i);
    else
        bitSigTemp = [];
    end
    
    clear foo;

    % if no data, continue
    if isempty(xTemp)
        continue
    end
    
    % prepare color string
    if iscell(color)
        colorTemp = color{i};
    else
        % RGB vector
        colorTemp = color(i,:);
    end
    % convert color abbreviation to full name
    if ischar(colorTemp) && length(colorTemp)==1
        colorTemp = colorConvert(colorTemp);
    end
    
    % find outliers in X and Y, set them to max non-outlier values, and flag
    % them for special plotting:
    xTempForStats = xTemp;
    yTempForStats = yTemp;
    bitOutlierXMax = zscore(xTemp) >= outlierZThresh;
    bitOutlierXMin = zscore(xTemp) <= -outlierZThresh;
    bitOutlierYMax = zscore(yTemp) >= outlierZThresh;
    bitOutlierYMin = zscore(yTemp) <= -outlierZThresh;
    % set outliers to max non outlier value:
    xTemp(bitOutlierXMax) = max(xTemp(~bitOutlierXMax & ~bitOutlierXMin));
    xTemp(bitOutlierXMin) = min(xTemp(~bitOutlierXMax & ~bitOutlierXMin));
    yTemp(bitOutlierYMax) = max(yTemp(~bitOutlierYMax & ~bitOutlierYMin));
    yTemp(bitOutlierYMin) = min(yTemp(~bitOutlierYMax & ~bitOutlierYMin));

    % scatter and hist
    % if an onClick callback function is specified, make each datapoint it's
    % own plot object with unique call to the callback function.  Handles for
    % each data point are compiled in H as a vector.
    if ~isempty(onClickFunction) || ~isempty(bitSigTemp) ...
            || any(any(bitOutlierXMax)) || any(any(bitOutlierYMax)) ...
            || any(any(bitOutlierXMin)) || any(any(bitOutlierYMin))
        for j = 1:length(xTemp)
            
            h{i}(j) = plot(xTemp(j),yTemp(j),'.');
            set(h{i}(j),'Color',colorTemp);

            % make markers circles if significance is to be plotted
            if ~isempty(bitSigTemp)
                set(h{i}(j),'Marker','o');
                
                % color in significant values
                if bitSigTemp(j)
                    set(h{i}(j),'MarkerFaceColor',colorTemp);    
                end
            end
            
            % make markers square if outlier:
            if bitOutlierXMax(j) || bitOutlierYMax(j) || ...
                    bitOutlierXMin(j) || bitOutlierYMin(j)
                set(h{i}(j),'Marker','s','MarkerSize',4);
            end
            
            if ~isempty(onClickFunction) 
                % Attach custom datatip to specific point
                setappdata(h{i}(j),'obsind',obsind{i}(j));
                if ~isempty(dtLabelTemp)
                    setappdata(h{i}(j),'datatipLabel',dtLabelTemp{j});
                end

                % if data is included to be passed to the callback function...
                if ~isempty(onClickDataTemp)
                    set(h{i}(j),'ButtonDownFcn',{onClickFunction,onClickDataTemp{j}});
                else
                    set(h{i}(j),'ButtonDownFcn',onClickFunction);
                end
            else
                % Attach custom datatip to specific point
                setappdata(h{i}(j),'obsind',obsind{i}(j));
                if ~isempty(dtLabelTemp)
                    setappdata(h{i}(j),'datatipLabel',dtLabelTemp{j});
                end
            end
        end
    else
        h{i} = plot(xTemp,yTemp,'.');
        set(h{i},'Color',colorTemp);

        % Attach custom datatips to set of points
        setappdata(h{i},'obsind',obsind{i});
        if ~isempty(dtLabelTemp)
            setappdata(h{i},'datatipLabel',dtLabelTemp);
        end
    end

    if ~isempty(alpha)
        if bitStatOnSig
            foo = bitSig;
        else
            foo = true(size(xTempForStats));
        end
        
        % fit line:
        if strcmpi(statTest,'corrPear')
            [rho(i),P(i)] = corr(xTempForStats(foo),yTempForStats(foo),'type','Pearson','rows','complete');
            % fits line to data excluding outliers
%             [lineParam(i,:)] = polyfit(xTemp,yTemp,1);
            % fits line to data including outliers
            [lineParam(i,:)] = polyfit(xTempForStats(foo),yTempForStats(foo),1);
        % rank-absed correlation
        elseif strcmpi(statTest,'corrSpear')
            [rho(i),P(i)] = corr(xTempForStats(foo),yTempForStats(foo),'type','Spearman','rows','complete');
            % fits line to data excluding outliers
%             [lineParam(i,:)] = polyfit(xTemp,yTemp,1);
            % fits line to data including outliers
%             [lineParam(i,:)] = polyfit(xTempForStats(foo),yTempForStats(foo),1);
        % compare medians
        elseif strcmpi(statTest,'signrank')
            % quick stat test:
            P(i) = signrank(xTempForStats(foo),yTempForStats(foo),'Method','exact'); % Wilcoxon signed-rank test
            %     P = ranksum(xTempForStats,yTempForStats); % Wilcoxon rank sum for *independent* samples (not this)
            %     P = signtest(xTempForStats-yTempForStats); % no rank. uses binomial distribution in "exact"
            %                          mode. Equivalent to a "coin flip" test.
            %     [~,P] = ttest2(xTempForStats,yTempForStats);
            % set line properties if significant:
        % compare medians
        elseif strcmpi(statTest,'signtest')
            % quick stat test:
            P(i) = signtest(xTempForStats(foo),yTempForStats(foo),'Method','exact'); % simple sign test
        % compare means
        elseif strcmpi(statTest,'ttest')
            [~,P(i)] = ttest2(xTempForStats(foo),yTempForStats(foo)); % t-test
        % total linear regression
        elseif strcmpi(statTest,'totregress')
            lineParam(i,:) = fliplr(totregress(yTempForStats(foo),xTempForStats(foo)));
            rho(i) = NaN;
            P(i) = NaN;
        end
        
        clear foo

    end

    % prepare probability text
    if ischar(colorTemp)
        colorStr = ['{',colorTemp,'}'];
    else
        colorStr = ['[rgb]{',num2str(colorTemp),'}'];
    end
    pTxt{i,1} = ['{\color',colorStr,'p=',sprintf('%0.2g',P(i))];
    % add correlation coeff
    if strcmpi(statTest,'corrPear')
        pTxt{i,1} = [pTxt{i,1},', r=',sprintf('%0.2g',rho(i))];
    elseif strcmpi(statTest,'corrSpear')
        pTxt{i,1} = [pTxt{i,1},', rho=',sprintf('%0.2g',rho(i))];
    elseif strcmpi(statTest,'totregress')
        pTxt{i,1} = [pTxt{i,1},', m=',sprintf('%0.2g',lineParam(i,1))];
    end
    pTxt{i,1} = [pTxt{i,1},'}'];
    
    % plot slope of correlation (Removed command to plot line for Spearman
    % correlation)
    if strcmpi(statTest,'corrPear') || strcmpi(statTest,'totregress') % ...
%             || strcmpi(statTest,'corrSpear')
        xLine = minmax(xTemp);
        yLine = polyval(lineParam(i,:),xLine);
        hL{i} = plot(xLine,yLine,'-');
        set(hL{i},'Color',colorTemp);
        if P(i) <= alpha
            set(hL{i},'LineWidth',2)
        end
    end
    
    % store mins
    xyTempLim(i,:) = [min(min(xTemp),min(yTemp)), max(max(xTemp),max(yTemp))];
    xTempLim(i,:) = minmax(xTemp');
    yTempLim(i,:) = minmax(yTemp');
end

% calculate limits
xyLim = [min(xyTempLim(:,1)), max(xyTempLim(:,2))];
xLim = [min(xTempLim(:,1)) max(xTempLim(:,2))];
yLim = [min(yTempLim(:,1)) max(yTempLim(:,2))];

% "round" to the nearest 1/10th:
xyLim = [floor(xyLim(1)*10) ceil(xyLim(2)*10)]/10;
xLim = [floor(xLim(1)*10) ceil(xLim(2)*10)]/10;
yLim = [floor(yLim(1)*10) ceil(yLim(2)*10)]/10;

% if not plotting correlation, plot unity line
if ~strcmpi(statTest,'corrPear') && ~strcmpi(statTest,'totregress')...
        && ~strcmpi(statTest,'corrSpear') 
    if bitEqualAxes 
        % plot unity line according to equal axes
        hL{i} = line(xyLim,xyLim);
    else
%         % plot unity line according to unequal axes
%         hL{i} = line(xLim,xLim);
    end
    
    % highlight line if significant for any column
    if any(P <= alpha)
        set(hL{i},'Color','c','LineWidth',2);
    end
end

if ~isempty(alpha)
    if bitEqualAxes
        % write probability text for equal axes
%         hT = text(xyLim(2),xyLim(2),pTxt);
        hT = text(1,1,pTxt,'Units','normalized');
    else
        % write probability text for unequal axes
%         hT = text(xLim(2),yLim(2),pTxt);
        hT = text(1,1,pTxt,'Units','normalized');
    end
    set(hT,'FontSize',pTextFontSize,'HorizontalAlignment','left','VerticalAlignment','bottom');
    % make italic for significance
    if any(P <= alpha)
        set(hT,'FontAngle','italic')
    end
end


if bitEqualAxes
    % make limits on x and y equal (if values were plotted)
    if ~any(isnan(xyLim))
        xlim(xyLim);
        ylim(xyLim);
    end
    axis square
else
    % make limits for unequal axes
    if ~any(isnan(xlim))
        xlim(xLim);
    end
    if ~any(isnan(ylim))
        ylim(yLim);
    end
end

% for each line, move it to back of axes so that points can be selected
for i = 1:length(hL)
    foo = get(gca,'Children');
    % eliminate line handles
    foo(ismember(foo,hL{i})) = [];
    % add back at end of list
    if size(hL{i},1) < size(hL{i},2)
        hL{i} = hL{i}';
    end
    foo = [foo; hL{i}];
    set(gca,'Children',foo);
end

% return axes to original overwrite mode
set(gca,'NextPlot',nextPlot);

% repackage output vars to be cells only if necessary
if length(h) == 1
    h = h{1};
end
if length(hL) == 1
    hL = hL{1};
elseif isempty(hL)
    hL = [];
end

% turn off data cursor mode
datacursormode off;

% -----------------------------
%% callback function for custom text on datatips
function datatipTxt = scatterDatatipCallback(obj,evt)

target = get(evt,'Target');
i = get(evt,'DataIndex');
pos = get(evt,'Position');

% retrieve attached data
% not all elements in plot necessarily have data attached
if isappdata(target,'obsind')
    obsind = getappdata(target,'obsind');
else
    obsind = [];
end
if isappdata(target,'datatipLabel')
    datatipLabel = getappdata(target,'datatipLabel');
else
    datatipLabel = {};
end

datatipTxt = {...
    ['x: ' num2str(pos(1))]...
    ['y: ' num2str(pos(2))]
    };
if ~isempty(obsind)
    datatipTxt{end+1} = ['Obs: ' num2str(obsind(i))];
end
if ~isempty(datatipLabel) 
    if iscell(datatipLabel)
        datatipTxt{end+1} = datatipLabel{i};
    else
        datatipTxt{end+1} = datatipLabel;
    end
end

% % get figure handle
% hA = get(hL,'Parent');
% hF = get(hA,'Parent');
% % open new figure if control button is held during click:
% if strcmp(get(hF,'SelectionType'),'normal')
%     figure('Position',[0 0 100 100]);
% end
% --------------------------------
% function h = openFigure(obj,evt)
% 
% 
% figure('Position',[100 100 100 100]);
