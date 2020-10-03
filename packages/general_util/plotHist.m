function [h,hL,bitObsPlotted,P] = plotHist(data,varargin)
% [h,hL,bitObsPlotted,P] = plotHist(DATA, ... [name,value pairs]) 
% 
% Plots "point" histogram of DATA in current axes.  This is distinct from a
% bar histrogram in that the number of observations per bin is represented
% as a stack of individual data points.  This enables the user to click an
% individual data point and trigger a custom callback routine specific to
% that one data point. 
%
% DATA is a M x N matrix with M observations across N conditions. Multiple
% conditions are plotted as separate histograms in distinct colors.
%
% RETURNS:
% H = vector of handles to indvidual data points. If N > 1, H is a 1xN cell
% array with each cell containing the vector relvant to column n.
%
% HL = vector of handles to the vertical lines plotted at the null
% hypothesis (see below) for each column in DATA.
%
% bitObsPlotted = M x N logical matrix specifying whether data point in
% DATA was plotted and consequently has a handle in H.
%
% P = 1 x N vector of p values of the probability the distribution in a
% column of DATA differs from its null hypothesis (see below) by chance.
%
% ACCEPTS the following as optional name,value pairs:
%
% bitSig = M x N logical matrix about whether individual observation was
% signficant and should be plotted by a filled circle.  All signficant
% observations are plotted before (below) non-significant ones.
%
% datatipLabel = M x N cell array of custom datatip labels for individual
% data points.
%
% onClickFunction = function handle of callback function to be called when
% user clicks a data point. 
%
% onClickData = M x N cell array of custom data to be passed to
% onClickFunction when user clicks data point.  
%
% color = Either 1 x N cell array of color names or N x 3 matrix of RGB
% values used to plot each column in DATA.  If only one color is provided,
% it will be used for all columns.  If no color is specified, colors will
% be drawn from the axes ColorOrder property.
% 
% bins = 1 x B vector specifying edges of B bins into which to place DATA.
% If not provided, bins are computed based on nBin
%
% nBin = Number of bins.  Overwritten by BINS.  Defaults to 10.
%
% bitDoubleDir = logical about whether second condition should be plotted
% upside down. Only applies when N = 2. (default = true).
%
% outlierZThresh = threshold for outliers in z-scores, beyond which data 
% points will be excluded from the computation of bin edges. Can be used
% concurrently with outlierNMean -- one does not affect the other's
% behavior.(default = Inf).
%
% outlierNMean = threshold for outliers in multiples of the mean, beyond
% which data points will be excluded from the computation of bin edges. Can
% be used concurrently with outlierZThresh -- one does not affect the
% other's behavior. (default = Inf). 
% 
% outlierPrctile = 2 element vector specify lower and upper percentile (0
% to 100) at which to truncate the data. Data points below or above,
% respectively, these percentiles will be counted as outliers.  Pass NaN as
% either or both element to not invoke a lower or upper bound. This
% variable is the most stringent of outlier thresholds because it
% guarantees that a certain percentage of values will be eliminated.
% Default [NaN NaN].
% 
% nullH = 1 x N vector specifying the null hypothesis for each column
% in DATA, that is, the number to which the distribution should be
% compared for statistical tests.  Specifying nullH will place a vertical
% line at nullH in the color of the condition.
%
% testType = string with the type of test to use for comparing distrubtion
% to nullH.  Options are 'signrank', 'ttest', 'signtest'.  If specified,
% probability will be reported in return var P, printed on figure, and, if
% less than alpha, result in a boldface line for null hypothesis.
%
% alpha = level of probability needed to reject null hypothesis. Defaults
% to 0.05.
%
% markerSize = size of marker used to plot individual data points (defaults
% to 5).
%
% bitPlotBar = logical on whether to plot histogram as traditional bars
% (instead of markers).  onClickFunction will not work in this mode.  Only
% tested for single condition plotting as of now (12 Nov 2010).
%
% bitSymmetricX = logical on whether to set xaxis as symmetrical on either
% side of the null hypothesis.  Only supported when plotting a single
% condition.
%
% bitPlotMean = logical on whether to draw arrow and write text of mean of
% each column in DATA. Overrides bitPlotMedian.
%
% bitPlotMedian = same as bitPlotMean, but plots median instead.
% Over-ridden by bitPlotMean.
%
% bit4Paper = logical on whether to use colors suitable for white paper.
% default = false;
%
% Daniel Kimmel, February 6, 2009

%% default values:
alpha = [];
bitSig = false(size(data));
datatipLabel = {};
onClickFunction = @highlightPoint;
onClickData = {};
color = [];

bins = []; % specify bin edges
nBin = 10; % number of bins.  overwritten by BINS.
bitDoubleDir = true; % logical about whether second condition (column in X
                     % and Y) should be plotted upside down (default =
                     % true).
outlierZThresh = Inf; % threshold for outliers in z-scores, beyond which data 
                    % points will be excluded from the computation of bin
                    % edges.
outlierNMean = Inf; % see above
outlierPrctile = [NaN NaN]; % see above
                    
nullH = []; % vector as long as columns in X and Y specifying the null 
            % hypothesis for each column, that is, the number to which the
            % distribution should be compared.
testType = []; % string with the type of test to use for comparing 
               % distrubtion to nullH.  Options are 'signrank', 'ttest',
               % 'signtest'.
alpha = 0.05; % level of probability needed to reject null hypothesis.

markerSize = 4; 

pFontSize = 9; % font size for probability text.

bitPlotBar = false;

bitSymmetricX = false;

bitPlotMean = false; % draws arrow and writes text of mean of each column in DATA
bitPlotMedian = false;

bit4Paper = false; % logical on whether to use colors suitable for white paper.
%% collect optional name,value pairs
warnopts(assignopts(who, varargin));

%% checks

% data must be in column form
if any(size(data)) == 1 && size(data,1) == 1 && size(data,2) > 1
    error('DATA input must be in column form');
end

% bitSig must be in column form
if ~isempty(bitSig) && ~isequal(size(data),size(bitSig))
    error('bitSig input must be in column form and match size of data');
end

% currently, we only tested bitPlotBar for a single condition
if bitPlotBar && size(data,2) > 1
    error('bitPlotBar only tested for single condition currently')
end

% bitSig might be passed as empty
if isempty(bitSig)
    bitSig = false(size(data));
end

%% prepare and plot

% determine colors
if bit4Paper
    textColor = 'k';
    highlightColor = 'k';
    faceColorDefault = [0.5 0.5 0.5];
else
    textColor = 'w';
    highlightColor = 'm';
    faceColorDefault = 'yellow';
end    

% save axes overwrite property
nextPlot = get(gca,'NextPlot');
hold on;

% handle for attaching custom datatips:
hdt = datacursormode;
set(hdt,'UpdateFcn',{@scatterDatatipCallback});

% instatiate:
bitObsPlotted = false(size(data));
obsind = cell(1,size(data,2));
h = cell(1,size(data,2));
hL = NaN(1,size(data,2));
P = NaN(1,size(data,2));

% if data is all NaNs, return
if all(all(isnan(data)));
    return
end

pTxt = cell(size(data,2),1);
xyTempLim = NaN(size(data,2),2);
mean_color = cell(size(data,2),2);

% prepare color
if iscell(color) && ~isempty(color)
    % repeat color for all columns
    if length(color) < size(data,2)
        foo = color;
        clear color;
        color(1:size(data,2)) = foo(1);
        clear foo
    end
elseif ischar(color)
    foo = color;
    clear color
    % repeat color for all columns
    color(1:size(data,2)) = {foo};
    clear foo
elseif isempty(color)
    % use default color if only one
    if size(data,2) == 1
        color = faceColorDefault;
    else
        goo = get(gca,'ColorOrder');
        color = goo(1:size(data,2),:);
        clear goo
    end
end
    
% find range in data:
% optionally exclude outliers
if outlierZThresh < Inf | outlierNMean < Inf | any(~isnan(outlierPrctile))
    % linearize all data into single vector
    dataLin = reshape(data,prod(size(data)),1);
    % compute mean 
    dataMean = nanmean(dataLin);
    % find percentile cutoffs for outliers
    if isfinite(outlierPrctile(1))
        pCutoff(1) = prctile(dataLin,outlierPrctile(1));
    else
        pCutoff(1) = -Inf;
    end
    if isfinite(outlierPrctile(2))
        pCutoff(2) = prctile(dataLin,outlierPrctile(2));
    else
        pCutoff(2) = Inf;
    end
    % transform in to z-scores:
    dataLinT = (dataLin - dataMean) ./ nanstd(dataLin);
    % transform into multiples of the mean
    dataLinT2 = dataLin ./ dataMean;
    % find outliers
    bitRemove = false(size(dataLin));
    if any(dataLinT > outlierZThresh) || any(dataLinT2 > outlierNMean) ||...
            any(dataLinT > pCutoff(2))
        % flag that a large outlier was removed
        bitOutlierMax = true;
        bitRemove = bitRemove | dataLinT > outlierZThresh | dataLinT2 > outlierNMean | ...
            dataLin > pCutoff(2);
    else 
        bitOutlierMax = false;
    end
    if any(dataLinT < -outlierZThresh) || any(dataLinT2 < -outlierNMean) ||...
            any(dataLinT < pCutoff(1))
        % flag that a small outlier was removed
        bitOutlierMin = true;
        bitRemove = bitRemove | dataLinT < -outlierZThresh | dataLinT2 < -outlierNMean | ...
            dataLin < pCutoff(1);
    else
        bitOutlierMin = false;
    end
    % remove outliers
    dataLin(bitRemove) = [];
    
    % find min and max
    dataMin = min(dataLin);
    dataMax = max(dataLin);
    clear dataLin dataLinT bitRemove
else
    dataMin = min(min(data));
    dataMax = max(max(data));
    bitOutlierMax = false;
    bitOutlierMin = false;
end

% prepare null hypothesis.
% if only one is given, assume it is the same for all conditions.
if length(nullH) == 1 && size(data,2) > 1
    nullH = ones(1,size(data,2)) * nullH;
end

% define bins
if isempty(bins)
    bins = linspace(dataMin,dataMax,nBin);
    % if datamin and datamax are equal, there will only be one value for
    % bins, repeated nBin times.  the Bar() function doesn't like this.
    if length(unique(bins))==1
        bins = bins(1);
    end
    
    % if null hypothesis is given and it is the same for all columns in X
    % and Y, AND if the nullH is within the data range AND if nullH not
    % already present in bins, then make sure nullH is one of the edges in
    % bins. 
    if ~isempty(nullH) && all(nullH == nullH(1)) ...
            && nullH(1) > dataMin && nullH(1) < dataMax ...
            && ~ismember(nullH,bins)
        % take rough linspace and find entry just below nullH
        [pos] = find(bins < nullH(1),1,'last');
        % add the difference between this entry and nullH to bins, now bins
        % includes nullH
        bins = bins + nullH(1) - bins(pos);
        % now put back the first bin to include dataMin. Note that in one
        % instance this step caused a rounding error in which the original
        % bins(1) contained dataMin, but after subtracting binwidth from
        % the new bins(1), the new new bins(1) was e-17 > than dataMin.
        % This particular issue happened when bins already contained nullH,
        % and so was handled by only carrying out this section when bins
        % does not already include nullH.
        bins = [bins(1) - (bins(2)-bins(1)), bins];
    end
    
    % add Inf bins for outliers, but only if outlier would be outside of
    % existing bins
    if bitOutlierMax && max(max(data)) > bins(end)
        % make one additional bin at normal spacing, so that the previous
        % bin has the proper width.  Then add a final bin at Inf so that
        % the new additional bin captures all outliers.
        bins = [bins,bins(end)+(bins(end)-bins(end-1)),Inf];
    else
        bitOutlierMax = false;
        dataMax = max(max(data));
    end
    if bitOutlierMin && min(min(data)) < bins(1)
        bins = [-Inf, bins];
    else
        bitOutlierMin = false;
        dataMin = min(min(data));
    end
end

% compute bin spacing
if isinf(bins(1)) && length(bins)>2
    binSpace = bins(3)-bins(2); % don't use 2 - 1, since 1 could be inf.
elseif length(bins)>1 && ~isinf(bins(1))
    binSpace = bins(2)-bins(1); 
else
    binSpace = 0;
end
% if there is only one value for all of bins, then binSpace would be zero.
% Make it 1 just to be nice for plotting:
if binSpace == 0
    binSpace = 1;
end

% decide spacing away from bin edge to plot "bar"
% This makes points plot in middle of bin (which is unlike traditional
% histc() with bar() plotting, where bins are plotted at the lower edge of
% the bin.
xPlot = bins + (binSpace)/2;
% This makes points plot at the lower edge of the bin (as in histc() with
% bar()):
% xPlot = bins;

% if bin is to hold outlier, create special point and label
if bitOutlierMax 
    % if last bin is inf, it will have zero freq, and it's the second to
    % last bin that should be moved over.
    if isinf(xPlot(end))
        xPlot(end-1) = xPlot(end-2) + 2 * binSpace;
    else
        error('With a high-end outlier, the last bin should be Inf')
    end
end
if bitOutlierMin 
    xPlot(1) = xPlot(2) - 2 * binSpace;
end

% decide additional offset away from bar "center" to slide
% individual conditions.
if size(data,2) == 1
    xPlotOffset = 0;
    offsetInt = binSpace; % here used only for bar width computation below.
% for only two conditions, we'll plot second condition upside-down, so
% offset is not necessary
elseif size(data,2) == 2 && bitDoubleDir
    xPlotOffset = zeros(1,size(data,2));
    offsetInt = binSpace; % here used only for bar width computation below.
% for even number of conditions
elseif mod(size(data,2),2) == 0
    % offset interval
    offsetInt = (binSpace) ./ (size(data,2) + 3);
    % offset
    xPlotOffset = [-size(data,2)/2:size(data,2)/2] .* offsetInt;
    % remove zero point
    xPlotOffset(xPlotOffset==0) = [];
% for odd number of conditions
else
    % offset interval
    offsetInt = (binSpace) ./ (size(data,2) + 1);
    % offset
    xPlotOffset = [ceil(-size(data,2)/2):floor(size(data,2)/2)] .* offsetInt;
end


% DATA can be a matrix with multiple columns to be plotted in different
% colors.
for i = 1:size(data,2)
        
    % prepare data, eliminate nans
    foo = isnan(data(:,i));
    dataTemp = data(~foo,i);
    bitSigTemp = logical(bitSig(~foo,i));
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
    obsind{i} = 1:length(dataTemp);
    bitObsPlotted(:,i) = ~foo;
    clear foo;
    
    % index of data order used for later sorting handles into correct order
    dataOrder = [1:length(dataTemp)]';
    dataOrderAfter = []; % will compiles changed order here

    % if no data, continue
    if isempty(dataTemp)
        continue
    end
    
    % prepare color string
    if iscell(color)
        colorTemp = color{i};
    else
        % RGB vector
        colorTemp = color(i,:);
    end
    
    % decide direction that histogram "bars" will be plotted.  If only two
    % conditions (columns in X and Y), plot one up and one down.
    if size(data,2) == 2 && i == 2 && bitDoubleDir
        barDir = -1;
    else
        barDir = 1;
    end
    
    % compute histogram:
    [freq,posData] = histc(dataTemp,bins);
    
    %%%%%%%%%%%%%%%%%%%%%
    % IF NOT PLOTTING SINGLE POINTS, JUST PLOT THE VARS HERE
    if bitPlotBar
        
        %         h{1} = bar(bins',freq,'histc');
        % plot instead with xPlot (instead of bins) because it handles the
        % position of outliers.  Note that xPlot can still contain Inf, and
        % so must be checked.
        % Use regular bar command (without 'histc') so that bars don't span
        % the gap when an outlier is present.  Set width to 1 so that bars
        % touch each other.
        bitInf = isinf(xPlot);
        h{1} = bar(xPlot(~bitInf)',freq(~bitInf),1);
        set(h{1},'FaceColor',colorTemp,'EdgeColor',get(gca,'Color'));
                
        % find hist of significant points
        [freqSig] = histc(dataTemp(bitSigTemp),bins);
%         h{2} = bar(bins',freqSig,'histc');
        if ~isempty(freqSig)
            % see above for change to xPlot and 'hist'
            h{2} = bar(xPlot(~bitInf)',freqSig(~bitInf),1);
            set(h{2},'FaceColor',highlightColor,'EdgeColor',get(gca,'Color'));
        else
            h{2} = [];
        end
        % outliers. replot bars with just the outlier bars, and color their
        % edges red.
        if bitOutlierMax
            % skip the last point, since it's Inf
            h{4} = bar(xPlot(1:end-1)',[NaN(length(xPlot)-2,1);freq(end-1)],1);
            set(h{4},'FaceColor','none','EdgeColor','r','LineWidth',1.5);
%             set(get(h{end},'BaseLine'),'LineStyle','none');
        end
        hold on;
        if bitOutlierMin
            bitInf = isinf(xPlot);
            foo = [freq(1); NaN(length(xPlot)-1,1)];
            h{3} = bar(xPlot(~bitInf)',foo(~bitInf),1);
            set(h{3},'FaceColor','none','EdgeColor','r','LineWidth',1.5);
            clear foo bitInf
        end
        
    else
    
        % creating counting var for how many points have been added to each bin:
        freqCount = zeros(size(bins));

        % loop through each bin
        for binN = 1:length(bins)
            % find the data that belong to this bin.
            % combine data and bitSig flag into single matrix:
            dataSigBin = [dataTemp(posData==binN) bitSigTemp(posData==binN)];

            % only sort data by significance if some are significant
            if ~isempty(dataSigBin)
                % sort data on bitSig, then data value
                [dataSigBin sortOrder] = sortrows(dataSigBin,[-2 1]);
            else
                sortOrder = [];
            end

            % apply sort to observation number
            obsindBin = obsind{i}(posData==binN);
            obsindBin = obsindBin(sortOrder);

            % apply sort to datatip labels and onClickData
            if ~isempty(dtLabelTemp)
                dtLabelBin = dtLabelTemp(posData==binN);
                dtLabelBin = dtLabelBin(sortOrder);
            else 
                dtLabelBin = {};
            end
            if ~isempty(onClickDataTemp)
                onClickDataBin = onClickDataTemp(posData==binN);
                onClickDataBin = onClickDataBin(sortOrder);
            else
                onClickDataBin = {};
            end

            % apply sort to index of data points:
            dataOrderBin = dataOrder(posData==binN);
            dataOrderBin = dataOrderBin(sortOrder);
            % store changed order:
            dataOrderAfter = [dataOrderAfter; dataOrderBin];

            % make sure number of data points assigned to this bin matches
            % actual freq 
            if freq(binN) ~= size(dataSigBin,1)
                error('Manual process of adding points to graph did not result in same frequencies as returned by HISTC()')
            end

            % plot zero point if freq == 0 for bin
            if freq(binN) == 0
                foo = plot(xPlot(binN) + xPlotOffset(i), 0, '.');
                % color
                set(foo,'Color',colorTemp,'MarkerSize',markerSize*1.1);
                clear foo
            end

            % loop through data for this bin, plotting each point separately
            for j = 1:size(dataSigBin,1)
                % for last point, draw stem first, but do not include handle
    %             if j == size(dataSigBin,1)
    %                 foo = stem(xPlot(binN) + xPlotOffset(i), barDir * j,
    %                 'o');
    %                 foo = line(ones(1,2) * (xPlot(binN) + xPlotOffset(i)),barDir * [0 j]);
    %                 set(foo,'Color',colorTemp);
    %                 foo = bar(xPlot(binN) + xPlotOffset(i), barDir * j,offsetInt);
    %                 set(foo,'FaceColor','none');
    %                 clear foo
    %             end


                % plot point with handle. 
                h{i}(end+1) = plot(xPlot(binN) + xPlotOffset(i), barDir * j, 'o');
                clear foo
                % color
                set(h{i}(end),'Color',colorTemp,'MarkerSize',markerSize);
                % fill in if sig:
                if dataSigBin(j,2)
                    set(h{i}(end),'MarkerFaceColor',colorTemp);
                end
                % make square if outlier. High-end outliers live in the
                % second to last bin, since the last bin is Inf.
                if (binN == 1 && bitOutlierMin) || ...
                        (binN == length(bins)-1 && bitOutlierMax)
                    set(h{i}(end),'Marker','s');
                end

                % Attach custom datatip to specific point
                setappdata(h{i}(end),'obsind',obsindBin(j));
                setappdata(h{i}(end),'obsValue',dataSigBin(j,1));
                if ~isempty(dtLabelBin)
                    setappdata(h{i}(end),'datatipLabel',dtLabelBin{j});
                end

                % if an onClick callback function is specified, attach it to
                % the data point:
                if ~isempty(onClickFunction)
                    % if data is included to be passed to the callback function...
                    if ~isempty(onClickDataBin)
                        set(h{i}(end),'ButtonDownFcn',{onClickFunction,onClickDataBin{j}});
                    else
                        set(h{i}(end),'ButtonDownFcn',onClickFunction);
                    end
                end
            end
        end

        % data were not plotted in order that they were passed and the
        % handles list is therefore not in that order. Return the list to
        % the proper order:
        [~,foo] = sort(dataOrderAfter);
        h{i} = h{i}(foo);
        clear foo
    end
    
    % Compute probability
    if ~isempty(nullH) && ~isempty(testType)
        
        % quick stat test:
        switch testType
            case 'signrank'
                P(i) = signrank(dataTemp,nullH(i),'Method','exact'); % Wilcoxon signed-rank test
                %     P = ranksum(xTemp,yTemp); % Wilcoxon rank sum for *independent* samples (not this)
                %     P = signtest(xTemp-yTemp); % no rank. uses binomial distribution in "exact"
                %                          mode. Equivalent to a "coin flip" test.
                %     [~,P] = ttest2(xTemp,yTemp);
            case 'ttest'
                [~,P(i)] = ttest(dataTemp,nullH(i));
            case 'signtest'
                P(i) = signtest(dataTemp,nullH(i),'Method','exact');
        end
        
        
        % prepare probability text
        if P(i) < alpha
            emphasisStr = sprintf('\\fontsize{%f}\\it',pFontSize);
        else
            emphasisStr = '';
        end
        if ischar(colorTemp)
            colorStr = ['{',colorTemp,'}'];
        else
            colorStr = ['[rgb]{',num2str(colorTemp),'}'];
        end
        pTxt{i,1} = ['{\color',colorStr,emphasisStr,'p=',num2str(P(i),2),'}'];
        
    end
    
    % compute mean, store color
    if bitPlotMean 
        mean_color{i,1} = mean(dataTemp);
        mean_color{i,2} = colorTemp;
    elseif bitPlotMedian
        mean_color{i,1} = median(dataTemp);
        mean_color{i,2} = colorTemp;
    end
end

% set limits
if bitSymmetricX && size(data,2) == 1
    % plot an equal amount of space on either side of the nullH for the
    % x-axis.  only works when plotting single condition.
    
    % find max distance from nullH
    if bitOutlierMax % last bin in Inf
        maxD = max(abs([xPlot(1) - binSpace/2, xPlot(end-1) + binSpace/2] - nullH));
    else
        maxD = max(abs([xPlot(1) - binSpace/2, xPlot(end) + binSpace/2] - nullH));
    end
    % apply
    xlim([nullH-maxD, nullH+maxD]);
    
elseif bitSymmetricX
    error('bitSymmetricX only works when plotting a single condition')
else
    if bitOutlierMax % last bin in Inf
        xlim([xPlot(1) - binSpace/2, xPlot(end-1) + binSpace/2]);
    else
        xlim([xPlot(1) - binSpace/2, xPlot(end) + binSpace/2]);
    end
end
% get limits
xLim = get(gca,'XLim');
yLim = get(gca,'YLim');

% make white line along xaxis and fix tic mark direction
if bitPlotBar
    xLim = get(gca,'Xlim');
    line(xLim,[0 0],'Color',get(gca,'XColor'));
    set(gca,'TickDir','out');
end

% plot nullH line
if ~isempty(nullH)
    for i = 1:size(data,2)
        hL(i) = line([nullH(i) nullH(i)],yLim);
        % prepare color string
        if iscell(color)
            colorTemp = color{i};
        else
            % RGB vector
            colorTemp = color(i,:);
        end
        set(hL(i),'Color',colorTemp,'LineStyle',':');

        % highlight line if significant for any column
        if P(i) < alpha
            set(hL(i),'LineWidth',2,'LineStyle','--');
        end
        
        % special when plotting bars
        if bitPlotBar
            set(hL(i),'Color','r','LineStyle','-','LineWidth',2);
            if ~bit4Paper
                set(hL(i),'Color','r');
            else
                set(hL(i),'Color',textColor);
            end
        end
    end
end

% write probability text
hT = text(xLim(2),yLim(2),pTxt);
set(hT,'FontSize',pFontSize,'HorizontalAlignment','right','VerticalAlignment','top');

% plot means
if bitPlotMean || bitPlotMedian
    for i = 1:size(mean_color,1)
        % plot mean of distribution
        arrowLength = range(yLim) * 1/14;
    %     hQ = quiver(nanmean(data),yLim(2)+arrowLength,0,-arrowLength);
    %     set(hQ,'MaxHeadSize',100,'LineWidth',2)
        foo.dist = 0;
        foo.dir = 270;
        foo.length = arrowLength;
        foo.width = range(get(gca,'XLim')) * .03;
        foo.lengthratio = 0.6;
        foo.widthratio = 0.5;
        arrow = makearrow(foo);
        plot(arrow(1,:)+mean_color{i,1},arrow(2,:)+yLim(2)+arrowLength,'Color',mean_color{i,2});

        % write mean:
        hT = text(mean_color{i,1},yLim(2)+arrowLength,num2str(mean_color{i,1},2));
        set(hT,'HorizontalAlignment','center','VerticalAlignment','bottom',...
            'FontSize',pFontSize,'Color',mean_color{i,2});
        clear foo hT
    end
end

% return axes to original overwrite mode
set(gca,'NextPlot',nextPlot);

% repackage output vars to be cells only if necessary
if length(h) == 1
    h = h{1};
end

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
if isappdata(target,'obsValue')
    obsValue = getappdata(target,'obsValue');
else
    obsValue = [];
end
if isappdata(target,'datatipLabel')
    datatipLabel = getappdata(target,'datatipLabel');
else
    datatipLabel = {};
end

datatipTxt = {...
    ['binC: ' num2str(pos(1))]...
    ['count: ' num2str(pos(2))]
    };
if ~isempty(obsind)
    datatipTxt{end+1} = ['Obs N: ' num2str(obsind(i))];
end
if ~isempty(obsValue)
    datatipTxt{end+1} = ['Obs Value: ' num2str(obsValue(i))];
end
if ~isempty(datatipLabel) 
    if iscell(datatipLabel)
        datatipTxt{end+1} = datatipLabel{i};
    else
        datatipTxt{end+1} = datatipLabel;
    end
end
