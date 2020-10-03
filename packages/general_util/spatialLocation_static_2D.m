function [figH,siteCorr,siteCorrP] = spatialLocation_static_2D(sig,loc,varargin)
% [figH] = spatialLocation_static_2D(sig,loc,varargin)
%
% Plot the spatial location of a signal of interest. The strength of signal
% SIG is represented by color and it's location is plotted in a
% grid-within-a-grid. At the higher level grid (G1), the signal is plotted as
% per the vertical and horizontal location of the signal LOC. Within each
% position of G1, is a subgrid (G2) that allows for multiple observations
% from the same location to be plotted. The subgrid consists of individual
% pixels, one for each observation of SIG. These pixels are ordered by the
% value of SIG.
%
% Multiple signals can be passed in SIG. In this case, a separate map is
% made for each signal. The pixels in each subgrid are plotted in order of
% the *first* signal. The pixels for subsequent signal maitain the order of
% the first signal, so that a given pixel position consistent refers to a
% single observation of the signal.
%
% This function assumes a 2D spatial map and a static signal that does not
% change over time.
%
% ACCEPTS
% SIG -- n x s matrix that includes the value of each signal s across
% multiple observations n.
%
% LOC -- n x 2 matrix of the location in vertical (column 1) and horizontal
% (column 2) dimensions of each observation n.  Expects integer units.
%
% RETURNS
% FIGH -- handle to figures created by function
%
% SITECORR -- correlation coefficient of mean site response (across sites)
% as correlated pairwise between signals. Returned as an S x S matrix where
% each element is the corr coeff for signal r x signal c, where r and c are
% the row and column numbers respectively.
%
% SITECORRP -- p value associated with corr coeff as provided in SITECORR.
%
% Daniel Kimmel, January 21, 2014

%% default values for optional params

% name of colormap with which to plot colors
% colormapName = 'summer';
colormapName = 'jet';

% 'asc' or 'desc' order in which to plot single observations within each
% subgrid g2. 'ascend' results in ascending order. 'descend' results in
% descending order 
g2order = 'descend'; 

% signal name. cell array of strings for each column in SIG
sigName = {};

% G1 row and column names
rowNameG1 = {};
colNameG1 = {};

% logical on whether to plot individual unit responses in grid format or
% plot mean response of site:
bitSiteMean = false;

% logical on whether to plot a scatter plot of the mean sRA component per
% site for signal 1 vs signal 2.
bitPlotSiteMeanScatter = false;

%%%%%%%%%%%%
% appearance
%%%%%%%%%%%%

bgColor = 'k';
textColor = 'w';

%% collect optional parameter values:
warnopts(assignopts(who, varargin));

%% checks

if size(loc,1) ~= size(sig,1)
    error('rows of LOC and SIG must be equal in number')
end

% remove observations that contain a NaN for location
foo = any(isnan(loc),2);
if any(foo)
    warning('Some location values (n=%d) were NaN and entire observation will be removed',sum(foo));
    
    loc(foo,:) = [];
    sig(foo,:) = [];
    
end
clear foo

% convert loc to integer
if ~isinteger(loc)
    try
        loc = int16(loc);
    catch
        error('Coordinates in LOC must be integers')
    end
end


%% define grid g1

% find unique locations of signals
[locG1,~,iG1] = unique(loc,'rows');

% dimensions of G1
nG1 = size(locG1,1);
[rowG1] = unique(locG1(:,1));
[colG1] = unique(locG1(:,2));
nRowG1 = length(rowG1);
nColG1 = length(colG1);

% define matrix of G1 mean and SD
meanG1 = NaN(size(sig,2),nRowG1,nColG1);
sdG1 = meanG1;

%% define G2 within each G1

nG2 = [];

for i = 1:nG1
    nG2(i) = sum(iG1==i);
    
    % dimensions of G2
    % for now, assume square subgrid
    nRowG2(i) = ceil(sqrt(nG2(i)));
    nColG2(i) = nRowG2(i);

end

% find position of max G2
[nG2Max, nG2MaxPos] = max(nG2);


%% establish sort order for each G2

for i = 1:nG1
    % sort order for each G2
    [~,iG2{i}] = sort(sig(iG1==i,1),1,g2order);
    
end

%% create pixel maps for each G2

% here we have to loop through signals as well as grids

for s = 1:size(sig,2)
    for i = 1:nG1
        % define NaN
        img{s,i} = NaN(nRowG2(i),nColG2(i));
        % sort signal
        foo = sig(iG1==i,s);
        foo = foo(iG2{i});
        % fill in signals
        img{s,i}(1:nG2(i)) = foo;
        % transpose (since we want to order across rows, then down columns)
        img{s,i} = img{s,i}';
        clear foo
        
        %%% Separately compute mean for each G1
        % find G1 position in row/col coord
        row = find(rowG1 == locG1(i,1));
        col = find(colG1 == locG1(i,2));
        meanG1(s,row,col) = nanmean(sig(iG1==i,s));
        sdG1(s,row,col) = nanstd(sig(iG1==i,s));
        clear row col
    end
end

%% find correlation of site mean between signals

siteCorr = NaN(size(sig,2));
siteCorrP = siteCorr;

for s1 = 1:size(sig,2)
    for s2 = s1+1:size(sig,2)
        [foo, goo] = corr([reshape(meanG1(s1,:,:),numel(meanG1(s1,:,:)),1),...
            reshape(meanG1(s2,:,:),numel(meanG1(s2,:,:)),1)],...
            'rows','complete',...
            'type','Pearson');
        siteCorr(s1,s2) = foo(1,2);
        siteCorrP(s1,s2) = goo(1,2);
        clear foo goo
            
    end
end
    
%% plotting

% loop through signals
for s = 1:size(sig,2)

    % make figure
    figH(s) = figure();
    whitebg(gcf,bgColor);
    set(gcf,'Color',bgColor);
    
    if bitSiteMean
        foo = 'MEAN ';
    else
        foo = '';
    end
    
    if ~isempty(sigName)
        figName = sprintf('Location of %s signal %s',sigName{s},foo);
    else
        figName = sprintf('Location of %s signal %d',s,foo);
    end
    clear foo
        
    set(gcf,'Name',figName);
    
    % define colors of each pixel
%     eval(sprintf('colorset = %s(%d);',colormapName,max(nG2)));
    
    %%% FOR SITE MEAN
    if bitSiteMean
        plotMean();
        
    
    %%% FOR INDIVIDUAL UNITS
    else
        plotGrid();
    end
    
    % add annotation
    figTitleHand = annotation('textbox',[0 0 1 .05]);
    figAnnotStr = {[figName];
        [datestr(now),', ',mfilename,'.m']
        };
    set(figTitleHand,'String',figAnnotStr,'Interpreter','none','FontSize',10);
    set(figTitleHand,'VerticalAlignment','top','HorizontalAlignment','right');
    set(figTitleHand,'LineStyle','none');
    set(figTitleHand,'Color','r');
    
end

%% scatter of site means

if bitPlotSiteMeanScatter
    
    % number of signals
    S = size(meanG1,1);
    % number of axes
    A = nchoosek(S,2);
    nRow = 1;
    nCol = A;
   
    % make figure
    figName = 'Scatter of mean signal per site';
    figH(s) = figure();
    whitebg(gcf,bgColor);
    set(gcf,'Color',bgColor,'Name',figName);
    
    axN = 0;
    for s1 = 1:S-1
        for s2 = s1+1:S
            axN = axN + 1;
            subplot(nRow,nCol,axN)
            
            [h,~,~,~,hT] = plotScatterN(reshape(meanG1(s1,:,:),numel(meanG1(s1,:,:)),1),...
                reshape(meanG1(s2,:,:),numel(meanG1(s2,:,:)),1),...
                'statTest','corrPear','pTextFontSize',12,'alpha',0.05,...
                'bitEqualAxes',true);
            set(h,'Marker','o','MarkerSize',7);
            xlabel(sigName{s1},'FontSize',14);
            ylabel(sigName{s2},'FontSize',14);
            set(gca,'FontSize',14);
            set(hT,'Position',[0 1],'Units','normalized'); % probability text location

        end
    end

    % add annotation
    figTitleHand = annotation('textbox',[0 0 1 .05]);
    figAnnotStr = {[figName];
        [datestr(now),', ',mfilename,'.m']
        };
    set(figTitleHand,'String',figAnnotStr,'Interpreter','none','FontSize',10);
    set(figTitleHand,'VerticalAlignment','top','HorizontalAlignment','right');
    set(figTitleHand,'LineStyle','none');
    set(figTitleHand,'Color','r');
end

%% SUBFUNCTION: plot grid of individual unit responses by site
    function plotGrid()
        % loop through row and cols of G1
        axN = 0;
        for r = 1:nRowG1
            for c = 1:nColG1
                axN = axN + 1;
                
                i = find(ismember(locG1,[rowG1(r),colG1(c)],'rows'));
                
                subplot(nRowG1,nColG1,axN);
                colormap(colormapName);
                
                % only plot if G1 at that position
                if ~isempty(i)
                    
%                     h = imagesc(img{s,i},minmax(sig(:,s)')); % scale specific for signal
                    h = imagesc(img{s,i},minmax(sig(:)')); % scale across all signals
                    
                    % scale axes
                    axPos = get(gca,'Position');
                    scalarRow = nRowG2(i) / nRowG2(nG2MaxPos);
                    scalarCol = nColG2(i) / nColG2(nG2MaxPos);
                    set(gca,'Position',...
                        [axPos(1)...
                        axPos(2)...
                        axPos(3)*scalarCol...
                        axPos(4)*scalarRow]);
                    
                    % appearance
                    set(gca,'XTick',[],'YTick',[]);
                    set(gca,'box','off');
                    set(gca,'XColor',bgColor,'YColor',bgColor);
                    
                    % hide NaN's
                    set(h,'AlphaData',~isnan(img{s,i}));
                    set(gca,'Color',bgColor)
                else
                    % set axis visibility
%                     set(gca,'Visible','off')
                    set(gca,'box','off','XColor',bgColor,'YColor',bgColor)
                end
                
                % label rows and columns
                if r==1
                    if ~isempty(colNameG1)
                        title(colNameG1{c},'FontSize',16);
                    else
                        title(sprintf('%d',colG1(c)),'FontSize',16);
                    end
                end
                if c == 1
                    if ~isempty(rowNameG1)
                        ylabel(rowNameG1{r},'FontSize',16,'Color',textColor);
                    else
                        ylabel(sprintf('Row %d',rowG1(r)),'FontSize',16,'Color',textColor);
                    end
                end

            end
        end
        
        % add colorbar as additional axes
        axH2 = axes('Position',[0 0.88 0.0 0.1]);
        % plot fake image
%         h = imagesc([1],minmax(sig(:,s)')); % scale specific for signal
        h = imagesc([1],minmax(sig(:)')); % scale across all signals
        set(h,'Visible','off')
        % add colorbar
        hc = colorbar(axH2);
        % hide object
        set(axH2,'Visible','off')

    end

%% SUBFUNCTION: plot mean reponse of site

    function plotMean()
        %%% MEAN
        subplot(1,2,1);
%         colormap(colormapName);        
        h = imagesc(colG1,rowG1,squeeze(meanG1(s,:,:)));
        title('Mean','FontSize',16)
        xlabel('Grid Column','FontSize',14);
        ylabel('Grid Row','FontSize',14);
        set(gca,'XTick',colG1,'YTick',rowG1,'box','off');
        % replace numeric grid and col labels with text labels if provided
        if ~isempty(rowNameG1)
            set(gca,'YTickLabel',rowNameG1)
        end
        if ~isempty(colNameG1)
            set(gca,'XTickLabel',colNameG1)
        end
        % hide NaN's
        set(h,'AlphaData',~isnan(squeeze(meanG1(s,:,:))));
        set(gca,'Color',bgColor)
        axis square
        % add colorbar
        colorbar;
        
        %%% STD
        subplot(1,2,2);
%         colormap(colormapName);        
        h = imagesc(colG1,rowG1,squeeze(sdG1(s,:,:)));
        title('Standard Deviation','FontSize',16)
        set(gca,'XTick',colG1,'YTick',rowG1,'box','off');
        % replace numeric grid and col labels with text labels if provided
        if ~isempty(rowNameG1)
            set(gca,'YTickLabel',rowNameG1)
        end
        if ~isempty(colNameG1)
            set(gca,'XTickLabel',colNameG1)
        end        
        % hide NaN's
        set(h,'AlphaData',~isnan(squeeze(sdG1(s,:,:))));
        set(gca,'Color',bgColor)
        axis square
        % add colorbar
        colorbar;
        
    end
end

