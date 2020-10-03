%%% This code plots the figures in Kimmel et al., Nat Comm, 2020, based on
%%% the processed results obtained either by running oTDR against the raw
%%% data or by loading the preprocessed results, which are available along
%%% side the raw data in the separate data repository
%%% (https://github.com/danielkimmel/Kimmel_NatComm_2020).

%% add oTDR to path
%%%
%%% For this block to run, file must be executed from command line.
%%% Once path is loaded by calling startup.m in root directory of package,
%%% this block can be skipped
%%%

% name of current file (for later testing purposes)
curfilename = 'test_oTDR.m';

% path to oTDR package
otdr_package_path = '';
if isempty(otdr_package_path)
    otdr_package_path = input('Enter absolute path to oTDR package: ','s');
end
if isempty(otdr_package_path)
    error('Must provide path to oTDR package');
end

% get current directory 
curdir = pwd;

% get directory of current file (only works if file is being executed from
% command line, not if executing by line or code block)
filedir = fileparts(mfilename('fullpath'));

try
    % navigate to top-level directory of oTDR package
    cd(otdr_package_path);
    % run startup file there
    startup;
catch
    error('Run startup.m in root directory of oTDR package. Alternatively, execute %s from command line.',curfilename);
end
% return to original directory
cd(curdir);

% confirm that path was appropriately updated
if ~exist(curfilename,'file')
    error('MATLAB path was not successfully updated. Run startup.m in root directory of oTDR package. Alternatively, execute %s from command line.',curfilename);
end

%% define data, colors

% set animal name ('Norris', 'Kirby')
animalName = 'Norris';

% logical on whether to limit analysis to single units
bitSingleUnit = false;
% logical on whether to limit analysis to units that have the singleton
% condition present.
bitSingleton = false;
% logical on whether previous offer epoch extends until the current offer
% (TRUE) or until the end of the current trial (FALSE). This effects the
% medium-D subspace used for noise-reduction (i.e., the subspace is defined
% by the top D PCs, which in turn depend on the observations fed to PCA).
bitPP2CurrentOffer = true;

% 'present trial', 'single unit', 'prev trial', 'present trial in prev
% trial' -- This only works for recalling DK's RDR results. See below to
% access animals's different results.
analName = 'present trial'; 

condLabel = {'Benefit','Choice','Expected Reward'};

% set whether to denoise data (by projecting it into med-D space before
% computing RAs).
bitDenoise = false;

% set whether to use all units ('allUnits') or single units ('singleUnits')
unitType = 'allUnits';

% set whether to remove time varying signal (common across conditions)
% before doing PCA and TDR 
bitRemoveTS = true;

% define alpha
alpha = 0.05;

%%% REMOVE THE FOLLOWING because now all probabilities will be based on
%%% gamma fits to the chance distributions
% % set whether to plot deviation from chance in terms of sigma or
% % probability
% bitSigma = true;
%
% % define chance level in terms of sigma. 
% chanceLevel = norminv(1 - alpha); % for p = 0.05, 1 tailed

%%% appearance
% bgColor = 'k';
% lineColor = 'w';
% textColor = 'w';
bgColor = 'w';
lineColor = 'k';
textColor = 'k';

% define field names in Summary struct
% tdr = sprintf('TDRSummary%s',animalName(1));
% otdr = sprintf('oTDRSummary%s',animalName(1));
tdr = 'TDRSummary';
tdrPP = 'TDRSummary_PP';
tdrPRPC = 'TDRSummary_presRespPrevCond';
otdr = 'oTDRSummary';
otdrPP = 'oTDRSummary_PP';
otdrSerial = 'oTDRSummary_serial';
bs = 'BSSummary';
bsPP = 'BSSummary_PP';


%% load data for public release
%%% Assumes that separate .mat files are saved for each output struct per
%%% testTDRoTDR.m with filename: "results_monkey[initial]_[structName].mat"
%%% Each struct will become new field in master struct "Summary"

% set path to results .mat files
dataPath = ''; 
if isempty(dataPath)
    dataPath = input('Enter absolute path to separate data repository: ','s');
end

% path for saving results
dataPath_results = fullfile(dataPath,'results');

% list of structs to load
structs_to_load = {'oTDRSummary','oTDRSummary_serial','oTDRSummary_PP',...
    'TDRSummary','TDRSummary_PP','TDRSummary_presRespPrevCond','BSSummary',...
    'BSSummary_PP','angleProfile'};

% loop through structs 
for s = structs_to_load
    path_temp = fullfile(dataPath_results,sprintf('results_monkey%s_%s.mat',animalName(1),s{:}));
    try
        foo = load(path_temp);
        Summary.(s{:}) = foo.(s{:});
        clear foo
        fprintf('LOADED %s from %s\n',s{:},path_temp);
    catch
        fprintf('FAILED to load %s from %s\n',s{:},path_temp);
    end
end


%% distributions of sRAs

% logical on whether to use previous offer sRAS
bitPP = false;

%%%%%%%%%%%%%%%%%%%%%%
otdrTemp = otdr;
if bitPP
    % for previous offer sRAs
    otdrTemp = [otdrTemp,'_PP'];    
    condLabelTemp = cellfun(@(c) ['Prev ',c],condLabel,'UniformOutput',0);
    analName = 'Previous trial';
else
    condLabelTemp = condLabel;
    analName = 'Present trial';
end

figName = ['Distribution of non-orthogonal sRAs - ',analName];
figure;
whitebg(gcf,bgColor)
set(gcf,'Color',bgColor,'Name',figName);

sRA = Summary.(otdrTemp).sRAStar.RA;
S = size(sRA,2); % number of regression axes
D = size(sRA,1); % number of dimensions (neurons)
B = max(10, round(D/10)); % number of bins
sRAMin = min(sRA(:)); % min sRA coeff values
sRAMax = max(sRA(:)); % max sRA coeff values
bins = linspace(sRAMin,sRAMax,B); % bins
binwidth = diff(bins(1:2));

for i = 1:S
    subplot(1,S,i);
    title(condLabelTemp{Summary.(otdrTemp).sRAStar.param4RA(i)});
    hold on

    plotHist(sRA(:,i),'bins',bins,'testType','ttest','alpha',0.05,...
        'bitPlotBar',1,...
        'bitPlotMean',1,'bit4Paper',1,'nullH',0,'pFontSize',14);
    
    % test if distribution is normally distributed. By using zscore, we
    % transform the data to have mean=0 and std=1.
    [~,P] = kstest(zscore(sRA(:,i)));
    
    % plot KSTest results
    text(1,1,sprintf('P(non-normal) = %0.2g',P),'Units','normalized',...
        'VerticalAlignment','top','HorizontalAlignment','right',...
        'FontSize',14);
    
    % fit guassian
    [mu,sigma] = normfit(sRA(:,i));
    [freq] = hist(sRA(:,i),bins);
    x = linspace(sRAMin,sRAMax+binwidth,100);
    y = normpdf(x,mu,sigma);
    % scale pdf by total obs * bin width
    y = y * size(sRA,1) * diff(bins(1:2));
    % plot
    plot(x,y,'-','Color',textColor,'LineWidth',2);
    
    % appearance
    if i == 1
        ylabel('Number of neurons','FontSize',14);
    end
    if i == ceil(S/2)
        xlabel('sRA coefficient','FontSize',14);
    end
    
    set(gca,'FontSize',14);
end

%%% annotation
figTitleHand = annotation('textbox',[0 0 1 .05]);
    figAnnotStr = {figName
        [datestr(now),', ',mfilename,'.m']
        };    
set(figTitleHand,'String',figAnnotStr,'Interpreter','none','FontSize',12);
set(figTitleHand,'VerticalAlignment','top','HorizontalAlignment','right');
set(figTitleHand,'LineStyle','none');
set(figTitleHand,'Color','r');

clear otdrTemp condLabelTemp

%% pairwise scatter OR histogram of differences of sRA componenet

% set whether to plot scatter (true) vs. difference (false) between pairs
% of sRAs
bitScatter = true;

% set whether to plot absolute value (useful when comparing strength of
% encoding between variables so that the difference A - B is the same for 
% a - b and a - (-b)).
bitAbsMag = false;

% logical on whether to use previous offer sRAS
bitPP = false;
% logical on whether to compare present and previous offer sRAs --
% implies bitPP
bitPresVsPP = false;

%%%%%%%%%%%%%%%%%%

% if comparing present and previous sRAs, implies that previous sRAs are
% included
if bitPresVsPP
    bitPP = true;
end

otdrTemp = otdr;
bsTemp = bs;
if bitPP
    % for previous offer sRAs
    otdrTemp = [otdrTemp,'_PP'];
    bsTemp = [bsTemp,'_PP'];
    analName = 'Previous trial';
else
    analName = 'Present trial';
end

% if comparing present and previous sRAs, include present offer sRAs:
if bitPresVsPP 
    sig = Summary.(otdr).sRAStar.RA;
    condLabelTemp = condLabel;
    % include all pairs of present vs previous sRAs for comparison
%         sig2plot = combvec([1:size(sig,2)],[1:size(Summary.(otdrTemp).sRAStar.RA,2)]+size(sig,2))';
    % (below code is an alternative to MATLAB's combvec(), which is
    % only available with the Deep Learning Toolbox)
    [A,B] = meshgrid([1:size(sig,2)], [1:size(Summary.(otdrTemp).sRAStar.RA,2)]+size(sig,2));
    c = cat(2,A',B');
    sig2plot = reshape(c,[],2);
    clear A B c

    % compute bootstrapped mean and std for present offers
    if isfield(Summary,bs)
        sRA_mean = squeeze(mean(Summary.(bs).sRA,2))';
        sRA_std = squeeze(std(Summary.(bs).sRA,[],2))';
    else
        sRA_mean = [];
        sRA_std = [];
    end

else
    % comparisons only within an epoch (present or previous)
    sig2plot = [1 2; 1 3; 2 3];
    sig = [];
    condLabelTemp = {};
    sRA_mean = [];
    sRA_std = [];
end

% append sRAs (may be appending to empty array when not comparing
% between epochs
sig = [sig, Summary.(otdrTemp).sRAStar.RA];
if bitPP
    condLabelTemp = [condLabelTemp, cellfun(@(c) ['Prev ',c],condLabel,'UniformOutput',0)];
else
    condLabelTemp = condLabel;
end

% compute and append bootstrapped mean and std
if isfield(Summary,bsTemp)
    sRA_mean = [sRA_mean, squeeze(mean(Summary.(bsTemp).sRA,2))'];
    sRA_std = [sRA_std, squeeze(std(Summary.(bsTemp).sRA,[],2))'];
end           
    

% if isfield(Summary,bs)
%     sRA_mean = squeeze(mean(Summary.(bs).sRA,2))';
%     sRA_std = squeeze(std(Summary.(bs).sRA,[],2))';
%     
%     % add previous offer if requested
%     if bitPP && isfield(Summary,bsTemp)
%         sRA_mean = [sRA_mean, squeeze(mean(Summary.(bsTemp).sRA,2))'];
%         sRA_std = [sRA_std, squeeze(std(Summary.(bsTemp).sRA,[],2))'];
%     end   
% else
%     sRA_mean = [];
%     sRA_std = [];
% end

S = size(sig,2); % number of signals
D = size(sig,1); % number of dimensions (neurons)
B = max(10, round(D/10)); % number of bins
nAx = size(sig2plot,1);

% intantiate to store differences
d = cell(size(sig2plot,1),1);
dTitle = d; % to store name of difference

if ~bitScatter && ~bitAbsMag
    warning('When plotting distributions of differences, you should generally compare the absolute strength of encoding using bitAbsMag');
end

figure;
whitebg(gcf,bgColor)
if bitScatter
    figType = 'scatter';
else
    figType = 'difference';
end
figName = sprintf('sRA pairwise %s - %s',figType,analName);
set(gcf,'Color',bgColor,'Name',figName);

if bitAbsMag
    genericLabel = '|%s|';
else
    genericLabel = '%s';
end

% determine subplot layout
if nAx <= 3
    % bias toward 1 row
    nRow = 1;
    nCol = nAx;
else
    foo = get(gca,'Position');
    [nRow,nCol] = idealSubplotDim(nAx,foo(4)/foo(3),false);
    clear foo
end

for i = 1:nAx
    strX = condLabelTemp{sig2plot(i,1)};
    strY = condLabelTemp{sig2plot(i,2)};
    
    subplot(nRow,nCol,i);
    
    x = sig(:,sig2plot(i,1));
    y = sig(:,sig2plot(i,2));
    foo = 'corrPear';
    if bitAbsMag
        x = abs(x);
        y = abs(y);
        foo = [];
    end
    
    % store difference
    d{i} = x - y;
    % name of difference
    if bitAbsMag
        dTitle{i} = sprintf('|%s| - |%s|',strX,strY);
    else
        dTitle{i} = sprintf('%s - %s',strX,strY);
    end
    
    if bitScatter
        [h,~,~,~,hT] = plotScatterN(x,y,'bitEqualAxes',true,...
            'statTest',foo,'alpha',0.05,'pTextFontSize',14,'color',textColor);
        set(h,'Marker','o','MarkerSize',5);
        set(hT,'Position',[0.02 1],'Units','normalized',...
            'VerticalAlignment','top')
        hold on
        
        % add angle text 
        
        % index of signals within a set of sRAs (i.e., present vs present
        % or previous vs previous).
        goo = [1 2; 1 3; 2 3]; % order of signal pairs in .sRAStar.angleAnalysis
        if all(sig2plot(i,:) >= 4)
            goo = goo + 3;
        end
        [posI,posJ] = find(ismember(goo,sig2plot(i,:),'rows'));
        if ~isempty(posI)
            % set field to use
            angleField2Use = 'angleAnalysis';
        elseif bitPresVsPP
            
            % match with index of signals between sets of sRAs (i.e.,
            % present vs. previous). Because sig2Plot refers to present
            % offer in col 1 and previous offer in col 2 and yet matricies
            % in angleAnalysis_withVect4Angle are arranged with previous
            % offer in rows and present offer in columns, we have to
            % reverse the order of the columns. Also, we have to re-index
            % previous offer, since, in the matricies, it occupies its own
            % dimension (not concatenated to the end of present offer).
            posI = sig2plot(i,2) - size(Summary.(otdr).sRAStar.RA,2); % previous offer, re-indexed
            posJ = sig2plot(i,1); % present offer
            
            % set field to use
            angleField2Use = 'angleAnalysis_withVect4Angle';
            
        end
        
        % if matching angle data found, print text:
        if ~isempty(posI)
            text(1,1,sprintf('folded \\theta=%0.0f\\circ\np(H_1:\\theta=0\\circ)=%0.2g, p(H_1:\\theta=90\\circ)=%0.2g',...
                Summary.(otdrTemp).sRAStar.(angleField2Use).angleFolded(posI,posJ),...
                Summary.(otdrTemp).sRAStar.(angleField2Use).pAngleFolded(posI,posJ),...
                Summary.(otdrTemp).sRAStar.(angleField2Use).pAngleFolded_large(posI,posJ)),...
                'Units','normalized','FontSize',14,'Color',textColor,...
                'HorizontalAlignment','right','VerticalAlignment','top');
        end        
        clear foo goo
        
        %%% error bars showing std dev of estimates
        if ~isempty(sRA_mean)
            % horizontal bars
            xH = [x - sRA_std(:,sig2plot(i,1)), x + sRA_std(:,sig2plot(i,1))]';
            yH = repmat(y',2,1);
            line(xH,yH,'Color',textColor)
            
            xV = repmat(x',2,1);
            yV = [y - sRA_std(:,sig2plot(i,2)), y + sRA_std(:,sig2plot(i,2))]';
            line(xV,yV,'Color',textColor)
            
            strX = sprintf('%s (\\pm s.d.)',sprintf(genericLabel,strX));
            strY = sprintf('%s (\\pm s.d.)',sprintf(genericLabel,strY));
            
%             % fill in if x-var is significant
%             bitSig = Summary.(bsTemp).sRA_p(:,sig2plot(i,1))<0.05;
%             hFill = plot(x(bitSig),y(bitSig),'o','color',textColor,...
%                 'MarkerSize',5,'MarkerFaceColor',textColor);
            
%             % fill in if y-var is significant
%             bitSig = Summary.(bsTemp).sRA_p(:,sig2plot(i,2))<0.05;
%             hFill = plot(x(bitSig),y(bitSig),'o','color',textColor,...
%                 'MarkerSize',5,'MarkerFaceColor',textColor);
        end
        
        
        %%% horizontal and vertical meridians
        xLim = get(gca,'XLim');
        yLim = get(gca,'YLim');
        line([0 0],yLim,'Color',textColor','LineStyle','--');
        line(xLim,[0 0],'Color',textColor','LineStyle','--');
        
        %%% Labels
        xlabel(strX,'FontSize',16);
        ylabel(strY,'FontSize',16);
        set(gca,'FontSize',16);
    else
        % PLOT HISTOGRAM OF DIFFERENCE
        plotHist(x - y,'nBin',B,'testType','ttest','alpha',0.05,...
            'bitPlotBar',1,...
            'bitPlotMean',1,'bit4Paper',1,'nullH',0,'pFontSize',14);
        
        % appearance
        title(dTitle{i});
        if i == 1
            ylabel('Number of neurons','FontSize',16);
        end
        if i == ceil(S/2)
            xlabel('Difference of sRA coefficient','FontSize',16);
        end
        set(gca,'FontSize',16,'box','off');

        
    end
    clear x y
end

% add annotation
dt = round(diff(Summary.(otdrTemp).sRA.t4RA{1}(1:2)),2);
figTitleHand = annotation('textbox',[0 0 1 .05]);
if isfield(Summary.(otdrTemp).sRA,'magThresh') && ~isempty(Summary.(otdrTemp).sRA.magThresh)
    figAnnotStr = {sprintf('%s of sRA components, present trial for times where R2 > %g',...
        figType,Summary.(otdrTemp).sRA.magThresh);
        [datestr(now),', ',mfilename,'.m']
        };    
else
    figAnnotStr = {sprintf('%s of sRA components, present trial for times: Benefit [%g %g], Choice [%g %g], Exp Rwd [%g %g]',...
        figType,Summary.(otdrTemp).sRA.t4RA{1}(1)-dt/2,Summary.(otdrTemp).sRA.t4RA{1}(end)+dt/2,Summary.(otdrTemp).sRA.t4RA{2}(1)-dt/2,Summary.(otdrTemp).sRA.t4RA{2}(end)+dt/2,...
        Summary.(otdrTemp).sRA.t4RA{3}(1)-dt/2,Summary.(otdrTemp).sRA.t4RA{3}(end)+dt/2);
%         Summary.(otdrTemp).sRA.t1(1)-dt/2,Summary.(otdrTemp).sRA.t1(end)+dt/2,Summary.(otdrTemp).sRA.t2(1)-dt/2,Summary.(otdrTemp).sRA.t2(end)+dt/2,...
%         Summary.(otdrTemp).sRA.t3(1)-dt/2,Summary.(otdrTemp).sRA.t3(end)+dt/2);
        [datestr(now),', ',mfilename,'.m']
        };
end
set(figTitleHand,'String',figAnnotStr,'Interpreter','none','FontSize',12);
set(figTitleHand,'VerticalAlignment','top','HorizontalAlignment','right');
set(figTitleHand,'LineStyle','none');
set(figTitleHand,'Color','r');

%%% Compute difference of differences

% instantiate
dd = NaN(nchoosek(size(d,1),2),2);

% loop through all pairs
ddN = 0;
for i = 1:size(d,1)-1
    for j = i+1:size(d,1);
        ddN = ddN + 1;
%         % compute mean and p-value difference of differences
%         dd(ddN,1) = mean(d{i} - d{j});
%         [~,dd(ddN,2)] = ttest(d{i} - d{j});

        % compute median and p-value difference of differences
        dd(ddN,1) = median(d{i} - d{j});
        dd(ddN,2) = signrank(d{i} - d{j});
    
        % display
        fprintf('\nDifference of differences dd = (%s) - (%s):\n',dTitle{i},dTitle{j})
        fprintf('... median(dd) = %0.2f, p(median(dd) == 0) = %0.2f\n',dd(ddN,1),dd(ddN,2))
        
    end
end
clear otdrTemp bsTemp condLabelTemp

%% sRA angles 

% logical on whether to use previous offer sRAS
bitPP = false;

%%%%%%%%%%%%%%%%%%%%%%
otdrTemp = otdr;
if bitPP
    % for previous offer sRAs
    otdrTemp = [otdrTemp,'_PP'];    
    condLabelTemp = cellfun(@(c) ['Prev ',c],condLabel,'UniformOutput',0);
    analName = 'Previous trial';
else
    condLabelTemp = condLabel;
    analName = 'Present trial';
end

figure;
whitebg(gcf,bgColor)
figName = ['sRA angles - ',analName];
set(gcf,'Color',bgColor,'Name',figName);

nRA = size(Summary.(otdrTemp).sRA.RA,2); % number of regression axes

% build labels
crossLabel = {};
count = 0;
for i = 1:nRA-1
    for j = i+1:nRA
        count = count + 1;
        crossLabel{count} = sprintf('%s x %s',...
            condLabelTemp{Summary.(otdrTemp).sRAStar.param4RA(i)},...
            condLabelTemp{Summary.(otdrTemp).sRAStar.param4RA(j)});
    end
end


% unfolded angles
subplot(1,2,1)
% center at 90
h = bar(90 - Summary.(otdrTemp).sRAStar.angleAnalysis.angle);
% bar colors
set(h(1),'FaceColor','none','EdgeColor',lineColor)
% balance range
yLim = get(gca,'YLim');
ylim(max(abs(yLim)) * [-1 1]);
yLim = get(gca,'YLim');

% relabel to center at 90
yTickLabel = get(gca,'YTickLabel');
set(gca,'YTickLabel',num2str(cellfun(@str2num,yTickLabel)+90));
ylabel('Angle between sRA (deg)','FontSize',16)

% appearance
set(gca,'XTick',[1:length(crossLabel)],'XTickLabel',crossLabel)
set(gca,'FontSize',12,'box','off');
title('Un-folded between-signal angles','FontSize',16);
% add p-values
for i = 1:length(Summary.(otdrTemp).sRAStar.angleAnalysis.pAngle)
    text(i,yLim(2),sprintf('P(angle > 90 | \nangle < 90) = %0.2g',Summary.(otdrTemp).sRAStar.angleAnalysis.pAngle(i)),...
        'VerticalAlignment','top','HorizontalAlignment','center','FontSize',10);
end

%%%%
% FOLDED angles
subplot(1,2,2)
% center at 90
h = bar(Summary.(otdrTemp).sRAStar.angleAnalysis.angleFolded - 90);
% bar colors
set(h(1),'FaceColor','none','EdgeColor',lineColor)
yLim = get(gca,'YLim');
% ylim([ yLim(2) + 5]);
yLim = get(gca,'YLim');
ylabel('Angle between sRA (deg)','FontSize',16)

% relabel to center at 90
yTickLabel = get(gca,'YTickLabel');
set(gca,'YTickLabel',num2str(cellfun(@str2num,yTickLabel)+90));
ylabel('Angle between sRA (deg)','FontSize',16)

% appearance
set(gca,'XTick',[1:length(crossLabel)],'XTickLabel',crossLabel)
set(gca,'FontSize',12,'box','off');
title('Folded between-signal angles','FontSize',16);
% add p-values
for i = 1:length(Summary.(otdrTemp).sRAStar.angleAnalysis.pAngleFolded)
    text(i,yLim(2),sprintf('p(angle > 0) = %0.2g \np(angle < 90) = %0.2g',...
        Summary.(otdrTemp).sRAStar.angleAnalysis.pAngleFolded(i),...
        Summary.(otdrTemp).sRAStar.angleAnalysis.pAngleFolded_large(i)),...
        'VerticalAlignment','top','HorizontalAlignment','center','FontSize',10);
end

%%% annotation
figTitleHand = annotation('textbox',[0 0 1 .05]);
    figAnnotStr = {figName
        [datestr(now),', ',mfilename,'.m']
        };    
set(figTitleHand,'String',figAnnotStr,'Interpreter','none','FontSize',12);
set(figTitleHand,'VerticalAlignment','top','HorizontalAlignment','right');
set(figTitleHand,'LineStyle','none');
set(figTitleHand,'Color','r');

clear otdrTemp condLabelTemp

%% comparing sRA values to bootstrapped values 

% logical on whether to use previous offer sRAS
bitPP = false;

%%%%%%%%%%%%%%%%%%%%%%
otdrTemp = otdr;
bsTemp = bs;
if bitPP
    % for previous offer sRAs
    otdrTemp = [otdrTemp,'_PP'];    
    bsTemp = [bsTemp,'_PP'];    
    condLabelTemp = cellfun(@(c) ['Prev ',c],condLabel,'UniformOutput',0);
else
    condLabelTemp = condLabel;
end

% compute mean and std of bootstrapped values
sRA_mean = squeeze(mean(Summary.(bsTemp).sRA,2))';
sRA_std = squeeze(std(Summary.(bsTemp).sRA,[],2))';

% plot scatter of average bootstrapped value to full-trial value
figName = 'Bootstrapped sRA vs. full-trial sRA coefficients';
figure('Name',figName);
whitebg(gcf,bgColor);
set(gcf,'Color',bgColor);

nSRA = size(sRA_mean,2);
for i=1:nSRA
    subplot(1,nSRA,i);
    [hS,~,~,~,hT] = plotScatterN(Summary.(otdrTemp).sRAStar.RA(:,i),sRA_mean(:,i),...
        'bitEqualAxes',true,...
        'statTest','ttest','alpha',0.05,'pTextFontSize',14,'color',textColor);
    set(hS,'Marker','o','MarkerSize',7);
    set(gca,'FontSize',16);
    title(condLabelTemp{i})
    if i==1
        ylabel('sRA bootstrapped coefficient (\pm s.d.)','FontSize',18);
    end
    if i == ceil(nSRA/2)
        xlabel('sRA full-trial coefficient','FontSize',18);
    end
    set(hT,'Position',[0.02 1],'Units','normalized',...
        'VerticalAlignment','middle')
    
    
    % add error bars representing STD of bootstrapped value
    x = repmat(Summary.(otdrTemp).sRAStar.RA(:,i)',2,1);
    y = [sRA_mean(:,i) + sRA_std(:,i), sRA_mean(:,i) - sRA_std(:,i)]';
    h = line(x,y,'Color',textColor);
end

clear otdrTemp bsTemp condLabelTemp

%% confusion maxtrix of significant sRA value

%%% Computes number of neurons with significant encoding for pairwise
%%% conjunction of variables i and j, where signifiance is determined as
%%% the 2-tailed probability that the distribution of bootstrapped sRA
%%% coefficients includes 0 is less than alpha

% set p-value threshold for significant encoding
alphaSRA = 0.05;

% logical on whether to use previous offer sRAS
bitPP = false;

%%%%%%%%%%%%%%%%%%%%%%
otdrTemp = otdr;
bsTemp = bs;
if bitPP
    % for previous offer sRAs
    otdrTemp = [otdrTemp,'_PP'];    
    bsTemp = [bsTemp,'_PP'];    
    condLabelTemp = cellfun(@(c) ['Prev ',c],condLabel,'UniformOutput',0);
else
    condLabelTemp = condLabel;
end


[~,S] = size(Summary.(bsTemp).sRA_p); % number of neurons and sRAs

% extract sRA p-value
sRA_p = Summary.(bsTemp).sRA_p;
% eliminate any neuron that does not have p-defined for all sRAs
sRA_p(~all(isfinite(sRA_p),2),:) = [];
numNeurons = size(sRA_p,1);

% determine signficance for individual neurons and sRAs
bitSig = sRA_p < alphaSRA;
% compute total number of significant per sRA
nSig = sum(bitSig);

% build index of pair of sRAs
nPair = nchoosek(S,2);
ind = NaN(nPair,2);
foo = 0;
for i = 1:S-1
    for j = i+1:S
        foo = foo+1;
        ind(foo,:) = [i j];
    end
end

% instantiate confusion matrix (pairs x 2 x 2)
CM = NaN(nPair,2,2);
CM_exp = CM; % expected
CM_simple = NaN(nPair,2);
CM_simple_exp = NaN(nPair,2);

% loop through pairs computing CM for each
for i = 1:nPair
    % jointly significant
    CM(i,1,1) = sum(all(bitSig(:,ind(i,:)),2));
    
    % jointly non-significant
    CM(i,2,2) = sum(all(~bitSig(:,ind(i,:)),2));
    
    % significant for row, but not column
    CM(i,1,2) = sum(bitSig(:,ind(i,1)) & ~bitSig(:,ind(i,2)));

    % significant for column, but not row
    CM(i,2,1) = sum(bitSig(:,ind(i,2)) & ~bitSig(:,ind(i,1)));
    
    % confirm that marginals are as expected
    if sum(CM(i,1,:)) ~= nSig(ind(i,1))
        error('Row marginal (=%d) did not match number signficant (=%d) for sRA #%d',sum(CM(i,1,:)),nSig(ind(i,1)),ind(i,1))
    end
    if sum(CM(i,2,:)) ~= numNeurons-nSig(ind(i,1))
        error('Row marginal (=%d) did not match number NON-signficant (=%d) for sRA #%d',sum(CM(i,2,:)),numNeurons-nSig(ind(i,1)),ind(i,1))
    end
    if sum(CM(i,:,1)) ~= nSig(ind(i,2))
        error('Column marginal (=%d) did not match number signficant (=%d) for sRA #%d',sum(CM(i,:,1)),nSig(ind(i,2)),ind(i,2))
    end
    if sum(CM(i,:,2)) ~= numNeurons-nSig(ind(i,2))
        error('Column marginal (=%d) did not match number NON-signficant (=%d) for sRA #%d',sum(CM(i,:,2)),numNeurons-nSig(ind(i,2)),ind(i,2))
    end
    
    % compute expected
    row_pSig = sum(CM(i,1,:)) / numNeurons;
    row_pNS = sum(CM(i,2,:)) / numNeurons;
    % confirm
    if row_pSig + row_pNS - 1 > eps
        error('proportion sig and non-sig should sum to 1')
    end
    CM_exp(i,:,:) = repmat([nSig(ind(i,2)) numNeurons-nSig(ind(i,2))],2,1) .* repmat([row_pSig; row_pNS],1,2);
    
    %%% 
    %%% compute simplified stats on jointly significant vs. not jointly significant
    %%%
    CM_simple(i,1) = CM(i,1,1); % jointly significant
    CM_simple(i,2) = numNeurons - CM(i,1,1); % not jointly significant
    
    CM_simple_exp(i,1) = CM_exp(i,1,1);
    CM_simple_exp(i,2) = numNeurons - CM_exp(i,1,1);
end

% compute chi2 on (obs-exp)^2 / exp
foo = (CM - CM_exp).^2 ./ CM_exp;
% sum within each pair
chi2 = NaN(nPair,1);
chi2_p = chi2;
for i = 1:nPair
    chi2(i) = sum(sum(foo(i,:,:)));
    % compute chi2 p-value based on 1 df
    chi2_p(i) = chi2cdf(chi2(i),1,'upper');
end

% I'm not sure this is kosher, but we will attempt to compute an omnibus
% p-value by summing chi2 across all pairs and, because we're treating each
% pair independently, compute a grand chi2 based on the df's are equal to
% the number of pairs - 1
chi2_p_grand = chi2cdf(sum(chi2),nPair-1,'upper');

%%%
%%% compute simplified stats on jointly significant vs. not jointly significant
%%%
% compute chi2 on (obs-exp)^2 / exp
foo = (CM_simple - CM_simple_exp).^2 ./ CM_simple_exp;
% sum within each pair
chi2_simple = NaN(nPair,1);
chi2_simple_p = chi2_simple;
for i = 1:nPair
    chi2_simple(i) = sum(sum(foo(i,:,:)));
    % compute chi2 p-value based on 1 df
    chi2_simple_p(i) = chi2cdf(chi2_simple(i),1,'upper');
end

% grand p-value
chi2_simple_p_grand = chi2cdf(sum(chi2_simple),nPair-1,'upper');
%%%


% display
fprintf('\n\n****************\n')
for i = 1:S
    fprintf('%s: %d of %d (%0.0f%%) significant coefficients \n',condLabelTemp{i},nSig(i),numNeurons,100*nSig(i)/numNeurons)
end
for i = 1:nPair
    fprintf('\nFOR %s x %s:\n',condLabelTemp{ind(i,1)},condLabelTemp{ind(i,2)})
    fprintf('Number of neurons with significant / non-significant coefficient (p < %0.2f):\n',alphaSRA)
    disp(squeeze(CM(i,:,:)))
    fprintf('Proportion of neurons with significant / non-significant coefficient (p < %0.2f):\n',alphaSRA)
    disp(squeeze(CM(i,:,:))/numNeurons);

    fprintf('EXPECTED Number of neurons with significant / non-significant coefficient (p < %0.2f):\n',alphaSRA)
    disp(squeeze(CM_exp(i,:,:)))
    fprintf('EXPECTED Proportion of neurons with significant / non-significant coefficient (p < %0.2f):\n',alphaSRA)
    disp(squeeze(CM_exp(i,:,:))/numNeurons);

    fprintf('Chi2 stat = %0.2f, df = %d\n',chi2(i),1)
    fprintf('Probability that joint significance was different than independent assortment = %0.2g\n',chi2_p(i))
end

fprintf('\nGRAND probability across pairs that joint significance was different than independent assortment = %0.2g\n',chi2_p_grand)

%%%
%%%  simplified stats on jointly significant vs. not jointly significant
%%%
fprintf('\n\n****************\n')
fprintf('Simplified stats based on number of jointly significant vs. non-jointly significant (not further separated by jointly non-sig, sig for 1 and not other, etc')
fprintf('Number of neurons by variable (rows) with significant / non-significant coefficient (columns) (p < %0.2f):\n',alphaSRA)
disp(CM_simple)
fprintf('Proportion of neurons by variable (rows) with significant / non-significant coefficient (columns) (p < %0.2f):\n',alphaSRA)
disp(CM_simple/numNeurons);

fprintf('EXPECTED Number of neurons by variable (rows) with significant / non-significant coefficient (p < %0.2f):\n',alphaSRA)
disp(CM_simple_exp)
fprintf('EXPECTED Proportion of neurons by variable (rows) with significant / non-significant coefficient (p < %0.2f):\n',alphaSRA)
disp(CM_simple_exp/numNeurons);

fprintf('Chi2 stat on jointly vs non-jointly significant:\n')
disp(chi2_simple)
fprintf('Probability that joint vs non-joint significance was different than independent assortment:\n')
disp(chi2_simple_p)

fprintf('GRAND probability across pairs that joint vs. non-joint significance was different than independent assortment = %0.2g\n',chi2_simple_p_grand)

clear condLabelTemp


%% Statistics of within- and between-sRA correlations

% logical on whether to use previous offer sRAS
bitPP = false;
 
%%%%%%%%%%%%%%%%%%%%%%
bsTemp = bs;

if bitPP
    % for previous offer sRAs
    bsTemp = [bsTemp,'_PP'];    
    condLabelTemp = cellfun(@(c) ['Prev ',c],condLabel,'UniformOutput',0);
else
    condLabelTemp = condLabel;
end

% for convenience
s = Summary.(bsTemp);

%%%% PLOTTING
figName = 'sRA correlations between and within signals - HISTOGRAMS';
figure('Name',figName);
whitebg(gcf,bgColor)
set(gcf,'Color',bgColor);

nSignal = size(s.corrWithin,1);
nSPair = size(s.corrBetween,1);
nAx = nSignal*2;
nCol = nSignal;
nRow = 2;

%%%% WITHIN SIGNALS
for i = 1:nSignal
    subplot(nRow,nCol,i)
    histogram(s.corrWithin(i,:),100,'EdgeAlpha',0,'FaceColor',[0.4 0.4 0.4]);
    title(condLabelTemp{i},'FontSize',16);
    xlim([0 1]);
    set(gca,'FontSize',16,'box','off');
    ylabel('Count','FontSize',16);

    text(0.05,0.95,sprintf('mean = %0.2f, p = %g',s.stats.within.r(i),s.stats.within.p(i)),'Units','normalized','FontSize',16);
end

%%%% BETWEEN SIGNALS
for i = 1:nSPair
    subplot(nRow,nCol,i+nSignal)
    condN = s.corrBetweenInd(i,:);    
    
    %%% Instead of plotting a scalar value for predicted r
    %%% (s.stats.nullModel.rPred) and a distribution of boostrapped
    %%% empirical r (s.stats.between.r), we now plot a distribution of
    %%% predicted r (s.stats.nullModel.rPred_dist) and a scalar empirical r
    %%% as measured from the all-trials sRAa (s.stats.between.r_allTrials).

    % store sign of observed correlation:
    rSign = sign(s.stats.between.r_allTrials(i));
    if rSign == 0
        error('sign of between signal correlation not expected to = 0')
    end
    
    %%% NEW WAY based on distribution of predicted R
    histogram(rSign .* s.stats.nullModel.rPred_dist(i,:),100,'EdgeAlpha',0,'FaceColor',[0.4 0.4 0.4]);
%     foo = mean(rSign .* s.stats.nullModel.rPred_dist(i,:));
    foo = rSign * s.stats.nullModel.rPred_dist_mean(i);
    yLim = get(gca,'YLim');
    text(1,1,...
        sprintf('pred r = %0.2g, p(obs-pred) = %0.2g',...
        foo, s.stats.nullModel.pVal_alt_ttest(i)),...
        'HorizontalAlignment','right','VerticalAlignment','top',...
        'Units','normalized','FontSize',12,'Color','k');
    
    % Plot the measured R value from sRAStar
    line(s.stats.between.r_allTrials(i)*[1 1],[yLim(1), 0.9*range(yLim)+yLim(1)],'LineStyle','--','Color','r');
    text(0.02,0.85,sprintf('r_{all} = %0.2g, p(mu = 0) = %0.2g',...
        (s.stats.between.r_allTrials(i)),s.stats.between.p_allTrials(i)),...
        'Units','normalized','FontSize',12,'Color','r');        
    
    xlim(sort([0 rSign]));    
    
    %%% OLD WAY based on scalar predicted R
%     hist(s.corrBetween(i,:),10);
%     text(0.02,0.92,sprintf('mean = %0.2f, p(mu = 0) = %g',...
%         s.stats.between.r(i),s.stats.between.p(i)),'Units','normalized','FontSize',12);    
%     
%     %%% Plot the expected r value based on perfect correlation corrupted by
%     %%% noise
%     yLim = get(gca,'YLim');
%     line(s.stats.nullModel.rPred(i)*[1 1],yLim,'LineStyle','--','Color','r');
%     text(s.stats.nullModel.rPred(i),yLim(2),...
%         sprintf('pred r = %0.2g, p(obs-pred) = %0.2g',...
%         s.stats.nullModel.rPred(i),...
%         s.stats.nullModel.pVal(i)),...
%         'HorizontalAlignment','right','VerticalAlignment','top',...
%         'FontSize',12,'Color','r');
%     xlim([-1 1]);    

    %%% APPEARANCE
    set(gca,'FontSize',16,'box','off');   
    s.stats.between.condLabelTemp{i} = sprintf('%s x %s',condLabelTemp{condN(1)},condLabelTemp{condN(2)});
    title(s.stats.between.condLabelTemp{i},'FontSize',16);
    ylabel('Count','FontSize',16);
    if i == ceil(nSignal/2)
        xlabel('Pearson''s r','FontSize',16);
    end
    
end
%%% annotation
figTitleHand = annotation('textbox',[0 0 1 .05]);
    figAnnotStr = {figName
        [datestr(now),', ',mfilename,'.m']
        };    
set(figTitleHand,'String',figAnnotStr,'Interpreter','none','FontSize',12);
set(figTitleHand,'VerticalAlignment','top','HorizontalAlignment','right');
set(figTitleHand,'LineStyle','none');
set(figTitleHand,'Color','r');




%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUMMARY
figName = 'sRA correlations between and within signals - SUMMARY';
figure('Name',figName);
whitebg(gcf,bgColor)
set(gcf,'Color',bgColor);

%%% BETWEEN SIGNAL
subplot(1,2,1)

% %%%%%% BASED ON MEAN/CI OF BOOTSTRAPPED AND SCALAR VALUE OF R_PRED
% % plot
% errTemp = [s.stats.between.r'-s.stats.between.r_ci(:,1), s.stats.between.r_ci(:,2)-s.stats.between.r'];
% h = barwitherr(errTemp,s.stats.between.r','LineWidth',2);
% % bar colors
% set(h(1),'FaceColor','none','EdgeColor',lineColor)
% 
% %%%
% % add horizontal line showing predicted correlation under "perfect
% % correlation corrupted by noise" model
% % Get bar location and width
% barX = get(h(1),'XData');
% barW = get(h(1),'BarWidth');
% % loop through each bar
% for i = 1:length(barX)
%     line([-barW/2, barW/2]'+barX(i),...
%         mean(s.stats.nullModel.rPred_dist(i,:))*sign(s.stats.between.r(i)) * [1 1]',...
%         'Color',textColor,'LineStyle','--');
% end

%%%%%% BASED ON R_ALLTRIALS AND MEAN/CI OF BOOTSTRAPPED R_PRED
% plot
h = bar(s.stats.between.r_allTrials,'LineWidth',2);
hold on
% bar colors
set(h(1),'FaceColor','none','EdgeColor',lineColor)

%%%
% add shadded area showing CI of predicted correlation under "perfect
% correlation corrupted by noise" model
% Get bar location and width
barX = get(h(1),'XData');
barW = get(h(1),'BarWidth');

% loop through each bar
for i = 1:length(barX)
    fooSign = sign(s.stats.between.r_allTrials(i));
    fooY = s.stats.nullModel.rPred_dist_CI(i,:) * fooSign;
    fooY = [fooY, fliplr(fooY)];
    fooX = [-barW/2, -barW/2, barW/2, barW/2] + barX(i);
    fill(fooX,fooY,[0.7 0.7 0.7],'EdgeColor','none');
    % line for mean
    line([-barW/2, barW/2]' + barX(i),...
        s.stats.nullModel.rPred_dist_mean(i) * [1 1] * fooSign,...
        'Color',[0.4 0.4 0.4],'LineWidth',2,'LineStyle','--');
end

% appearance
set(gca,'XTick',[1:nSignal],'XTickLabel',s.stats.between.condLabelTemp,'XTickLabelRotation',-45)
ylabel(sprintf('Pearson r, %0.1f to %0.1f%% CI',s.stats.between.ci(1),...
    s.stats.between.ci(2)),'FontSize',16);
ylim([-1 1]);
set(gca,'FontSize',16,'box','off');
title('Between-signal correlation','FontSize',16);

% pad either side of the bar plot
xLim = get(gca,'XLim');
xlim([xLim(1)-0.5 xLim(2)+0.5]);

% move bar series to top
foo = get(gca,'Children');
set(gca,'Children',circshift(foo,1));

clear foo fooX fooY fooSign

%%% WITHIN SIGNAL
subplot(1,2,2)

% plot
errTemp = [s.stats.within.r'-s.stats.within.r_ci(:,1), s.stats.within.r_ci(:,2)-s.stats.within.r'];
h = barwitherr(errTemp,s.stats.within.r','LineWidth',2);
% bar colors
set(h(1),'FaceColor',[0.7 0.7 0.7])

% appearance
set(gca,'XTick',[1:nSignal],'XTickLabel',condLabelTemp,'XTickLabelRotation',-45)
ylabel(sprintf('Pearson r, %0.1f to %0.1f%% CI',s.stats.within.ci(1),...
    s.stats.within.ci(2)),'FontSize',16);
ylim([-1 1]);
set(gca,'FontSize',16,'box','off');
title('Within-signal reliability','FontSize',16);

% pad either side of the bar plot
xLim = get(gca,'XLim');
xlim([xLim(1)-0.5 xLim(2)+0.5]);

%%% annotation
figTitleHand = annotation('textbox',[0 0 1 .05]);
    figAnnotStr = {figName
        [datestr(now),', ',mfilename,'.m']
        };    
set(figTitleHand,'String',figAnnotStr,'Interpreter','none','FontSize',12);
set(figTitleHand,'VerticalAlignment','top','HorizontalAlignment','right');
set(figTitleHand,'LineStyle','none');
set(figTitleHand,'Color','r');

clear s bsTemp condLabelTemp


%% oTDR Variance explained by sRAs

% This section plots time courses of variance explained. To recreate the
% figures in Kimmel et al., Nat Comm, 2000, use the following parameter
% values:
%   Figure 4c,d and 6a,b (right panels): 
%         eN = 1;
%         var2Plot =  {'RSV','ISV'}; 
%         varType = 2;
%         bitSumDim = false; 
%         bitMatchSignal2SRA = true;
%         bitUsePC = false;
%         bitUsePC_CC = false;
%         bitSepAx = false;
%         bitTimeOfSRA = true;
%         bitTimeOfEvent = false;
%         bitRSVPartial = true;
%         maxVar2Plot = 'none';
%         bitRSVBound = false;
%         bitRSVBoundDiff = false;
%         fn = otdr;
%         totVarFN = 'totalVar_t';
%         bitLine4PThresh = true;
%         bit4EPS = false; 
%   Figure 4e,f and 6c,d (right panels) -- same as Figure 4c,d EXCEPT:
%         varType = 4;
%   Figure 6a,b (left panels) -- same as Figure 4c,d EXCEPT:
%         eN = 2;
%         var2Plot =  {'RSV'};
%         bitTimeOfEvent = true;
%   Figure 6c,e (left panels) -- same as Figure 4e,f EXCEPT:
%         eN = 2;
%         var2Plot =  {'RSV'};
%   Figure 6e,f (left panels) -- same as Figure 4e,f EXCEPT:
%         eN = 3;
%         var2Plot =  {'RSV'};
%         bitTimeOfEvent = true;
%   Figure 6e,f (right panels) -- same as Figure 4e,f EXCEPT:
%         eN = 5;
%         var2Plot =  {'RSV'};
%   Figure 6g,h (left panels) -- same as Figure 6e,f (left panels) EXCEPT:
%         varType = 4;
%   Figure 6g,h (right panels) -- same as Figure 6e,f (right panels) EXCEPT:
%         varType = 4;
%   Supp Fig 12 -- same as Figure 4c-f, but with altnerative
%         implementations of oTDR.
%   Supp Fig 13a,b -- same as Figure 4c,d EXCEPT
%         var2Plot =  {'V'};
%         bitLine4PThresh = false;
%         maxVar2Plot = 'dPC';
%   Supp Fig 13c,d -- same as Supp Fig 13a,b EXCEPT
%         varType = 4;
%         maxVar2Plot = 'none';
%   Supp Fig 14a,b -- Same as Figure 4c,d EXCEPT
%         var2Plot =  {'RSV'};
%         bitMatchSignal2SRA = false;
%         bitRSVBound = true;
%         bitSepAx = true;
%         (Also, colored dashed lines manually removed and line for
%         on-target variable thickened manually)
%   Supp Fig 14c,e -- Same as Supp Fig 14a,b EXCEPT
%         varType = 4;
%   Supp Fig 16a,b -- Same as Supp Fig 13a,b EXCEPT
%         bitUsePC = true;
%   Supp Fig 16c,d and 16e,f -- Same as Supp Fig 14a,b and 14c,d,
%         respectively, EXCEPT:
%         bitUsePC = true;
%   Supp Fig 17j,m -- Same as Supp Fig 16a,b EXCEPT
%         bitUsePC_CC = true;
%         maxVar2Plot = 'none';
%   Supp Fig 18a,b (left panel) -- Same as Figure 4c,d EXCEPT
%         var2Plot =  {'V'};
%         bitSumDim = true; 
%      and (middle panel): 
%         bitUsePC = true; 
%      and (right panel):
%         bitUsePC_CC = true;
%   Supp Fig 18c,d -- Same as Supp Fig 18c,d EXCEPT
%         var2Plot =  {'RSV'};
%         bitMatchSignal2SRA = false;
%   Supp Fig 26a-d -- Same as Supp Fig 14a-d EXCEPT
%         eN = 3;
%         bitTimeOfEvent = true;


plotParamSet = {};

% set which temporal epoch to use by passing integer value:
%   1: present-trial epoch projected onto present-trial sRAs
%   2: previous-trial epoch projected onto present-trial sRAs
%   3: previous-trial epoch projected onto previous-trial sRAs
%   4: present-trial epoch coded WRT PRESENT-trial conditions projected
%       onto previous-trial sRAs.
%   5: present-trial epoch coded WRT PREVIOUS-trial conditions projected
%       onto previous-trial sRAs. 
eN = 1;

% set variance (V), relevant signal variance (RSV), irrelevant
% signal variance (ISV). Set to 'R2' to plot the R2_tk field from
% TDRSummary.runTDRSummary, which measures the contribution of the
% signal-specific term to the total estimated firing rate. Set to 'normdRA'
% to plot the norm (magnitude) of the non-unit vector dRA(t). This is an
% estimate of the magnitude of signal at time t. Include multiple entries
% to plot multiple signals (except for R2 and normdRA, which can only take
% 1 signal)
var2Plot =  {'RSV','ISV'}; % {'V'}; {'R2'}; {'RSV','ISV','ARSV'}; {'normdRA'}
 
% set type of variance to plot by passing integer value:
%   1: absolute variance 
%   2: variance explained expressed as percentage relative to total variance
%   3: (DEFUNCT)
%   4: p-value of variance explained based on comparing veridical value to
%       gamma distribution fit to empirical null distribution
varType = 4;

% logical on whether to sum across dimensions/signals and plot a single
% line
bitSumDim = false; 
 
% logical on whether to plot RSV/ISV only for the signal matched to the sRA
% (true) (e.g., benefit RSV on the benefit sRA), or plot RSV for all
% signals on each sRA (false). When the latter, bitSepAx is forced to be
% true. Note that when plotting PCs, bitMatchSignal2SRA=false is implied,
% as the RSV/ISV of each signal is plotted for each PC. When summing across
% dimensions, each signal is summed across dimensions (sRAs or PCs) and
% plotted separately.
bitMatchSignal2SRA = true;

% logical on whether to plot V/RSV/ISV for the PCs instead of the sRAs
bitUsePC = false;

% logical on whether to use PCs from common condition response. When TRUE,
% forces bitUsePC = TRUE. Default FALSE.
bitUsePC_CC = false;

% logical on whether to plot each signal on a separate axis. When plotting
% RSV/ISV for PCs, this is the default and overrides setting here. Likewise
% when plotting variance of PCs, there should just be one axis.
bitSepAx = false;

% set whether to plot horizontal bars representing times from which sRA was
% built
bitTimeOfSRA = true;

% set whether to plot square marker with horizontal error bars representing
% median and 95% timing of key task events, respectively.
bitTimeOfEvent = false;

% logical on whether to plot RSV/ISV for off-target signals based on
% PARTIAL CORRELATIONS 
bitRSVPartial = true;

% set type of max variance to plot
%       'dPC' -- maximum variance possible at each time step (i.e., dynamic
%          first PC computed independently at each moment in time). This is
%          forced when variance or variance explained is plotted (var2Plot
%          = 'V')
%       'dRA' -- variance explained by dynamic RA computed from TDR. This
%          is only available when plotting RSV/ISV (var2Plot = 'RSV' or
%          'ISV')
%       'serialSRA' -- variance explained by sRA computed via oTDR
%          independently at each time bin. Can be used for any type of
%          variance.
%       'none' -- no max variance plotted (DEFAULT).
maxVar2Plot = 'none';

% logical on whether to plot bounds for RSV based on VE (upper bound for
% on-target signal) and based on correlation between on- and off-target
% predictors (upper bound of non-leak for off-target signals). Can only
% plot RSV bounds when plotting max variance for the on-target signal (see
% maxVar2Plot)
bitRSVBound = false;

% logical on whether to plot difference between off-target RSV and RSV
% upper-bound (due to correlation between signals). Only applies when
% plotting RSV with bitMatchSignal2SRA == false. Forces bitRSVBound =
% false.
bitRSVBoundDiff = false;

% set analysis (tdr vs otdr). To work for TDR, one has to add code to
% select which TDR vector to plot.
fn = otdr;
% fn = 'oTDRSummary';

% set which total variance to normalize by: totalVar_t or totalVarMedD_t
totVarFN = 'totalVar_t';

% logical on whether to use horizontal line (instead of fill patch) for
% p-value treshold. When true, overrides bit4EPS.
bitLine4PThresh = true;

% logical on whether to change patches to transparent given matlab's
% strange behavior when exporting filled patches. Irrelevant when
% bitLine4PThresh == TRUE.
bit4EPS = false; 

%%%%%%%%%%% PLOTTING
% if plotParamSet is not empty, plot all combinations of param values in
% plotParamSet
if ~isempty(plotParamSet)
    
    % intantiate
    paramSetPos = ones(size(plotParamSet,1),1);
    
    % loop through all param values, determining max number in each row
    clear paramSetPosAll contingencySet
    for i = 1:size(plotParamSet,1)
        paramSetPosAll{i} = 1:length(plotParamSet{i,2});
    end
    
    % The following code to make contingencySet was suggested by John
    % Cunningham, 23 Aug 2007:
    nParamToSplit = length(paramSetPosAll);
    nRow = 1;
    for i = 1:nParamToSplit
        nRow = nRow * length(paramSetPosAll{i});
    end
    sublength(nParamToSplit) = 1; % last column changes value at every row
    for i = nParamToSplit-1:-1:1
        sublength(i) = sublength(i+1) * length(paramSetPosAll{i+1});
    end
    for i = 1:nParamToSplit
        nCopies = nRow / sublength(i) / length(paramSetPosAll{i});
        foo = repmat(paramSetPosAll{i},sublength(i),1);
        foo = foo(:);
            contingencySet(:,i) = repmat(foo,nCopies,1);
    end
    clear foo
    %%%%%%% end Cunningham
    
    % make loop index
    indSet = 1:size(contingencySet);
else
    indSet = 1;
end

% loop through rows of contingency set
for i = indSet
    
    if ~isempty(plotParamSet)
        % loop through variables
        for j = 1:size(contingencySet,2)
            
            % set parameter depending on type
            if ischar(plotParamSet{j,2}{contingencySet(i,j)})
                eval(sprintf('%s = %s;',plotParamSet{j,1},plotParamSet{j,2}{contingencySet(i,j)}));
                %                 sprintf('%s = %s',plotParamSet{j,1},plotParamSet{j,2}{contingencySet(i,j)})
            else
                eval(sprintf('%s = %d;',plotParamSet{j,1},plotParamSet{j,2}{contingencySet(i,j)}));
                %                 sprintf('%s = %d',plotParamSet{j,1},plotParamSet{j,2}{contingencySet(i,j)})
            end
            
        end
    end
    
    % skip certain parameter combinations that are intentionally
    % illegal. Only apply when looping automatically through parameters
    if ~isempty(plotParamSet) && ...
            (all(ismember(var2Plot,{'V'})) && ~bitMatchSignal2SRA) || ...
            (bitSumDim && bitSepAx) || ...
            (bitUsePC && bitMatchSignal2SRA && any(ismember(var2Plot,{'RSV','ISV'}))) || ... % this setting is tolerated, but redundant with ~bitMatchSignal2SRA
            (bitSumDim && bitMatchSignal2SRA && ~bitUsePC && any(ismember(var2Plot,{'RSV','ISV'}))) % this setting is tolerated, but redundant with ~bitMatchSignal2SRA
%             (bitSumDim && ~bitMatchSignal2SRA && any(ismember(var2Plot,{'RSV','ISV'}))) || ...
        warning('Parameter combination not allowed: var=%s, varType=%d, bitSumDim=%d, bitMatchSignal=%d, bitUsePC=%d, bitSepAx=%d',...
        cell2mat(var2Plot),varType,bitSumDim,bitMatchSignal2SRA,bitUsePC,bitSepAx);
        continue
    end
    
    % plot
    figH = oTDR_plotSRAVariance(Summary,...
        'eN',eN,...
        'var2Plot',var2Plot,...
        'varType',varType, ...
        'bitSumDim',bitSumDim,...
        'bitMatchSignal2SRA',bitMatchSignal2SRA,...
        'bitUsePC',bitUsePC,...
        'bitUsePC_CC',bitUsePC_CC,...
        'bitSepAx',bitSepAx,...
        'bitTimeOfEvent',bitTimeOfEvent,...
        'bitTimeOfSRA',bitTimeOfSRA,...
        'bitRSVPartial',bitRSVPartial,...
        'maxVar2Plot',maxVar2Plot,...
        'bitRSVBound',bitRSVBound,...
        'bitRSVBoundDiff',bitRSVBoundDiff,...
        'fn',fn,...
        'totVarFN',totVarFN,...foo(bitTarget(i,:))
        'bgColor',bgColor,...
        'lineColor',lineColor,...
        'textColor',textColor,...
        'condLabel',condLabel,...
        'alpha',alpha,...
        'bitLine4PThresh',bitLine4PThresh,...
        'bit4EPS',bit4EPS);
    
    % name figure with parameter set number
    figH.Name = sprintf('Fig %d - var=%s, varType=%d, bitSumDim=%d, bitMatchSignal=%d, bitUsePC=%d, bitSepAx=%d',...
        figH.Number,cell2mat(var2Plot),varType,bitSumDim,bitMatchSignal2SRA,bitUsePC,bitSepAx);
    
    if ~isempty(plotParamSet)
        % save figure as EPS and .mat
        printFigure(figH,fullfile('/Users/Daniel/FIG',animalName,'TDR/variance_explained_by_sRA'),1,1,2);
        savefig(figH,fullfile('/Users/Daniel/FIG',animalName,'TDR/variance_explained_by_sRA',get(figH,'Name')));
        
        % close figure
        close(figH);
    end
end


%% Project data onto Regression axes (sRA, PC, CC PCs)

% define type of axes to use:
% 'sRA' -- for sRAs
% 'PCs' -- for data PCs
% 'PCs_CC' -- for common condition PCs
% 'dRA' -- for dRAs -- plots projection of data(t) onto dRA(t)
RA = 'sRA';

% logical on whether to use sRAs for present offer (FALSE) or previous
% offer (TRUE). NOTE: only supports RA = {sRA, PCs, PCs_CC} and
% bitSerialOTDR=false.
bitSRA_PP = false;

% logical to use the projection onto the serial sRA computed independently
% at each regression time
bitSerialOTDR = false; 

% define which data to project
% 'proj' -- for mean subtracted data; 
% 'projNoMeanSubtraction' -- for non-mean subtracted
% 'projLateReject' -- separate conditions for early vs. late rejections
% 'projRejectAlign' -- aligned to time of rejection (reject chocies only),
%   no mean subtraction
% 'projRejectAlign_wLateReject' -- aligned to time of rejection (reject
%   chocies only), no mean subtraction, AND separate conditions for early
%   vs. late rejections 
% 'projPP' -- projects previous offer responses onto present offer sRAs
%       (when bitSRA_PP == false), coded WRT previous offer conditions.
% 'projPresent' -- projects present trial responses onto previous offer
%       sRAS (when bitSRA_PP == true), coded WRT present offer conditions.
% 'projPresRespPrevCond' -- projects present offer responses onto previous
%       offer sRAs (when bitSRA_PP == true), coded WRT previous offer
%       conditions.
d2p = 'proj'; 

% logical to remove axes and plot scale bar
bitRemoveAxes = true;

% logical on whether to plot each offer on separate axes
bitSepAxPerOffer = false;

% logical on whether to include a legend
bitLegend = true;

%%%%%% MOVIE PARAMS
bitMovie = false; % make movie of projections evolving over time

movie_save_path = ''; % absolute path to directory for saving movie files

% movie params
mov.dimN = [1 2 3]; % dimensions to include
mov.offer = []; % which offers to include. leave empty for all
mov.choice = [0]; % which choices to include. leave empty for all
mov.offerPre = [0:8]; % which offers to include prior to plotting movie. leave empty for none
mov.choicePre = [1]; % which choices to include prior to plotting movie. leave empty for none
mov.framerate = 5; % frame rate -- leave empty for default
mov.kernal_width = 5; % width of gaussian kernal for smoothing
mov.nTrail = [10]; % number of trailing time points to preserve. Leave empty or Inf to plot all.
% transparency of trailing time points earlier than those in nTrail. Set to
% 1 to have inf tail. Set to 0 to have tail as long as nTrail. Also sets
% transparency of any trajectories plotted prior to plotting movie.
mov.trailTrans = 0.4; 

%%%%%%%%%%%%%%%%
% checks
if ~bitSRA_PP && (strcmp('projPresent',d2p) || strcmp('projPresRespPrevCond',d2p))
    error('When setting d2p to "projPresent" or "projPresRespPrevCond", you must use the previous-trial sRAs (i.e., bitSRA_PP = true). Indeed these projections are defined as projecting present-trial responses onto the previous-trial sRAs.')
end

if bitSRA_PP && strcmp('projPP',d2p)
    error('When setting d2p to "projPP", you must use the present-trial sRAs (i.e., bitSRA_PP = false). Indeed these projections are defined as projecting previous-trial responses onto the present-trial sRAs.')
end

if bitSRA_PP && any(strcmp(d2p,{'projLateReject','projRejectAlign','projRejectAlign_wLateReject'}))
    error('The projection d2p = "%s" is not defined for the previous-trial sRAs (i.e., bitSRA_PP = true). Set bitSRA_PP = False.',d2p);
end

if bitMovie && isempty(movie_save_path)
    movie_save_path = input('Enter absolute path at which to save movie files: ','s');
end

% mean-center each column (I don't think this is necessary)
% data = bsxfun(@minus,Summary.(otdr).data,mean(Summary.(otdr).data));
% data = Summary.(otdr).data;

% set summary field name to use
if bitSRA_PP
    % previous offer sRAs
    sfn = [otdr,'_PP'];
    % add prefix for sRA labels
    condLabelTemp = cellfun(@(c) ['Prev ',c],condLabel,'UniformOutput',0);
    analName = 'Previous trial';
else
    % present offer sRAs (or dRAs)
    condLabelTemp = condLabel;

    if strcmp(RA,'dRA')
        sfn = tdr;
        % turn off bitSerialOTDR since it does not apply to dRAs
        bitSerialOTDR = false;
    elseif bitSerialOTDR
        sfn = otdrSerial;
    else
        sfn = otdr;
    end
    analName = 'Present trial';
end

% determine set of conditions
if ismember(d2p,{'projLateReject','projRejectAlign_wLateReject'}) 
    condSet = Summary.(otdr).Predictors_withLateReject;    
else
    condSet = Summary.(otdr).Predictors;
end
offerSetMaster = unique(condSet(:,1)); % necessary in case we eliminate conditions below

% special case -- when plotting rejection aligned data, we only have reject
% conditions, so we have to pare down condSet
if ismember(d2p,{'projRejectAlign','projRejectAlign_wLateReject'})
    condSet(condSet(:,2)~=0,:) = [];
end
% nTime = length(Summary.(otdr).Times.all_times);
% nCond = size(condSet,1);
offerSet = unique(condSet(:,1));
nOffer = numel(offerSet);
% clear p
% p{1} = reshape(data * Summary.(sfn).sRA.sRA1,nTime,nCond)';
% p{2} = reshape(data * Summary.(sfn).sRA.sRA2,nTime,nCond)';
% p{3} = reshape(data * Summary.(sfn).sRA.sRA3,nTime,nCond)';

% define colorset
if strcmp(bgColor,'w')
    colorSet = CB_offerColor();
else
    colorSet = jet(numel(offerSetMaster));
end
% eliminate color entries that are not present in offerSet
[~,foo] = setdiff(offerSetMaster,offerSet);
colorSet(foo,:) = [];

% find number of regression axes -- always base this on the sRAs, so that
% we never plot more PCs than sRAs
if strcmp('dRA',RA) 
    % note that TDR includes the constant as an additional axis
    nRA = size(Summary.(sfn).dRAs,3);
elseif bitSerialOTDR
    nRA = size(Summary.(sfn).sRA.RA,3);
else
    nRA = size(Summary.(sfn).sRA.RA,2);
end

% If making movie, open file
if bitMovie
    % make movie path
    if bitSerialOTDR
        foo = '_serial';
    else
        foo = '';
    end
    goo = '';
    if ~isempty(mov.offer)
        goo = [goo,'_offer(',num2str(mov.offer),')'];
    end
    if ~isempty(mov.choice)
        goo = [goo,'_choice(',num2str(mov.choice),')'];
    end
    mov.path = fullfile(movie_save_path,...
        sprintf('%s_%s%s_%s_%s%s',animalName(1),RA,foo,d2p,datestr(round(now),'yyyy_mm_dd'),goo));
    clear foo goo
    v = VideoWriter(mov.path,'MPEG-4');
    v.Quality = 100;    
    if ~isempty(mov.framerate)
        v.FrameRate = mov.framerate;
    end
    open(v);
    
end

%%%%%

figure;
whitebg(gcf,bgColor)
figName = sprintf('Projection onto %s for %s',RA,analName);
set(gcf,'Color',bgColor,'Name',figName);

axN = 0;
axH = [];
if bitSepAxPerOffer
    % separate axes for each offer
    nRow = nRA;
    nCol = nOffer;
else
    nRow = 1;
    nCol = nRA;
end

% % loop through signals of interest
% for a = 1:length(p)
%     axN = axN+1;
%     subplot(nRow,nCol,axN)
%     hold on;
%     for i = 1:size(condSet,1)
%         c = colorSet(ismember(offerSet,condSet(i,2)),:);
%         switch condSet(i,1)
%             case 0
%                 lw = 1;
%             case 1
%                 lw = 2;
%         end
%         plot(Summary.(sfn).Times.all_times,...
%             p{a}(i,:),'Color',c,'LineWidth',lw);
%     end
%     ylabel(condLabelTemp{a},'FontSize',18);
%     set(gca,'box','off','FontSize',16);
%     xlim(minmax(Summary.(sfn).Times.all_times));    
%     % line at x = 0
%     yLim = get(gca,'YLim');
%     line([0 0],yLim,'Color',lineColor);
% end


if bitMovie
    %%% MAKE MOVIE

    % prepare data
    % get time vector (not displayed)
    if ismember(d2p,{'projRejectAlign','projRejectAlign_wLateReject'})
        % different times for aligning to rejection
        tSet = Summary.(sfn).Times.all_times_rejectAlign;
    else
        tSet = Summary.(sfn).Times.all_times;
    end
    
    % get projections
    if strcmp('dRA',RA)
        % from dRA
        y = Summary.(sfn).(d2p);
    else
        y = Summary.(sfn).(RA).(d2p);
    end
    
    % pare down dimensions as requested
    if ~isempty(mov.dimN)
        y = y(:,:,mov.dimN);
    end
    % determine number of times and dims
    [T,~,D] = size(y);
    
    % correct number of trailing points -- when empty, set to inf so as to
    % include all
    if isempty(mov.nTrail)
        mov.nTrail = Inf;
    end
    
    %%% smooth projection
%     dt = tSet(2) - tSet(1); % time bin width
%     kern = gausswin(mov.nKernPoint,1/mov.kernSigma); % construct kernal
%     kern = kern / (sum(kern)*dt*mov.nKernPoint); % normalize kernal

    % filter along time dimension. permute necessary because smoothts
    % operates along the row (not columns) and time is in the columns.
    y = permute(smoothdata(permute(y(:,:,:),[2 1 3]),'gaussian',mov.kernal_width),[2 1 3]);

    % default values
    mov.cameraPos = []; 
    mov.cameraViewAngle = [];
    mov.cameraTarget = [];
%     % determine camera angle
%     switch RA
%         case 'dRA'
%             mov.cameraPos = [105.9498   54.8861   14.0527]; % default
%         case 'sRA'
%             mov.cameraPos = [105.9498   54.8861   14.0527]; % default
%         case 'PCs'
%             mov.cameraPos = [105.9498   54.8861   14.0527]; % default
%         case 'PCs_CC'
%             mov.cameraPos = [105.9498   54.8861   14.0527]; % default
%     end
%     mov.cameraPos = [61.7891   48.6960   18.6284]; 
%     mov.cameraPos = [117.8535   56.5459   15.4121];
%     mov.cameraPos = [];
    %     mov.cameraViewAngle = 8.3581;
%     mov.cameraTarget = [5.0000   -0.5000    2.5000];
    mov.viewStart = [180 3];
    mov.viewEnd = [92 3];
    % build view vector that remains on start position until end of benefit
    % period
    % find end of benefit position
    posBenEnd = find(tSet <= 0.5,1,'last');
    viewVect = bsxfun(@times,ones(T,2),mov.viewStart);
    viewVect(posBenEnd:end,:) = [linspace(mov.viewStart(1),mov.viewEnd(1),T-posBenEnd+1)'...
        linspace(mov.viewStart(2),mov.viewEnd(2),T-posBenEnd+1)'];
    
    % make axes
    subplot(1,1,1)
    axis tight manual 
    
    % figure size is important
    set(gcf,'WindowStyle','normal','Position',[1 1 800 600]);
    
    % axes appearance
    set(gca,'FontSize',18,'box','off');
%     set(gca,'XTick',[],'YTick',[],'ZTick',[]);
    set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
    grid on
    
%     % camera position, etc
    if ~isempty(mov.cameraPos)
        set(gca,'CameraPosition',mov.cameraPos);
    end
    if ~isempty(mov.cameraViewAngle)
        set(gca,'CameraViewAngle',mov.cameraViewAngle);
    end
    if ~isempty(mov.cameraTarget)
        set(gca,'CameraTarget',mov.cameraTarget);
    end
    
    
    % find and set limits, build list of axes labels and apply them.
    % This is important to do before paring down conditions so that
    % constant limits can be maintained regardless of specific set of
    % conditions used.
    projMinMax = NaN(D,2);
    axLabel = cell(D,1);
    for d = 1:D
        projMinMax(d,1) = min(min(y(:,:,d)));
        projMinMax(d,2) = max(max(y(:,:,d)));
        
        % axis label depending on type of regression axes
        switch RA
            case 'dRA'
                if mov.dimN(d) <= length(condLabelTemp)
                    axLabel{d} = condLabelTemp{mov.dimN(d)};
                else
                    axLabel{d} = 'Constant';
                end
            case 'sRA'
                axLabel{d} = condLabelTemp{Summary.(sfn).(RA).param4RA(mov.dimN(d))};
            case 'PCs'
                axLabel{d} = sprintf('PC_%d',mov.dimN(d));            
            case 'PCs_CC'
                axLabel{d} = sprintf('CC PC_%d',mov.dimN(d));            
        end
        if bitSerialOTDR
            axLabel{d} = horzcat('Serial ',axLabel{d});
        end
        
        switch d
            case 1
                xlim(ceilfix(projMinMax(d,:)));
                xlabel(axLabel{d});
            case 2
                ylim(ceilfix(projMinMax(d,:)));
                ylabel(axLabel{d});
            case 3
                zlim(ceilfix(projMinMax(d,:)));
                zlabel(axLabel{d});
        end
    end     
    
    % determine pre-plot conditions (will appear prior to movie). must
    % supple both offer and condition
    bitPre = false(size(condSet,1),1);
    if ~isempty(mov.offerPre) && ~isempty(mov.choicePre)
        bitPre = ismember(condSet(:,1),mov.offerPre) & ismember(condSet(:,2),mov.choicePre);
    end
    yPre = y(:,bitPre,:);
    condSetPre = condSet(bitPre,:);

    % pare down conditions as requested
    bitKeep = true(size(condSet,1),1);
    if ~isempty(mov.offer)
        bitKeep = bitKeep & ismember(condSet(:,1),mov.offer);
    end
    if ~isempty(mov.choice)
        bitKeep = bitKeep & ismember(condSet(:,2),mov.choice);
    end
    y(:,~bitKeep,:) = [];
    condSet(~bitKeep,:) = [];
    
    % determine number of conditions
    [~,C,~] = size(y);
    
    % loop through time points
    for t = 1:T
        
        % force new lines on axes
        set(gca,'nextplot','replacechildren');
        
        % set view angle
%         camorbit(10,0,'camera');
        view(viewVect(t,:));
        
        % add starting point, time, and trial period
        if tSet(t) < 0
            period = 'HOLD';
        elseif tSet(t) >= 0 && tSet(t) <= 0.5
            period = 'OFFER';
        elseif tSet(t) < tSet(end)
            period = 'WORK';
        elseif tSet(t) >= tSet(end)
            period = 'REWARD';
        else
            period = '';
        end
        
        % determine oldest time position to plot given trailing time points
        tStart = max(t - mov.nTrail+1,1);
        
        if D == 2
            plot(mean(y(1,:,1)),mean(y(1,:,2)),'o','g','MarkerSize',12,'MarkerFaceColor','g');
            text(1,1,sprintf('t = %0.2f s',tSet(t)),'FontSize',18,'Units','normalized','Color',textColor);
            text(0,1,period,'FontSize',24,'Units','normalized','Color',textColor);
        elseif D == 3
            plot3(mean(y(1,:,1)),mean(y(1,:,2)),mean(y(1,:,3)),'og','MarkerSize',12,'MarkerFaceColor','g');
            text(1,1,1,sprintf('t = %0.2f s',tSet(t)),'FontSize',18,'Units','normalized','Color',textColor);
            text(0,1,1,period,'FontSize',24,'Units','normalized','Color',textColor);
        end
        
        % hold lines until all conditions are plotted
        set(gca,'nextplot','add');
        
        % PRE-CONDITIONS -- plotted in entirety
        if ~isempty(yPre)
            for i = 1:size(condSetPre,1)
                % find color based on offer
                c = colorSet(ismember(offerSet,condSetPre(i,1)),:);
                switch condSetPre(i,2)
                    case 0
                        % rejection
                        lw = 1;
                    case 1
                        % accept
                        lw = 3;
                end
                
                if D == 2
                    plot(yPre(:,i,1),yPre(:,i,2),'Color',[c,mov.trailTrans],'LineWidth',lw);
%                     % add SQAURE to mark end point
%                     plot(yPre(end,i,1),yPre(end,i,2),'s','Color',c,'MarkerFaceColor',c,'MarkerSize',16);
                elseif D == 3
                    h = plot3(yPre(:,i,1),yPre(:,i,2),yPre(:,i,3),'Color',[c,mov.trailTrans],'LineWidth',lw);
%                     % add SQAURE to mark end point
%                     plot3(yPre(end,i,1),yPre(end,i,2),yPre(end,i,3),'s','Color',c,'MarkerFaceColor',c,'MarkerSize',16);
                else
                    error('number of dimensions (%d) cannot be plotted',D)
                end
            end
        end
        
        % loop through conditions
        for i = 1:C
            
            % find color based on offer
            c = colorSet(ismember(offerSet,condSet(i,1)),:);
            switch condSet(i,2)
                case 0
                    % rejection
                    lw = 1;
                case 1
                    % accept
                    lw = 3;
            end
            
            % MOVIE CONDITIONS
            % plot all conditions from beginning to current time step
            if D == 2
                plot(y(tStart:t,i,1),y(tStart:t,i,2),'Color',c,'LineWidth',lw);
                % add transparency of longer trail 
                if tStart > 1
                    plot(y(1:tStart,i,1),y(1:tStart,i,2),'Color',[c,mov.trailTrans],'LineWidth',lw);
                end
                % add SQAURE to mark end point
                plot(y(t,i,1),y(t,i,2),'s','Color',c,'MarkerFaceColor',c,'MarkerSize',16);
            elseif D == 3
                plot3(y(tStart:t,i,1),y(tStart:t,i,2),y(tStart:t,i,3),'Color',c,'LineWidth',lw);
                % add transparency of longer trail 
                if tStart > 1
                    plot3(y(1:tStart,i,1),y(1:tStart,i,2),y(1:tStart,i,3),'Color',[c,mov.trailTrans],'LineWidth',lw);
                end
                % add SQAURE to mark end point
                plot3(y(t,i,1),y(t,i,2),y(t,i,3),'s','Color',c,'MarkerFaceColor',c,'MarkerSize',16);
            else
                error('number of dimensions (%d) cannot be plotted',D)
            end
            
        end
        
        % capture frame and write to video
        frame = getframe(gcf);
        writeVideo(v,frame);
        
    end
    
    % close video
    close(v);    
    
else
    % plot projection on separate axes
    % loop through signals of interest
    for a = 1:nRA
        
        % if plotting all offers on same axes:
        if ~bitSepAxPerOffer
            axN = axN+1;
            axH(axN) = subplot(nRow,nCol,axN);
            hold on;
        end
        
        % instantiate separate for each panel
        dataLim = [Inf -Inf];
        hL = [];
        legText = {};
        
        % determine title depending on type of regression axes
        switch RA
            case 'dRA'
                if a <= length(condLabelTemp)
                    titleStr = condLabelTemp{a};
                else
                    titleStr = 'Constant';
                end
            case 'sRA'
                titleStr = condLabelTemp{Summary.(sfn).(RA).param4RA(a)};
            case 'PCs'
                titleStr = sprintf('PC_%d',a);            
            case 'PCs_CC'
                titleStr = sprintf('CC PC_%d',a);            
        end
        if bitSerialOTDR
            titleStr = horzcat('Serial ',titleStr);
        end        
        
        % loop through conditions
        for i = 1:size(condSet,1)
            c = colorSet(ismember(offerSet,condSet(i,1)),:);
            switch condSet(i,2)
                case 0
                    % rejection
                    lw = 1;
                    choiceTxt = 'reject';
                case 1
                    % accept
                    lw = 2;
                    choiceTxt = 'accept';
            end
            % if including early vs. late rejections
            if size(condSet,2) > 3
                switch condSet(i,4)
                    case 0
                        % non-late rejection (i.e. early rejection or accept)
                        lineSty = '-';
                        rejectTxt = '';
                    case 1
                        % late rejection
                        lineSty = '--';
                        rejectTxt = 'late ';
                end
            else
                rejectTxt = '';
                lineSty = '-';
            end                        
            if ismember(d2p,{'projRejectAlign','projRejectAlign_wLateReject'})
                % di.fferent times for aligning to rejection
                x = Summary.(sfn).Times.all_times_rejectAlign;
            elseif strcmp('projPP',d2p)
                x = Summary.(sfn).Times.all_times_PP;
            elseif strcmp('projPresent',d2p) || strcmp('projPresRespPrevCond',d2p) 
                x = Summary.(sfn).Times.all_times_present;
            else                
                x = Summary.(sfn).Times.all_times;
            end
            if strcmp('dRA',RA)
                % from dRA
                y = Summary.(sfn).(d2p)(:,i,a);
            else
                y = Summary.(sfn).(RA).(d2p)(:,i,a);
            end
                        
            % if separate axes for each offer, make subplot here (or change
            % focus to former subplot)
            yStd = [];
            if bitSepAxPerOffer
                axN = (a-1)*nOffer + find(offerSet==condSet(i,1));
                axH(axN) = subplot(nRow,nCol,axN);
                hold on;
                
                % title includes title for regression axes and offer
                title(sprintf('%s, %d reward',titleStr,condSet(i,1)));
                
                % special case: if plotting Late Rejections, also include
                % standard projections in background
                if ismember(d2p,{'projLateReject','projRejectAlign_wLateReject'}) 
                    
                    % get name of std projection depending on aligning to
                    % offer or rejection
                    if strcmp(d2p,'projLateReject')
                        foo = 'proj';
                    else
                        foo = 'projRejectAlign';
                    end
                        
                    % find matching condition (ignoring late vs early
                    % reject)
                    condSetStd = Summary.(sfn).Predictors; % std condition set                    
                    posCondStd = ismember(condSetStd,condSet(i,1:size(condSetStd,2)),'rows');
                    yStd = Summary.(sfn).(RA).(foo)(:,posCondStd,a);
                    plot(x,yStd,'Color',[0.7 0.7 0.7],'LineWidth',lw);
                end
            end
            
            % plot
            hL(end+1) = plot(x,y,'Color',c,'LineWidth',lw,'LineStyle',lineSty);
            
            % legend text
            legText{end+1} = sprintf('%d reward, %s%s',condSet(i,1),rejectTxt,choiceTxt);            
            
            % update limits
            dataLim(1) = min(dataLim(1),min([y;yStd]));
            dataLim(2) = max(dataLim(2),max([y;;yStd]));  
            clear y            

        end
        
        % add singleton (always for accepts) if available
        if ~strcmp('dRA',RA) && isfield(Summary.(sfn).(RA),'projExcluded') && ...
                ~all(isnan(Summary.(sfn).(RA).projExcluded(:,:,a))) && ...
                ~ismember(d2p,{'projRejectAlign','projRejectAlign_wLateReject'})
            c = colorSet(ismember(offerSet,8),:);
            lw = 2;
            y =squeeze(Summary.(sfn).(RA).projExcluded(:,:,a));
            hL(end+1) = plot(Summary.(sfn).Times.all_times,y,...
                'Color',c,'LineWidth',lw,'LineStyle','--');
            dataLim(1) = min(dataLim(1),min(y));
            dataLim(2) = max(dataLim(2),max(y));        
            legText{end+1} = '8 reward singleton, accept';
        end

        % draw title if not plotting separate axes
        if ~bitSepAxPerOffer
            title(titleStr,'FontSize',18);
        end
        
        % determine which on axes to label y axis
        if bitSepAxPerOffer
            axHPos_Col1 = [1:nCol:length(axH)];
            axHPos_Col1 = axHPos_Col1(end); % current row only
            
            % save as first axes for RA
            axHPos_4RA = axHPos_Col1;
        else
            axHPos_Col1 = 1;
            
            axHPos_4RA = a;
        end
        ylabel(axH(axHPos_Col1),'Projection magnitude (a.u.)','FontSize',16);
        
        % determine which on axes to label x axis
        foo = [];
        if bitSepAxPerOffer 
            if a==nRow
                foo = (a-1)*nOffer + ceil(nCol/2);
            end
        else
            foo = ceil(nCol/2);
        end        
        if ~isempty(foo) && foo <= length(axH)
            if ismember(d2p,{'projRejectAlign','projRejectAlign_wLateReject'})
                xlabel(axH(foo),'Time from rejection (s)','FontSize',16);
            else
                xlabel(axH(foo),'Time from offer (s)','FontSize',16);
            end
        end
        clear foo
        
        % get axes handle position for all axes for current RA
        if bitSepAxPerOffer
            axHPos = (a-1)*nOffer + [1:nCol];
        else
            axHPos = a;
        end
        set(axH(axHPos),'box','off','FontSize',16);
                
        % ylim, round to nearest Integer I
        I = 1;
        yLim(1) = dataLim(1) - mod(dataLim(1),I);
        yLim(2) = dataLim(2) + (I-mod(dataLim(2),I));
                
        % apply limits and add line to each axes for current RA
        for j = 1:length(axHPos)
            xlim(axH(axHPos(j)),minmax(x));    
            ylim(axH(axHPos(j)),yLim);
            
            % line at x = 0
            plot(axH(axHPos(j)),[0 0],yLim,'-','Color',lineColor);
        end
        
        % label line as offer or rejection (only in first column)
        if ismember(d2p,{'projRejectAlign','projRejectAlign_wLateReject'})
            foo = 'Rejection';
        elseif strcmpi(d2p,'projPP') || (strcmp(d2p,'proj') && bitSRA_PP)
            foo = 'Fixation';
        else
            foo = 'Offer';
        end
        axes(axH(axHPos_4RA));
        text(0,yLim(2),foo,'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','bottom');

        % legend (only if all offers on same axes)
        if ~bitSepAxPerOffer && a == 1 && bitLegend
            legend(hL,legText,'box','off','Location','NorthEast');
        end
        
        %%%% SCALE BARS
        if bitRemoveAxes 
            % remove axes and replace with scale bars
            for j = 1:length(axH)
                % do separately for x and y axis to preserve title
                set(get(axH(j),'XAxis'),'Visible','off');
                set(get(axH(j),'YAxis'),'Visible','off');
            end
            
            % focus on first axes for RA
            axes(axH(axHPos_4RA));

            % horizontal scalebar
            x = [0 1]-0.5;
            y = [1 1]*yLim(1)-0.04*range(yLim);
            h = line(x,y,'Color',textColor);
            hT = text(x(1) + range(x)/2,y(1),sprintf('%0.0f s',range(x)),'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','top');
            set([h hT],'Clipping','off');
            % vertical scalebar
            y = y(1) + [0 2];
            x = x(1) * [1 1];
            h = line(x,y,'Color',textColor);
            hT = text(x(1),y(1) + range(y)/2,sprintf('%0.0f a.u. ',range(y)),'FontSize',16,'HorizontalAlignment','right','VerticalAlignment','middle');
            set([h hT],'Clipping','off');
        end
    end
    
    %%% annotate
    annotateFig(figName,mfilename);
    
    % link a-axes
    linkaxes(axH,'x');
end


clear condLabelTemp

%% Plot anatomical/spatial location of signals

% this requires the PHYS struct from CB_plexon_allFile.m

% logical on whether to plot individual unit responses in grid format or
% plot mean response of site:
bitSiteMean = false;

% logical on whether to plot a scatter plot of the mean sRA component per
% site for signal 1 vs signal 2.
bitPlotSiteMeanScatter = false;

bgColor = 'w';
textColor = 'k';

%%%%%
% extract bit mask of cells to include
bitCell = Summary.(otdr).bitCellIncluded;

% % make sure that bitCell2Include vector matches the expected size of the
% % phys struct, since we are registering between DK's original code and
% % DK/GE's oTDR package
% if length(phys) ~= length(bitCell)
%     error('PHYS struct and Summary.(otdr).bitCellIncluded must match')
% end

% loc = [double([phys(bitCell).gridRow]')-64 [phys(bitCell).gridCol]'];
loc = [double(Summary.(otdr).siteInfo.row(bitCell))-64 Summary.(otdr).siteInfo.col(bitCell)];

sig = [Summary.(otdr).sRA.RA];
%%%%

% find row labels
foo = unique(loc(:,1));
% eliminate non letters
foo = foo(foo>=0);
% convert back to characters
foo = char(foo+64);
% convert to cell array of strings
foo = cellstr(foo);

[figH,siteCorr,siteCorrP] = spatialLocation_static_2D(sig,loc,'sigName',condLabel,...
    'bitSiteMean',bitSiteMean,'bitPlotSiteMeanScatter',bitPlotSiteMeanScatter,...
    'rowNameG1',foo,...
    'bgColor',bgColor,'textColor',textColor);
clear foo

disp('Correlation coefficient of site mean between signals:')
disp(siteCorr);
disp('P value associated with coefficient:');
disp(siteCorrP);

%% Common Condition Response projected onto RAs

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROJECTION OF CC onto CC PCs, PCs, and sRAs
figure;
whitebg(gcf,bgColor)
analName = 'Present trial';
figName = ['Common condition response - CC projection - ',analName];
set(gcf,'Color',bgColor,'Name',figName);

nRow = 1;
nCol = 3;

%%% projection onto CC PCs
subplot(nRow,nCol,1);
title('Common condition PCs')
hold on
h = plot(Summary.oTDRSummary.Times.all_times',squeeze(Summary.oTDRSummary.PCs_CC.projCC));
set(h,'Color',lineColor);
set(gca,'XLim',minmax(Summary.oTDRSummary.Times.all_times));
% line at x = 0
yLim = get(gca,'YLim');
line([0 0],yLim,'Color',lineColor);
% legend
legText = {};
for i = 1:length(h)
    set(h(i),'LineWidth',length(h) - (i-1));
    legText{i} = sprintf('Proj onto CC PC_%d (%0.1f%% temporal var exp)',i,...
        Summary.oTDRSummary.PCs_CC.varAnalysis.CC.VE_RA(i));
end
legend(legText,'box','off');
set(gca,'box','off','FontSize',16);
xlabel('Time from offer (s)','FontSize',16);
ylabel('Projection magnitude (a.u.)','FontSize',16);

%%% projection onto data PCs
subplot(nRow,nCol,2);
title('Data PCs')
hold on
h = plot(Summary.oTDRSummary.Times.all_times',squeeze(Summary.oTDRSummary.PCs.projCC));
set(h,'Color',lineColor);
set(gca,'XLim',minmax(Summary.oTDRSummary.Times.all_times));
% line at x = 0
yLim = get(gca,'YLim');
line([0 0],yLim,'Color',lineColor);
% legend
legText = {};
for i = 1:length(h)
    set(h(i),'LineWidth',length(h) - (i-1));
    legText{i} = sprintf('Proj onto data PC_%d (%0.1f%% temporal var exp)',i,...
        Summary.oTDRSummary.PCs.varAnalysis.CC.VE_RA(i));
end
legend(legText,'box','off');
set(gca,'box','off','FontSize',16);
xlabel('Time from offer (s)','FontSize',16);
ylabel('Projection magnitude (a.u.)','FontSize',16);

%%% projection onto sRA
subplot(nRow,nCol,3);
title('sRAs')
hold on
h = plot(Summary.oTDRSummary.Times.all_times',squeeze(Summary.oTDRSummary.sRA.projCC));
set(h,'LineWidth',2);
% line at x = 0
yLim = get(gca,'YLim');
line([0 0],yLim,'Color',lineColor);
set(gca,'XLim',minmax(Summary.oTDRSummary.Times.all_times));
legText = {};
lc = colorCategorical(length(condLabel));
for i = 1:length(h)
    set(h(i),'Color',lc(Summary.(otdr).sRA.param4RA(i),:));
    legText{i} = sprintf('Proj onto %s sRA (%0.1f%% temporal var exp)',condLabel{Summary.(otdr).sRA.param4RA(i)},...
        Summary.oTDRSummary.sRA.varAnalysis.CC.VE_RA(i));
end
legend(legText,'box','off');
set(gca,'box','off','FontSize',16);
xlabel('Time from offer (s)','FontSize',16);
ylabel('Projection magnitude (a.u.)','FontSize',16);

%%% add annotation
annotateFig(figName);

%%%%%%%%%%%%%%%%%%%%%%
% Bar chart of temporal variance of CC response explained by the three low-dimensional
% representations: sRAs, PCs, CC PCs

nRA = length(Summary.oTDRSummary.sRA.varAnalysis.CC.VE_RA);
nSpace = 3;

% holding var for temporal variance explained: space x RA
tempVE = NaN(nSpace, nRA);
fn = cell(nSpace,1);

for i = 1:nSpace
    switch i
        case 1
            fn{i} = 'PCs_CC';
        case 2
            fn{i} = 'sRA';
        case 3 
            fn{i} = 'PCs';
        otherwise
            error('no more than 3 spaces expected')
    end
    
    tempVE(i,:) = Summary.oTDRSummary.(fn{i}).varAnalysis.CC.VE_RA(1:nRA);
end

figure;
whitebg(gcf,bgColor)
figName = ['Temporal variance of CC response explained by 3D RAs - ',analName];
set(gcf,'Color',bgColor,'Name',figName);

% set colormap for sRA colors
colormap(colorCategorical(nRA,'CB'));
h = bar(tempVE);
set(gca,'XTickLabel',fn,'TickLabelInterpreter','none');
set(gca,'FontSize',16,'box','off');
ylabel('Temporal variance explained (%)','FontSize',16);


%%%%%%%%%%%%%%%%%%%%%%%%
% Bar chart of temporal variance per condition explained by the three 3D
% spaces: sRAs, PCs, and PCc_CC

figure;
whitebg(gcf,bgColor)
figName = ['Temporal variance of per condition response explained by 3D RAs - ',analName];
set(gcf,'Color',bgColor,'Name',figName);

data = [];
spaceDesc = {};

% loop through each set of RAs
for i=1:3
    switch i
        case 1
            spaceFn = 'sRA';
            spaceDesc{end+1} = 'sRAs';
        case 2
            spaceFn = 'PCs';
            spaceDesc{end+1} = 'PCs';
        case 3
            spaceFn  = 'PCs_CC';
            spaceDesc{end+1} = 'Common condition PCs';
    end    
    
    % collect sum of percent variance explained by all regression axes, per
    % condition
    dataTemp = sum(Summary.(otdr).(spaceFn).varAnalysis.V_RA_cond,2);
    % add row for mean
    dataTemp = [dataTemp; mean(dataTemp )];
    % compile across 3D spaces
    data = [data, dataTemp];
   
end

% generate list of condition labels
condDesc = {};
for i = 1:size(Summary.(otdr).Predictors,1)
    condDesc{end+1} = sprintf('%s = %d, %s = %d',...
        condLabel{1},...
        Summary.(otdr).Predictors(i,1),...
        condLabel{2},...
        Summary.(otdr).Predictors(i,2));
end
% add summary line
condDesc{end+1} = 'Mean across conditions';

% plot bars
h = bar(data);
% relabel
set(gca,'XTickLabel',condDesc,'XTickLabelRotation',-45);
legend(spaceDesc,'box','off','FontSize',16);
set(gca,'box','off','FontSize',16);
xlabel('Condition','FontSize',16);
ylabel('Absolute variance across time','FontSize',16);

%%% add annotation
annotateFig(figName);


%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot absolute temporal variance per neuron vs. cross-condition variance
% per neuron
figure;
whitebg(gcf,bgColor)
figName = ['Temporal vs. cross-condition variance - ',analName];
set(gcf,'Color',bgColor,'Name',figName);

% loop for absolute and z-score variance
for i=1:2
    if i == 1
        %%% ABSOLUTE variance
        % transform cross-condition and cross-time variance from z-scores to
        % absolute variance in Hz. A
        condVar = bsxfun(@times,Summary.(otdr).meta.preprocessingSummary.var_acrossCond,...
            1./Summary.(otdr).meta.preprocessingSummary.nrm.^2);
        tempVar = bsxfun(@times,Summary.(otdr).meta.preprocessingSummary.var_acrossTime,...
            1./Summary.(otdr).meta.preprocessingSummary.nrm'.^2);
        
        titleTxt = 'Variance of absolute firing rate';
        unitTxt = 'Hz';
    elseif i == 2
        %%% z-score variance
        % collect z-scored data 
        condVar = Summary.(otdr).meta.preprocessingSummary.var_acrossCond;
        tempVar = Summary.(otdr).meta.preprocessingSummary.var_acrossTime;

        titleTxt = 'Variance of z-scored firing rate';
        unitTxt = 's.d.';
    end
    
    % for cross-condition and cross-time variance and take average across
    % time and conditions, respectively
    condVar = median(condVar)';
    tempVar = median(tempVar,2);
    
    subplot(1,2,i)
    [h,hL,~,~,hP] = plotScatterN(tempVar, condVar,...
        'alpha',alpha,'statTest','signrank','bitEqualAxes',1,...
        'pTextFontSize',14,'color',textColor);
    hold on;
    
    % add mean and median
    text(0.05,1,sprintf('Temporal - condition variance:\nmedian = %0.2f \nmean = %0.2f',...
        median(tempVar - condVar),mean(tempVar - condVar)),'FontSize',16,...
        'Color',textColor,'VerticalAlignment','top','Units','normalized');
    
    title(titleTxt);
    set(h,'Marker','o','MarkerSize',8);
    set(hL,'Color',lineColor);
    set(hP,'FontSize',14);
    
    set(gca,'box','off','FontSize',16);
    xlabel(sprintf('Temporal variance, median across conditions (%s)',unitTxt),'FontSize',16);
    ylabel(sprintf('Cross-condition variance, median across time (%s)',unitTxt),'FontSize',16);
end

%%% add annotation
annotateFig(figName);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ANGLES between CC PCs and sRAs
nSRA = size(Summary.oTDRSummary.sRA.RA,2);
fprintf('\nSUBSPACE ANGLE between sRAs and top %d data PCs = %0.1f deg\n',...
    nSRA,rad2deg(subspace(Summary.oTDRSummary.PCs.RA(:,1:nSRA),Summary.oTDRSummary.sRA.RA)));
fprintf('SUBSPACE ANGLE between sRAs and top %d CC PCs = %0.1f deg\n',...
    nSRA,rad2deg(subspace(Summary.oTDRSummary.PCs_CC.RA(:,1:nSRA),Summary.oTDRSummary.sRA.RA)));
fprintf('SUBSPACE ANGLE between top %d data PCs and top %d CC PCs = %0.1f deg\n',...
    nSRA,nSRA,rad2deg(subspace(Summary.oTDRSummary.PCs_CC.RA(:,1:nSRA),Summary.oTDRSummary.PCs.RA(:,1:nSRA))));

fprintf('\nALIGNMENT INDEX between sRAs and top %d data PCs = %0.2f \n',...
    nSRA,align_ix(Summary.oTDRSummary.PCs.RA(:,1:nSRA),Summary.oTDRSummary.sRA.RA));
fprintf('ALIGNMENT INDEX between sRAs and top %d CC PCs = %0.2f \n',...
    nSRA,align_ix(Summary.oTDRSummary.PCs_CC.RA(:,1:nSRA),Summary.oTDRSummary.sRA.RA));
fprintf('ALIGNMENT INDEX between top %d data PCs and top %d CC PCs = %0.2f \n',...
    nSRA,nSRA,align_ix(Summary.oTDRSummary.PCs_CC.RA(:,1:nSRA),Summary.oTDRSummary.PCs.RA(:,1:nSRA)));

    
%% ANGLES

% whether to use previous offer TDR
bitPP = false;

% whether to use previous offer conditions on present trial responses
% (superceded by bitPP)
bitPRPC = false;

% set whether to plot autocorrelograms or cross correlograms
bitAuto = 1;

% optional string to append to figure titles
strTitleAppend = '';

bitReflectTheta = 1; % logical on whether to reflect theta about 90deg. 
        % Generally, if we reflect theta, then we should use the
        % "main_step" (cis) fits. If we do not reflect theta, then we should use
        % the "main_step_trans" fits.
bitRestrictRange = 1; % logical on whether to restrict range of theta to
        % 90 to 180 deg. This only applies when bitReflectTheta = false and
        % is used to isolate visualization of large angles.
bit2Step = false; % logical on whether to use fit from 2-step function. 
        
% Set type of p-value to use for pairwise angles
pVal2Use4Pair = 'Pe'; %'Pe' for empirical p-value; 'P' for model-based CDF
% Set type of p-value to use for stability metric
pVal2Use = 'Pe'; %'Pe' for empirical p-value; 'P' for model-based CDF

% when plotting stability metric heatmap, we can either plot the magnitude
% of the stable period (bitStabilityMag = 1) or the associated p-value
% (bitStabilityMag = 0). When plotting the magnitude, we can either
% mask-out the non-significant rows (bitStabilityMask = 1) or plot all
% values (bitStabilityMask = 0)
bitStabilityMag = 1;
bitStabilityMask = 1;

% Logical on whether to plot error of fit to surrogate distribution for
% pairwise angles
bitPlotDistFitErrPair = true;

% Logical on whether to plot error of fit to surrogate distribution for
% steps
bitPlotDistFitErrStep = true;

% Logical on whether to plot histograms of null theta with fits
bitPlotDistFit = false;

%%%%%%%%%%%
% set TDR field to use depending of whether using present or previous offer
tdrTemp = tdr;
if bitPP
    tdrTemp = [tdrTemp,'_PP'];
    
elseif bitPRPC
    tdrTemp = [tdrTemp,'_presRespPrevCond'];
end

% set which type of angle fit to use for stability metric
% Generally, if we reflect theta, then we should use the "main_step" TDRAlignDataAndRegressTimes (cis)
% fits. If we do not reflect theta, then we should use the
% "main_step_trans" fits.
if bitReflectTheta
    % 'main_decay' for exponential decay fits
    % 'main_step' for step function fits
    % 'main_step_2' for 2-step function fits
    % 'main_step_trans' for step function TRANS fit (angles that are anti-similar)
    s2u = 'main_step'; 
else
    s2u = 'main_step_trans'; 
end

% if using 2 step function, add "_2" suffix
if bit2Step
    s2u = [s2u,'_2'];
end
% append for using previous offer
if bitPP
    s2u = [s2u,'_PP'];
    analName = 'Previous trial';
elseif bitPRPC
    s2u = [s2u,'_presRespPrevCond'];
    analName = 'Prevoius trial conditions';
else
    analName = 'Present trial';
end

% Define upper limit of color range for p-value
maxPValStab = 0.01;
maxPValPair = 0.05;

% Define lower limit of color range for p-value depending on
% analysis, but common to the two animals
if bitAuto
    minPValStabPair = 10^-3; % for pairwise analysis of angles
    minPValStab = 10^-4; % for step-function
else
    minPValStabPair = 10^-3; % for pairwise analysis of angles
    minPValStab = 10^-3; % for step-function
end

% set alpha for stability metric to plot as bars alongside heatmaps. Set to
% empty to not plot bars.
if isfield(Summary.angleProfile(1).(s2u),sprintf('step_theta%s',pVal2Use))
    switch animalName
        case 'Norris'
            % define level to use for masking-out non-signficicant steps
            % per bitStabilityMask. When EMPTY, NO STABILITY panels are
            % plotted
            alphaStability = 10^-2.5;
            
            % Define lower limit of color range for p-value.
%             minPValStab = 10^-4;
            
            % tolerance for error in pairwise angle surrogate fits. It
            % needs to be high for monkey N, autocorrelogram, pairwise
            % angles.
            errThreshSigTestPairwise = 0.4;
            
        case 'Kirby'
            % define level to use for masking-out non-signficicant steps
            % per bitStabilityMask. When EMPTY, NO STABILITY panels are
            % plotted
            alphaStability = 10^-3;
            % Define lower limit of color range for p-value.
%             minPValStab = 10^-3;

            % tolerance for error in pairwise angle surrogate fits
            errThreshSigTestPairwise = 0.11;
        otherwise
            error('Animal name %s not recognized',animalName);
    end
else
    alphaStability = [];
    minPValStab = [];
    
end

%%% Override alphastability with the upper bound on p-values
alphaStability = maxPValStab;

% set mask color for regions of heat map without data
% and set colormap used for probability scale
if strcmp(bgColor,'w')
    maskColor = [0.7 0.7 0.7];
    nonSigColor = [0.3 0.3 0.3];
    pColormap = 'copper';
else
	maskColor = [.3 .3 .3];
    nonSigColor = [0.7 0.7 0.7];
    pColormap = 'summer';
end

% set colormap used for angle scale
aColormap = 'jet';

%%%%%%%%%%%%%%%%%%
strAnglePair = {};
strTitle = {};

% extract signal numbers
sigN = NaN(length(Summary.angleProfile),2);
bitTranspose = NaN(length(Summary.angleProfile),1);
for i = 1:length(Summary.angleProfile)
    sigN(i,:) = Summary.angleProfile(i).(s2u).signalN;
    bitTranspose(i) = Summary.angleProfile(i).(s2u).bitTranspose;
end

if bitAuto 
    strAngleType = 'Auto';
    
    sigN2Plot = find(diff(sigN,1,2)==0)';
%     for i = 1:3
%         strAnglePair{i} = [num2str(i),num2str(i)];
%         switch i
%             case 1
%                 strTitle{i} = 'Benefit';
%             case 2
%                 strTitle{i} = 'Choice';
%             case 3
%                 strTitle{i} = 'Expected Reward';
%             otherwise
%                 error('number of conditions exceeds max');
%         end
%     end
else
    strAngleType = 'Cross';
    
    sigN2Plot = find(diff(sigN,1,2)~=0 & ~bitTranspose)';
    
%     for i = 1:3
%         switch i
%             case 1
%                 strAnglePair{i} = '12';
%                 strTitle{i} = 'Benefit x Choice';
%             case 2
%                 strAnglePair{i} = '13';
%                 strTitle{i} = 'Benefit x Expected Reward';
%             case 3
%                 strAnglePair{i} = '23';
%                 strTitle{i} = 'Choice x Expected Reward';
%         end                
%     end
end
for i = sigN2Plot
    strTitle{end+1} = Summary.angleProfile(i).(s2u).signalDesc;
    strAnglePair{end+1} = [num2str(sigN(i,1)),num2str(sigN(i,2))];
end
   
figure;
whitebg(gcf,bgColor)
figName = [strAngleType,'-angle between regression vectors - ',analName,strTitleAppend];
set(gcf,'Color',bgColor,'Name',figName);

nCol = length(condLabel);
nRow = 2; % angles and probabilities
% add row for RSV magnitude for auto angles
if bitAuto
    nRow = nRow + 1;
end
% add row for stability metric, if available
if ~isempty(alphaStability)
    nRow = nRow + 1;
end
% add row for transposed stability metric
if ~bitAuto & any(diff(sigN,1,2)~=0 & bitTranspose)
    nRow = nRow + 1;
end
% add row for surrogate distribution fit error for pairwise angles
if bitPlotDistFitErrPair && strcmp('P',pVal2Use4Pair)
    nRow = nRow + 1;
end
% add row for surrogate distribution fit error for steps
if bitPlotDistFitErrStep && strcmp('P',pVal2Use)
    nRow = nRow + 1;
end

axN = 0;

% define range of p-value for stability metric
if ~isempty(alphaStability)
    cBarRangeStab = log10([minPValStab maxPValStab]);
end

% find min angle
minAngle = NaN(length(condLabel),1);
maxAngle = NaN(length(condLabel),1);
for i = 1:length(condLabel);

    foo = Summary.(tdrTemp).angleAnalysis.(['angle',strAnglePair{i}]);
    
    % set diagonal to 0
    foo = foo + eye(size(foo,1))*90 - diag(diag(foo));

    % transform theta such that values > 90 are reflected about 90deg.
    if bitReflectTheta
        foo = foo - 2*max(foo-90,0);
    end
    
    minAngle(i) = min(foo(:));
    maxAngle(i) = max(foo(:));
end

% round min/max angle to nearest multiple of 5
goo = [0:5:180];
if bitRestrictRange && ~bitReflectTheta
    minAngle4Plot = 90;
else
    foo = find(min(minAngle)>goo,1,'last');
    minAngle4Plot = goo(foo);
end
foo = find(max(maxAngle)<=goo,1,'first');
maxAngle4Plot = goo(foo);
clear goo foo

% set color ranges for plotting boxcar magnitude vs. p-value
mag.cBarRange = [minAngle4Plot maxAngle4Plot];
mag.cBarTick = [minAngle4Plot:15:maxAngle4Plot];
mag.cBarTick(end) = maxAngle4Plot;
magP.cBarRange = cBarRangeStab;
magP.cBarTick = [log10(minPValStab):1:log10(maxPValStab),log10(maxPValStab)];
% add point for alpha stability threshold, and use unique
% in case alphaStability is replicated
%             cBarTick = unique([cBarTick log10(alphaStability)]);
magP.cBarTick = unique(magP.cBarTick);
magP.cBarTick(end) = log10(maxPValStab);

%%% PLOT ANGLES
t = Summary.(tdrTemp).angleAnalysis.regressTimes;
dt = t(2) - t(1);
tRange = range(t);
for i = 1:length(condLabel)
    
    axN = axN + 1;
    axH = subplot(nRow,nCol,axN);
        
    foo = Summary.(tdrTemp).angleAnalysis.(['angle',strAnglePair{i}]);
    
%     % set diagonal to 0
%     foo = foo + eye(size(foo,1))*90 - diag(diag(foo));
    
    % transform theta such that values > 90 are reflected about 90deg.
    if bitReflectTheta
        foo = foo - 2*max(foo-90,0);
    end
    
    % set diagonal to NaN for auto angles
    if bitAuto
        foo(logical(eye(size(foo,1)))) = NaN;
    end
    
    % PLOT heatmap
    h = imagesc(t,t,foo,[minAngle4Plot maxAngle4Plot]);
    
    % set NaNs to gray
    set(h,'AlphaData',~isnan(foo));
    set(gca,'Color',maskColor)

    axis square
    title(strTitle{i},'FontSize',16);
    set(gca,'FontSize',14);
    
    h = colorbar('FontSize',14);

    if i ~= nCol
        set(h,'Visible','off')
    else
        goo = [minAngle4Plot:15:maxAngle4Plot];
        goo(end) = maxAngle4Plot;
        set(h,'YTick',goo)
    end    
        
    % colormap -- NOTE: we used to do this after creating the axes but
    % before plotting the heatmap. With MATLAB 2019a, it appears that
    % plotting the heatmap resets the colormap. So now we do this inversion
    % afterwards. I'm not sure what the behavior will now be with earlier
    % versions of MATLAB; it may require performing these steps before
    % plotting.
    colormap(gca,aColormap);
    
    % invert colormap only for reflected theta (so that small angles are
    % warm colors) 
    if bitReflectTheta
        goo = colormap(gca);
        colormap(gca,flipud(goo));
        clear goo
    end
    
    % add sidebar to show stability metric p-value for each dRAs 
    if ~isempty(alphaStability)
        % extract pValues for steps:
        pVal = Summary.angleProfile(sigN2Plot(i)).(s2u).(sprintf('step_theta%s',pVal2Use));
        % in some cases there are multiple p-values per row (1 for each
        % step), but we need a single p-value. In these cases take the 
        % smallest pVal
        pValOnly1 = min(pVal,[],2);

        % at times, there will be values of p=0. change these to min
        % p-value for display purposes
        pValOnly1(pValOnly1==0) = min(10.^cBarRangeStab);
        
        if bitStabilityMask
            % limit to the significant time points (rows)
            tSig = t(pValOnly1 <= alphaStability);
            pValOnly1(pValOnly1 > alphaStability | ~isfinite(pValOnly1)) = [];
        else
            tSig = t;
        end
        % for vertical lines
%         tSig = [tSig-dt/2; tSig+dt/2];
%         x = (t(1)-0.2*tRange) * ones(size(tSig));
        % for horizontal lines
        tSig = [tSig; tSig];
        x = [(t(1)-0.26*tRange) * ones(1,size(tSig,2));
            (t(1)-0.22*tRange) * ones(1,size(tSig,2))];
        
        % define colors for each p-value
        colorSet = eval(pColormap);
        % invert colormap for black background
        if strcmp(bgColor,'k')
            colorSet = flipud(colorSet);
        end    
        
        % build evenly space set of pvals
        pValSet = linspace(cBarRangeStab(1),cBarRangeStab(2),size(colorSet,1));
        % interpolate p-values to the set of pVal, finding the position of
        % the nearest pVal in the set:
        pos = interp1(pValSet,[1:length(pValSet)],...
            log10(pValOnly1)',...
            'nearest','extrap');
        
        hL = line(x,tSig,'LineWidth',6,'Clipping','off');
        for j = 1:length(hL)
            if isfinite(pos(j))
                set(hL(j),'Color',colorSet(pos(j),:));
            else
                % if p-value is not defined (i.e., no step was fit), then
                % hide line
                set(hL(j),'Visible','off');
            end
        end
        clear bitSig
                
        % For cross-signal angles, we need the "vertical" steps, that is,
        % steps fit to the transposed matrix of theta
        bitCousin = ismember(sigN,fliplr(Summary.angleProfile(sigN2Plot(i)).(s2u).signalN),'rows');
        if length(unique(Summary.angleProfile(sigN2Plot(i)).(s2u).signalN))>1 && ...
                ~Summary.angleProfile(sigN2Plot(i)).(s2u).bitTranspose && ...
                any(bitCousin)
            
            % find signalN that is the transposed version of current signal
            % number
            sigPos = find(bitCousin);
            
            % extract pValues for steps:
            pVal = Summary.angleProfile(sigPos).(s2u).(sprintf('step_theta%s',pVal2Use));
            % in some cases there are multiple p-values per row (1 for each
            % step), but we need a single p-value. In these cases take the
            % smallest pVal
            pValOnly1 = min(pVal,[],2);
            % at times, there will be values of p=0. change these to min
            % p-value for display purposes
            pValOnly1(pValOnly1==0) = min(10.^cBarRangeStab);
            
            
            if bitStabilityMask
                % limit to the significant time points (rows)
                tSig = t(pValOnly1 <= alphaStability);
                pValOnly1(pValOnly1 > alphaStability | ~isfinite(pValOnly1)) = [];
            else
                tSig = t;
            end
            % for vertical lines
            tSig = [tSig; tSig];
            y = [(t(end)+0.28*tRange) * ones(1,size(tSig,2));
                (t(end)+0.24*tRange) * ones(1,size(tSig,2))];
            
            % interpolate p-values to the set of pVal, finding the position of
            % the nearest pVal in the set:
            pos = interp1(pValSet,[1:length(pValSet)],...
                log10(pValOnly1)',...
                'nearest','extrap');
            
            hL = line(tSig,y,'LineWidth',6,'Clipping','off');
            for j = 1:length(hL)
                if isfinite(pos(j))
                    set(hL(j),'Color',colorSet(pos(j),:));
                else
                    % if p-value is not defined (i.e., no step was fit), then
                    % hide line
                    set(hL(j),'Visible','off');
                end            
            end
            clear bitSig
            
        end

        if bitStabilityMag && i==1
            % when plotting stability magnitude, make a fake axes with
            % pvalues so as to make p-value colorbar used to decode the
            % p-value sidebar
            
            % make new axes in same place as existing axes
            axTemp = copyobj(axH,gcf);
            % remove image from figure
            delete(get(axTemp,'Children'));
            
            % make fake image with correct range
            hTemp = imagesc(axTemp,t,t,[magP.cBarRange;magP.cBarRange],magP.cBarRange);
            axis square
            
            % set color map to probabilities -- NOTE: with MATLAB 2019a,
            % this must now be done after imagesc(). See above note.
            if strcmp(bgColor,'k')
                % flip order when on dark background
                foo = eval(pColormap);
                colormap(axTemp,flipud(foo));
                clear foo
            else
                colormap(axTemp,pColormap);
            end            
            
            hCB = colorbar(axTemp,'FontSize',14);
            set(hCB,'YTick',magP.cBarTick)
            goo = get(hCB,'YTick');
            clear soo
            for j = 1:length(goo)
                soo{j} = sprintf('10^{%0.2g}',goo(j));
            end
            soo{end} = sprintf('%g',maxPValStab);
            set(hCB,'YTickLabel',soo);
            clear soo

            % font size
            set(axTemp,'FontSize',14);
            % make square
            axis(axTemp,'image');
            % hide stability p-value axes
            set(axTemp,'Visible','off')
            % set to transparent
            set(hTemp,'AlphaData',0);
            clear hTemp axTemp
            
            % put original axes back in front
            axes(axH);
            
        end
        
    end
end

%%% PLOT Stability metric if available
if ~isempty(alphaStability)
    clear foo
    for i = 1:length(condLabel);
        
        axN = axN + 1;
        axH = subplot(nRow,nCol,axN);
                
        if bitStabilityMag
            cBarRange = mag.cBarRange;
            cBarTick = mag.cBarTick;
        else
            cBarRange = magP.cBarRange;
            cBarTick = magP.cBarTick;
        end
        
        % extract stability
        foo = NaN(size(Summary.angleProfile(sigN2Plot(i)).(s2u).b,1)); % for boxcars to be colored
        fooNon = NaN(size(Summary.angleProfile(sigN2Plot(i)).(s2u).b,1)); % for boxcars to be gray
        bitSig = Summary.angleProfile(sigN2Plot(i)).(s2u).(sprintf('step_theta%s',pVal2Use)) <= alphaStability;
        % loop through rows
        for j = 1:length(bitSig)
%             if bitStabilityMask && bitStabilityMag && ~bitSig(j)
%                 continue
%             end
            % loop through multiple steps
            for st = 1:size(Summary.angleProfile(sigN2Plot(i)).(s2u).b,3)
                bitStep = t >= Summary.angleProfile(sigN2Plot(i)).(s2u).b(j,2,st) & t <= sum(Summary.angleProfile(sigN2Plot(i)).(s2u).b(j,2:3,st),2);
                if bitStabilityMag
                    if bitStabilityMask && ~bitSig(j,st)
                        % if non-significant boxcar, color gray
                        fooNon(j,bitStep) = 0;
%                         foo(j,bitStep) = -Inf;
                    else
                        % fill in elements with magnitude of step, converted back to
                        % original theta from theta*:
                        % theta* = angleScale*(angleOffset - theta)
                        % theta = angleOffset - (theta* / angleScale)
                        foo(j,bitStep) = Summary.angleProfile(sigN2Plot(i)).(s2u).angleOffset - ...
                            Summary.angleProfile(sigN2Plot(i)).(s2u).b(j,1,st) / Summary.angleProfile(sigN2Plot(i)).(s2u).angleScale;
                    end
                else
                    % fill in elements with p-value
                    foo(j,bitStep) = log10(Summary.angleProfile(sigN2Plot(i)).(s2u).(sprintf('step_theta%s',pVal2Use))(j,st));
                end
            end
        end
        
        % PLOT
        h = imagesc(t,t,foo,cBarRange);
        
        % use angle color map for angles -- NOTE: With MATLAB 2019a, we now
        % need to set colormap after imagesc(). See note above.        
        if bitStabilityMag
            colormap(gca,aColormap);
            % invert colormap only for reflected theta (so that small angles are
            % warm colors)
            if bitReflectTheta
                colorMapTemp = flipud(colormap(gca));
                colormap(gca,colorMapTemp);
            end
        else
            % invert probability colormap for black background
            if strcmp(bgColor,'k')
                foo = eval(pColormap);
                colormap(gca,flipud(foo));
                clear foo
            else
                % probability colormap
                colormap(gca,pColormap);
            end
        end
        
        % set NaNs to gray
        set(h,'AlphaData',~isnan(foo));
        set(gca,'Color',maskColor);
                
        axis square
        set(gca,'FontSize',14);
        
        h = colorbar('FontSize',14);
        if i ~= nCol
            set(h,'Visible','off')
        else
            set(h,'YTick',cBarTick)

            % special plotting for probabilities
            if ~bitStabilityMag
                goo = get(h,'YTick');
                clear soo
                for j = 1:length(goo)
                    soo{j} = sprintf('10^{%0.2g}',goo(j));
                end
                soo{end} = sprintf('%g',maxPValStab);
                set(h,'YTickLabel',soo);
                clear soo
            end
        end
        
        clear foo bitSig

        if bitStabilityMag && i==1
            % when plotting stability magnitude, make a fake axes with
            % pvalues so as to make p-value colorbar used to decode the
            % p-value sidebar
            
            % make new axes in same place as existing axes
            axTemp = copyobj(axH,gcf);
            % link axes so they stay the same size
            linkprop([axH,axTemp],{'Position'});
            % remove image from figure
            delete(get(axTemp,'Children'));

            
            % make fake image with correct range
            hTemp = imagesc(axTemp,t,t,[magP.cBarRange;magP.cBarRange],magP.cBarRange);
            
            % set color map to probabilities -- NOTE: with MATLAB 2019a,
            % this must now be done after imagesc(). See above note.
            if strcmp(bgColor,'k')
                % flip order when on dark background
                foo = eval(pColormap);
                colormap(axTemp,flipud(foo));
                clear foo
            else
                colormap(axTemp,pColormap);
            end            
            
            hCB = colorbar(axTemp,'FontSize',14);
            set(hCB,'YTick',magP.cBarTick)
            goo = get(hCB,'YTick');
            clear soo
            for j = 1:length(goo)
                soo{j} = sprintf('10^{%0.2g}',goo(j));
            end
            soo{end} = sprintf('%g',maxPValStab);
            set(hCB,'YTickLabel',soo);
            clear soo
            
%             % font size
%             set(axTemp,'FontSize',14);
            % make square
            axis(axTemp,'square')            
            % hide stability p-value axes
            set(axTemp,'Visible','off')
            % set to transparent
            set(hTemp,'AlphaData',0);
            clear hTemp axTemp

            % put original axes back in front
            axes(axH);
            
        end 
        
        % add new axis to show non-significant boxcars 
        if any(~isnan(fooNon(:)))
            axNS = copyobj(axH,gcf);
            hN = imagesc(axNS,t,t,fooNon);
            colormap(axNS,[nonSigColor]);
            % set NaNs to transparent
            set(hN,'AlphaData',~isnan(fooNon));
            axis(axNS,'square');
            % link axes so they stay the same size
            linkprop([axH,axNS],{'Position'});            
            
%             % add fake color bar
%             hTemp = colorbar;
%             set(hTemp,'Visible','off');
            % turn off axes
            set(axNS,'Visible','off');
        end        
        
        clear fooNon
        
        %%%%%% When plotting stability magnitude
        % add sidebar to show stability metric p-value for each dRAs
        if ~isempty(alphaStability) && bitStabilityMag
            % extract pValues for steps:
            pVal = Summary.angleProfile(sigN2Plot(i)).(s2u).(sprintf('step_theta%s',pVal2Use));
            % in some cases there are multiple p-values per row (1 for each
            % step), but we need a single p-value. In these cases take the
            % smallest pVal
            pValOnly1 = min(pVal,[],2);
            % at times, there will be values of p=0. change these to min
            % p-value for display purposes
            pValOnly1(pValOnly1==0) = min(10.^cBarRangeStab);
            
            if bitStabilityMask
                % limit to the significant time points (rows)
                tSig = t(pValOnly1 <= alphaStability);
                pValOnly1(pValOnly1 > alphaStability | ~isfinite(pValOnly1)) = [];
            else
                tSig = t;
            end
            % for vertical lines
            %         tSig = [tSig-dt/2; tSig+dt/2];
            %         x = (t(1)-0.2*tRange) * ones(size(tSig));
            % for horizontal lines
            tSig = [tSig; tSig];
            x = [(t(1)-0.26*tRange) * ones(1,size(tSig,2));
                (t(1)-0.22*tRange) * ones(1,size(tSig,2))];
            
            % define colors for each p-value
            colorSet = eval(pColormap);
            % invert colormap for black background
            if strcmp(bgColor,'k')
                colorSet = flipud(colorSet);
            end
            
            % build evenly space set of pvals
            pValSet = linspace(cBarRangeStab(1),cBarRangeStab(2),size(colorSet,1));
            % interpolate p-values to the set of pVal, finding the position of
            % the nearest pVal in the set:
            pos = interp1(pValSet,[1:length(pValSet)],...
                log10(pValOnly1)',...
                'nearest','extrap');
            
            hL = line(x,tSig,'LineWidth',6,'Clipping','off');
            for j = 1:length(hL)
                if isfinite(pos(j))
                    set(hL(j),'Color',colorSet(pos(j),:));
                else
                    % if p-value is not defined (i.e., no step was fit), then
                    % hide line
                    set(hL(j),'Visible','off');
                end
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If cross-correlogram in need of supplemental transposed data, include
    % here
    
    % loop again through signals
    for i = 1:length(condLabel);
        
        bitCousin = ismember(sigN,fliplr(Summary.angleProfile(sigN2Plot(i)).(s2u).signalN),'rows');
        if length(unique(Summary.angleProfile(sigN2Plot(i)).(s2u).signalN))>1 && ...
                ~Summary.angleProfile(sigN2Plot(i)).(s2u).bitTranspose && ...
                any(bitCousin)
            clear foo
            
            % find signalN that is the transposed version of current signal
            % number
            sigPos = find(bitCousin);
            if length(sigPos) ~= 1
                error('More or less than 1 match to transposed stability analysis was found')
            end
            if ~bitTranspose(sigPos)
                error('Stability analysis found is not flagged as transposed')
            end
            
            
            axN = axN + 1;
            axH = subplot(nRow,nCol,axN);
            
            % set some ranges depending on plotting magnitude vs. p-value
            if bitStabilityMag
                cBarRange = [minAngle4Plot maxAngle4Plot];
                cBarTick = [minAngle4Plot:15:maxAngle4Plot];
                cBarTick(end) = maxAngle4Plot;
            else
                cBarRange = cBarRangeStab;
                cBarTick = [log10(minPValStab):5:log10(maxPValStab),log10(maxPValStab)];
                % add point for alpha stability threshold, and use unique
                % in case alphaStability is replicated
                cBarTick = unique([cBarTick log10(alphaStability)]);
                cBarTick(end) = log10(maxPValStab);
            end
            
            % extract stability
            foo = NaN(size(Summary.angleProfile(sigPos).(s2u).b,1));
            fooNon = NaN(size(Summary.angleProfile(sigPos).(s2u).b,1)); % for boxcars to be gray            
            bitSig = Summary.angleProfile(sigPos).(s2u).(sprintf('step_theta%s',pVal2Use)) <= alphaStability;
            
            % loop through rows
            for j = 1:length(bitSig)
                % loop through multiple steps
                for st = 1:size(Summary.angleProfile(sigN2Plot(i)).(s2u).b,3)
                    bitStep = t >= Summary.angleProfile(sigPos).(s2u).b(j,2,st) & t <= sum(Summary.angleProfile(sigPos).(s2u).b(j,2:3,st),2);
                                        
                    if bitStabilityMag
                        if bitStabilityMask && ~bitSig(j,st)
                            % if non-significant boxcar, color gray
                            fooNon(j,bitStep) = 0;
                        else
                            % fill in elements with magnitude of step, converted back to
                            % original theta from theta*:
                            % theta* = angleScale*(angleOffset - theta)
                            % theta = angleOffset - (theta* / angleScale)
                            foo(j,bitStep) = Summary.angleProfile(sigPos).(s2u).angleOffset - ...
                                Summary.angleProfile(sigPos).(s2u).b(j,1,st) / Summary.angleProfile(sigPos).(s2u).angleScale;
                        end
                    else
                        % fill in elements with p-value
                        foo(j,bitStep) = log10(Summary.angleProfile(sigPos).(s2u).(sprintf('step_theta%s',pVal2Use))(j,st));
                    end
                end
            end
                
            % transpose so that bars run vertically:
            foo = foo';
            fooNon = fooNon';
            
            % plot image
            h = imagesc(t,t,foo,cBarRange);
            
            % use angle color map for angles -- NOTE: In MATLAB 2019a, this
            % must be done after imagesc(). See above note.
            if bitStabilityMag
                colormap(gca,aColormap);
                % invert colormap only for reflected theta (so that small angles are
                % warm colors)
                if bitReflectTheta
                    goo = colormap(gca);
                    colormap(gca,flipud(goo));
                    clear goo
                end
            else
                % probability colormap
                colormap(gca,pColormap);
            end
            % invert colormap for black background
            if strcmp(bgColor,'k')
                goo = colormap(gca);
                colormap(gca,flipud(goo));
                clear goo
            end
            
            % set NaNs to gray
            set(h,'AlphaData',~isnan(foo));
            set(gca,'Color',maskColor)
            
            axis square
            set(gca,'FontSize',14);
            
            %%%% When plotting stability magnitude
            if bitStabilityMag
                % For cross-signal angles, we need the "vertical" steps, that is,
                % steps fit to the transposed matrix of theta
                bitCousin = ismember(sigN,fliplr(Summary.angleProfile(sigN2Plot(i)).(s2u).signalN),'rows');
                if length(unique(Summary.angleProfile(sigN2Plot(i)).(s2u).signalN))>1 && ...
                        ~Summary.angleProfile(sigN2Plot(i)).(s2u).bitTranspose && ...
                        any(bitCousin)
                    
                    % find signalN that is the transposed version of current signal
                    % number
                    sigPos = find(bitCousin);
                    
                    % extract pValues for steps:
                    pVal = Summary.angleProfile(sigPos).(s2u).(sprintf('step_theta%s',pVal2Use));
                    % in some cases there are multiple p-values per row (1 for each
                    % step), but we need a single p-value. In these cases take the
                    % smallest pVal
                    pValOnly1 = min(pVal,[],2);
                    % at times, there will be values of p=0. change these to min
                    % p-value for display purposes
                    pValOnly1(pValOnly1==0) = min(10.^cBarRangeStab);
                    
                    
                    if bitStabilityMask
                        % limit to the significant time points (rows)
                        tSig = t(pValOnly1 <= alphaStability);
                        pValOnly1(pValOnly1 > alphaStability | ~isfinite(pValOnly1)) = [];
                    else
                        tSig = t;
                    end
                    % for vertical lines
                    tSig = [tSig; tSig];
                    y = [(t(end)+0.28*tRange) * ones(1,size(tSig,2));
                        (t(end)+0.24*tRange) * ones(1,size(tSig,2))];
                    
                    % interpolate p-values to the set of pVal, finding the position of
                    % the nearest pVal in the set:
                    pos = interp1(pValSet,[1:length(pValSet)],...
                        log10(pValOnly1)',...
                        'nearest','extrap');
                    
                    hL = line(tSig,y,'LineWidth',6,'Clipping','off');
                    for j = 1:length(hL)
                        if isfinite(pos(j))
                            set(hL(j),'Color',colorSet(pos(j),:));
                        else
                            % if p-value is not defined (i.e., no step was fit), then
                            % hide line
                            set(hL(j),'Visible','off');
                        end
                    end
                    clear bitSig
                    
                end
            end
            
            h = colorbar('FontSize',14);
            if i ~= nCol
                set(h,'Visible','off')
            else
                set(h,'YTick',cBarTick)
                
                % special plotting for probabilities
                if ~bitStabilityMag
                    goo = get(h,'YTick');
                    clear soo
                    for j = 1:length(goo)
                        soo{j} = sprintf('10^{%0.2g}',goo(j));
                    end
                    soo{end} = sprintf('%g',maxPValStab);
                    set(h,'YTickLabel',soo);
                    clear soo
                end
            end
            
            % add new axis to show non-significant boxcars
            if any(~isnan(fooNon(:)))
                axNS = copyobj(axH,gcf);
                hN = imagesc(axNS,t,t,fooNon);
                colormap(axNS,[nonSigColor]);
                % set NaNs to transparent
                set(hN,'AlphaData',~isnan(fooNon));
                axis(axNS,'square');
                % link axes so they stay the same size
                linkprop([axH,axNS],{'Position'});
                
                %             % add fake color bar
                %             hTemp = colorbar;
                %             set(hTemp,'Visible','off');
                % turn off axes
                set(axNS,'Visible','off');
            end
            
            % add title because inverted -- actually, skip this. The data
            % was originally inverted so as to fit horizontal bars. But
            % then the bars were inverted back to the original orientation
            % for plotting
%             title(Summary.angleProfile(sigPos).(s2u).signalDesc,'FontSize',12);
            
            clear foo fooNon bitSig
        end
    end
end

%%% PLOT error in fitting distribution to surrogate steps
if bitPlotDistFitErrStep && strcmp('P',pVal2Use)
    yLimErr = [];
    axNErr = [];
    for i = 1:length(condLabel);
        
        axN = axN + 1;
        axNErr(end+1) = subplot(nRow,nCol,axN);
        
        % plot error
        plot(Summary.angleProfile(sigN2Plot(i)).(s2u).tMaster,...
            Summary.angleProfile(sigN2Plot(i)).(s2u).step_thetaSDistFitErr,'^',...
            'Color',lineColor);
        hold on
        
        if i == 1
            ylabel('Fit error step','FontSize',16);
        end
        
        axis square
        set(gca,'FontSize',14);
        h = colorbar('FontSize',14); % so that axis is same width as above axes
        set(h,'Visible','off')
        
        yLimErr(end+1,:) = get(gca,'YLim');
    end
    
    % If cross-correlogram in need of supplemental transposed data, include
    % here
    % loop again through signals
    for i = 1:length(condLabel)
        
        bitCousin = ismember(sigN,fliplr(Summary.angleProfile(sigN2Plot(i)).(s2u).signalN),'rows');
        if length(unique(Summary.angleProfile(sigN2Plot(i)).(s2u).signalN))>1 && ...
                ~Summary.angleProfile(sigN2Plot(i)).(s2u).bitTranspose && ...
                any(bitCousin)
            clear foo
                       
            % find signalN that is the transposed version of current signal
            % number
            sigPos = find(bitCousin);
            if length(sigPos) ~= 1
                error('More or less than 1 match to transposed stability analysis was found')
            end
            if ~bitTranspose(sigPos)
                error('Stability analysis found is not flagged as transposed')
            end
            
            % plot error
            plot(axNErr(i),Summary.angleProfile(sigPos).(s2u).tMaster,...
                Summary.angleProfile(sigPos).(s2u).step_thetaSDistFitErr,'v',...
                'Color',lineColor);
            
            yLimErr(end+1,:) = get(gca,'YLim');
            
        end
    end
    
    
    set(axNErr,'YLim',[min(yLimErr(:,1)) max(yLimErr(:,2))]);
    set(axNErr,'XLim',Summary.angleProfile(sigN2Plot(1)).(s2u).tMaster([1 end]));
end


%%% PLOT RSV of dRA at each timepoint for autocorrelograms only
if bitAuto
    clear foo
    for i = 1:length(condLabel);
        
        axN = axN + 1;
        subplot(nRow,nCol,axN);
        
        % get times
        tData = round(Summary.(tdrTemp).Times.all_times,3);
        tRegress = round(Summary.(tdrTemp).Times.regressTimes,3);
        dTRegress = mean(diff(tRegress));
        
        % accommodate legacy Summary struct in which RSV for the 3 signals
        % (not dRAs) were compiled in the second dimension
        [T dim2 nRA] = size(Summary.(tdrTemp).varAnalysis.RSVE_RA);
        rT = size(Summary.(tdrTemp).dRAs,1);
        nSignal = length(condLabel);
        if ndims(Summary.(tdrTemp).varAnalysis.RSVE_RA) == 3 ...
                && dim2 == rT * nSignal
            foo = reshape(Summary.(tdrTemp).varAnalysis.RSVE_RA,...
                T,rT,nRA,nSignal);
        else
            foo = Summary.(tdrTemp).varAnalysis.RSVE_RA;
        end
        % extract RSV where data times = regression times
        RSV = TDRAlignDataAndRegressTimes(foo,i,...
            'T',tData,'RT',tRegress,...
            'nBinRegress',Summary.(tdrTemp).Times.regressBins);
        clear foo
        
        % accommodate legacy Summary struct in which VE for the 3 signals
        % (not dRAs) were compiled in the second dimension
        [T dim2] = size(Summary.(tdrTemp).varAnalysis.VE_RA);
        if ndims(Summary.(tdrTemp).varAnalysis.VE_RA) == 2 ...
                && dim2 == rT * nSignal
            foo = reshape(Summary.(tdrTemp).varAnalysis.VE_RA,T,rT,nRA);
        else 
            foo = Summary.(tdrTemp).varAnalysis.VE_RA;
        end
        % convert to r^2, where r = corr(data,predictor), and RSV =
        % Variance * r^2.
        r2 = RSV ./ TDRAlignDataAndRegressTimes(foo,i,...
            'T',tData,'RT',tRegress,...
            'nBinRegress',Summary.(tdrTemp).Times.regressBins);
        clear foo
        % extract norm of dRA
        dRANorm = Summary.(tdrTemp).normdRAs(:,i);

        %%%%%%
        %%%%%%
        % The following code was necessary when we had to align data times
        % and regression times
%         tRSV = round(Summary.(tdrTemp).Times.all_times,3);
%         tdRA = round(Summary.(tdrTemp).Times.regressTimes,3);
%         dt_RSV = roundDec(tRSV(2) - tRSV(1),10);
%         dt_dRA = roundDec(tdRA(2) - tdRA(1),10);
%         tWin_RSV = [tRSV(1) - dt_RSV/2, tRSV(end) + dt_RSV/2];
%         tWin_dRA = [tdRA(1) - dt_dRA/2, tdRA(end) + dt_dRA/2];
%         
%         % eliminate RSV times that are outside of dRA times:
%         bitTime2Elim = tRSV < tWin_dRA(1) | tRSV > tWin_dRA(2);
%         foo(bitTime2Elim,:) = [];
%         tRSV(bitTime2Elim) = [];
%         
%         % if dRA and RSV time bins aren't equal, then we have to take average
%         % RSV across timebins corresponding to single dRA timebin
%         if ~isequal(dt_RSV,dt_dRA)
%             if dt_RSV > dt_dRA
%                 error('time bin for RSV %f cannot be greater than time bin for dRA %f',dt_RSV,dt_dRA)
%             elseif mod(dt_dRA,dt_RSV) ~= 0
%                 error('time bin for dRA must be an integer multiple of time bine for RSV')
%             else
%                 % the most careful way to do this is to loop through each dRA
%                 % time bin and find corresponding RSV time bins
%                 RSV = NaN(size(tdRA));
%                 
%                 for j = 1:length(tdRA)
%                     bitFoo = roundDec(tRSV-dt_RSV/2,10) >= roundDec(tdRA(j)-dt_dRA/2,10) & ...
%                         roundDec(tRSV+dt_RSV/2,10) <= roundDec(tdRA(j)+dt_dRA/2,10);
%                     
%                     % make sure number of matches was as expected
%                     if sum(bitFoo) ~= roundDec(dt_dRA/dt_RSV,10)
%                         error('Number of matching time bins for dRA #%d (%d) was not as expected (%d)',...
%                             j,sum(bitFoo),dt_dRA/dt_RSV)
%                     end
%                     % mean
%                     RSV(j) = mean(foo(bitFoo,j));
%                 end
%             end
%         else
%             RSV = diag(foo);
%         end
        
%         plot(tData,RSV,'-','Color',lineColor)
%         plot(tData,r2,'-','Color',lineColor)
        plot(tRegress,dRANorm,'-','Color',lineColor)
        axis square
        
        % limit to regression times
        set(gca,'FontSize',14,'XLim',[tRegress(1)-dTRegress/2, tRegress(end)+dTRegress/2],...
            'box','off'); %'xTick',[tWin_dRA(1):2:tWin_dRA(2)]
        
        h = colorbar('FontSize',14); % so that axis is same width as above axes
        set(h,'Visible','off')
        
        xlabel('Time from offer (s)','FontSize',14);
        
        if i == 1
            ylabel('dRA magnitude (a.u.)','FontSize',14);
        end
        
        clear RSV tData nBinRegress
    end
end


%%% PLOT PROBABILITIES
minPVal = NaN;

% we used to compute miPVal based on number of entries, i.e., Bonforonni
% correct:
% clear foo
% foo = Summary.(tdrTemp).angleAnalysis.(['pVal',strAnglePair{1}]);
% minPVal = 10^floor(log10(alpha / (size(foo,1)^2 - size(foo,1)) / 2)); %   %min(foo(:));
clear foo distFit pVal
for i = 1:length(condLabel);

    %%% Here we compute the probabilities denovo based on the surrogate
    %%% data and ignore the precomputed field
    %%% TDRSummary.angleAnalysis.pValij. The reason is that we want the
    %%% null distribution to adjust depending on the type of angles we're
    %%% plotting: folded vs. unfolded, restricted vs. full range.
    
    % extract angles of interest for real and surrogate data
    foo = Summary.(tdrTemp).angleAnalysis.(['angle',strAnglePair{i}]);
    fooS = Summary.(tdrTemp).angleAnalysis.(['surrAngle',strAnglePair{i}]);
    
    % optionally transform theta such that values > 90 are reflected about
    % 90deg. 
    if bitReflectTheta
        % when reflecting theta, the distribution tends to pile up at <90
        % deg. This is difficult to fit. And so instead we invert the
        % distribution as theta' = 90 - theta' with "offset" = 90 and
        % "scale" = -1. The resulting distribution is >0 with most values
        % around 0.
        foo = 90 - (foo - 2*max(foo-90,0));
        fooS = 90 - (fooS - 2*max(fooS-90,0));
        
        %%% SEE NOTE BELOW REGARDING POOR FITS
        
        % The distribution of reflected surrogate pairwise angles is difficult to fit
        % with a model-based distribution. We've tried gamma and
        % half-normal fits. In general, the gamma fits have lower mean and
        % median error, as well as lower median difference of errors
        % (difference calculated for each "pixel" in pValue heatmap). For
        % monkey K, gamma also has the smaller maximum error. For monkey N,
        % gamma had larger maximum error for Benefit and Choice, but
        % smaller for ER. All said, if we use model-based fits, we should
        % use Gamma.
        % One can compare the errors by saving the distFit cell array as
        % distFit_hg (half gaussian) and distFit_gamma, and running the
        % following for each level of i:
        % i=3; disp(minmax(distFit_hg(i).err(:)')); disp(minmax(distFit_gamma(i).err(:)')); disp(mean(distFit_hg(i).err(:)')); disp(mean(distFit_gamma(i).err(:)')); disp(median(distFit_hg(i).err(:)')); disp(median(distFit_gamma(i).err(:)')); disp(median(distFit_hg(i).err(:) - distFit_gamma(i).err(:)));

        dist2Use = 'gamma';        
%         dist2Use = 'halfgauss';        
        % Note that this will return parameter that describe a normal PDF,
        % that must then be reflected (i.e., multiplied by 2)
        
        % tail depends on type of analysis. When reflecting theta, we are
        % interested in the absoluate similarity of angles, regardless of
        % sign, which is conveyed as small values of theta. Thus we are
        % concerned with the lower tail of the null distribution, where p
        % is the proportion of null angles LESS than our observed angle.
        % However, because we inverted the distribution above, we become
        % interested in values of theta' >> 0, thus the upper tail of the
        % distribution
        tail = 'upper';
    else
        %%% SEE NOTE BELOW REGARDING POOR FITS

        % when using UNreflected theta (range [0 180]), the distribution
        % tends to look normal and centered around 90 deg. No need to
        % transform this distribution.
%         foo = foo;
%         fooS = fooS;
        dist2Use = 'gauss';
        
        % when using UNreflected theta (range [0 180]), we are
        % interested in sign reversals of representation, that is, when
        % theta is GREATER than chance. Thus we are concered with the upper
        % tail of the null distribtion, where p is the proportion of null
        % angles GREATER than our observed angle. Because we are not
        % transforming this distribution, we remain interested in the upper
        % tail
        tail = 'upper';
    end
    
%     dist2Use = 'gamma';
%     tail = 'lower';
 
    % use empirical distribution if requested
    if strcmp('Pe',pVal2Use4Pair)
        dist2Use = 'empirical';
    else
        %%% Model-based fits to surrogate angles have been very poor. We do
        %%% not recommend using model-based p-values and instead recommend
        %%% empirical values.
        warning('Model-based fits to surrogate angles are frequently inaccurate and not recommended.')
    end
    
    % initialize pVal
    pVal{i} = NaN(size(foo));
    distFit(i).param = NaN(size(foo,1),size(foo,2),2);
    distFit(i).type = dist2Use;
    distFit(i).tail = tail;
    distFit(i).freq = NaN(size(foo,1),size(foo,2),50);
    if bitReflectTheta
        distFit(i).binC = linspace(0,90,50);
    else
        distFit(i).binC = linspace(0,180,50);
    end

    % loop through angles
    for j = 1:size(foo,1)
        if bitAuto
            % all dRAs except jth dRA
            kset = setdiff(1:size(foo,2),j);
        else
            % all dRAs
            kset = 1:size(foo,2);
        end
        for k = kset
%             tic
            [pVal{i}(j,k),distFit(i).param(j,k,:),distFit(i).err(j,k)] = ...
                sigTest(foo(j,k),squeeze(fooS(j,k,:)),tail,dist2Use,...
                'errThresh',errThreshSigTestPairwise);
%             toc
            [distFit(i).freq(j,k,:), distFit(i).binC] = hist(squeeze(fooS(j,k,:)),distFit(i).binC);
        end
    end
    
    % extract p-values of interest -- DEFUNCT. See above
%     pVal{i} = Summary.(tdrTemp).angleAnalysis.(['pVal',strAnglePair{i}]);

    % set diagonal to NaN for auto angles
    if bitAuto
        pVal{i}(logical(eye(size(pVal{i},1)))) = NaN;
    end
    
    % set p = 0 to NaN
%     pVal{i}(pVal{i}==0) = NaN;
    
    clear foo fooS
    
    % update minimum if the min pvalue in the data is greater than preset
    % min. Round to nearest 0.5 of log10(P). Ignore p = 0 entries as this
    % causes log10(p) = Inf, which cannot be plotted.
    bitZero = pVal{i}(:) == 0;
    foo = [log10(minPValStabPair):0.5:0];
    fooPos = find(foo <= log10(min(pVal{i}(~bitZero))),1,'last');
    if isempty(fooPos)
        fooPos = 1;
    end
    minPVal = min(minPVal,10.^foo(fooPos));
    
end

%%%%%%%
% reloop through pairwise angle PROBABILITIES for plotting
distFitErrMinMax = [];
for i = 1:length(condLabel);
    axN = axN + 1;
    subplot(nRow,nCol,axN);

    % plot heatmap
    h = imagesc(t,t,log10(pVal{i}),log10([minPVal maxPValPair]));
    
    % set colormap -- NOTE: In MATLAB 2019a, this must be done after
    % imagesc(). See above note.
    colormap(gca,pColormap);
    % invert colormap for black background
    if strcmp(bgColor,'k')
        foo = colormap(gca);
        colormap(gca,flipud(foo));
        clear foo
    end    

    % set NaNs to gray
    set(h,'AlphaData',~isnan(pVal{i}));
    set(gca,'Color',maskColor)
    
    axis square
    set(gca,'FontSize',14);
    
    h = colorbar('FontSize',14);
    if i ~= nCol
        set(h,'Visible','off')
    else
        % set Ticks in log scale, then relabel with linear scale
%         set(h,'XTick',0);
        set(h,'YTick',unique([log10(minPVal):1:log10(maxPValPair),log10(maxPValPair)]));
        goo = get(h,'YTick');
        clear soo
        for j = 1:length(goo)
            soo{j} = sprintf('10^{%0.2g}',goo(j));
        end
        soo{end} = sprintf('%g',maxPValPair);
%         zoo = get(h,'XTickLabel');
        set(h,'YTickLabel',soo);
%         format_ticks(h,[''],soo);
    end    
    
    % save distribution fit error min and max
    distFitErrMinMax(end+1,:) = [min(distFit(i).err(:)) max(distFit(i).err(:))];
end

%%% PLOT DISTRIBUTION FIT ERROR FOR PAIRWISE ANGLES
if bitPlotDistFitErrPair && strcmp('P',pVal2Use4Pair)
    for i = 1:length(condLabel);
        
        axN = axN + 1;
        subplot(nRow,nCol,axN);
        
        % Plot heatmap
%         h = imagesc(t,t,distFit(i).err);
        h = imagesc(t,t,(distFit(i).err),([0 max(distFitErrMinMax(:,2))]));
        
        % Set colormap -- NOTE: In MATLAB 2019a, this must be done after
        % imagesc(). See above note.
        colormap(gca,pColormap);
        % invert colormap for black background
        if strcmp(bgColor,'k')
            foo = colormap(gca);
            colormap(gca,flipud(foo));
            clear foo
        end        
        
        % set NaNs to gray
        set(h,'AlphaData',~isnan(pVal{i}));
        set(gca,'Color',maskColor)
        
        axis square
        set(gca,'FontSize',14);
        
        h = colorbar('FontSize',14);
        if i ~= nCol
            set(h,'Visible','off')
        else
            % set Ticks in log scale, then relabel with linear scale
            %         set(h,'XTick',0);
%             set(h,'YTick',[log10(minPVal):5:log10(0.05),log10(0.05)]);
%             goo = get(h,'YTick');
%             clear soo
%             for j = 1:length(goo)
%                 soo{j} = sprintf('10^{%0.2g}',goo(j));
%             end
%             soo{end} = '0.05';
%             %         zoo = get(h,'XTickLabel');
%             set(h,'YTickLabel',soo);
%             %         format_ticks(h,[''],soo);
        end
        
        if i == 1
            ylabel('Fit error pairs','FontSize',16);
        end
        
    end
end

% clear pVal


%%% PLOT distribution fits used to compute p-values if requested
if bitPlotDistFit
    figure('Name','Fits to theta null distributions')
    
    for i = 1:length(condLabel);
        subplot(1,length(condLabel),i);
        hold on;
        title(condLabel{i});
        
        binWidth = mean(diff(distFit(i).binC));
        
        % loop through all rows and columns above diagonal
        for j = 1:size(distFit(i).param,1)-1
            for k = j+1:size(distFit(i).param,2)
                
                % plot rand fraction
                if rand < 0.9
                    continue
                end
                
                % plot histograms
                freq = squeeze(bsxfun(@rdivide,distFit(i).freq(j,k,:),sum(distFit(i).freq(j,k,:)))/binWidth);
                plot(distFit(i).binC,freq,'Color',[0.6 0.6 0.6]);
                
                % set PDF params
                if strcmp(distFit(i).type,'halfgauss')
                    % use normal distribution and scale by 2. Note that
                    % even when fit params were obtained using
                    % fitdist('half normal'), the params are the same and
                    % thus the below method works either way. If we wanted
                    % to be fancy, we could pass 'half normal' to PDF(),
                    % but this would require scaling by 1 and thus
                    % different operation depending on whether
                    % fitdist('half normal') was used.
                    foo.dt = 'normal';
                    foo.scale = 2;
                elseif strcmp(distFit(i).type,'gauss')
                    % use normal distribution and scale by 1
                    foo.dt = 'normal';
                    foo.scale = 1;                    
                else
                    % use specified distribution and scale by 1
                    foo.dt = distFit(i).type;
                    foo.scale = 1;                    
                end
   
                % plot fits
                if ~any(isnan(distFit(i).param(j,k,:)))
                    plot(distFit(i).binC,foo.scale * pdf(foo.dt,distFit(i).binC,distFit(i).param(j,k,1),distFit(i).param(j,k,2)),'Color',lineColor);
                end
                
%                 if any(freq > 0.08)
%                     error('stop')
%                 end
                clear foo freq
            end
        end
    end
end

% annotation
[figTitleHand] = annotateFig(figName,mfilename);

%% PLOT 'TOPO MAP' OF ANGLES

% set which class of angle fit
% 'decay' for exponential decay fits
% 'step' for step function fits
% 'step_trans' for step function TRANS fit (angles that are anti-similar)
fitClass2u = 'step'; 

% logical on whether to plot fits from 2-step fit functions
bit2Step = 0;

% logical whether to horizontally align each row of angles such that the
% time of the reference dRA is time = 0 (TRUE). When FALSE, horizontal
% position refers to the absolute time of the comparison dRA (and time = 0
% refers to the task event to which the original PTSHs were aligned).
bitAlignRows = 1; 

% set vertical offset (in degrees) between rows
offset = -30;

% logical on whether to plot auto-angle (within signal) or cross-angles
% (between signal)
bitAuto = 1;

% logical on whether to plot error shading as simple lines for output to
% illustrator (applies to average profile only).
bit4Illustrator = false;

% PLOT
stabilityPlot(Summary.angleProfile,'bitAlignRows',bitAlignRows,'offset',offset,...
    'fitClass2u',fitClass2u,'bit2Step',bit2Step,'bitAuto',bitAuto,...
    'bit4Illustrator',bit4Illustrator); 

%% PLOT ANGLE PROFILE PARAMETERS 

% Plots two figures showing the angle profile fit parameters for each pair
% of task variables. Generates TWO or THREE figures:
% (1) Fit parameter as function of time of the reference dRA
% (2) Scatter plots between all pairs of fit parameters.
% (3) (Optional) When angle profiles are fit to each surrogate dataset
% (i.e., field Summary.angleProfile.([analysis_name]).bS), distributions of
% surrogate fit parameters vs. veridical fit parameter are shown. 

% NOTE: these figures are not included in Kimmel et al., 2020.

% set which type of angle fit. See fieldnames(Summary.angleProfile) for all
% options.
% 'main_decay' for exponential decay fits
% 'main_step' for step function fits
% 'main_step_2' for 2-step function fits
% 'main_step_trans' for step function TRANS fit (angles that are anti-similar)
s2u = 'main_step'; 
f2u = 'fitO_step';

zThresh = 5; % zscore threshold above which to eliminate outliers when plotting histograms

alpha = 0.05; % significance threshold

pVal2Use = 'bPe'; %'bPe' for empirical p-value; 'bP' for model-based CDF

fontSize = 16;

bitPlotMinFromIter = 1; % whether to add markers for the minimum GOF across all iterations (i.e., make sure that the fit routine converges).

%%%%%%%%%%%
f1 = figure('Name','Angle profile parameters');
f2 = figure('Name','Angle profile parameter scatters');
if isfield(Summary.angleProfile(1).(s2u),'bS') && ~isempty(Summary.angleProfile(1).(s2u).bS)
    f3 = figure('Name','Surrogate Angle profile parameters');
end

t = Summary.angleProfile(1).(s2u).tMaster;
nSignal = length(Summary.angleProfile);

for i = 1:nSignal
    
%     switch i
%         case 1
%             strTitle = 'Benefit';
%         case 2
%             strTitle = 'Choice';
%         case 3
%             strTitle = 'Expected Reward';
%     end
    strTitle = Summary.angleProfile(i).(s2u).signalDesc;
    
    % determine type of fit
    switch Summary.angleProfile(i).(s2u).fitTypeName
        case 'exp1'
%             nParam = 2;
            paramName = {'A_1','\lambda_1'};
            param2Combine = {[1],[2]};
            param2Sum = {[],[]};
            limitCI = [-10 -100;
                       100 10];           % set lower and upper (rows) limits
                                          % at which to include confidence
                                          % intervals for each parameter
                                          % (columns)
        case 'exp2'
%             nParam = 4;
            paramName = {'A_1','\lambda_1','A_2','\lambda_2'};
            param2Combine = {[1 3],[2],[4]};
            param2Sum = {[1 3],[],[]};
            limitCI = [-10 -100 -10 -100;
                       100 10   100  10]; % set lower and upper (rows) limits
                                          % at which to include confidence
                                          % intervals for each parameter
                                          % (columns)
        case 'exp2_df3'
            paramName = {'A_1','\lambda_1','\lambda_2'};
            param2Combine = {[1],[2],[3]};
            param2Sum = {[],[],[]};
            limitCI = [-10 -100 -100;
                       100 10   10]; % set lower and upper (rows) limits
                                          % at which to include confidence
                                          % intervals for each parameter
                                          % (columns)
        case 'exp2_df2'
            paramName = {'\lambda_1','\lambda_2'};
            param2Combine = {[1],[2]};
            param2Sum = {[],[]};
            limitCI = [ Inf Inf;
                        -Inf   -Inf]; % set lower and upper (rows) limits
                                          % at which to include confidence
                                          % intervals for each parameter
                                          % (columns)
        case 'SquarePulse'
            paramName = Summary.angleProfile(i).(s2u).coeffNames; % {'A','Start','Width'};
            param2Combine = {[1],[2],[3]};
            param2Sum = {[],[],[]};
            ttemp = Summary.angleProfile(i).(s2u).tMaster;
            dttemp = ttemp(2) - ttemp(1);
            limitCI = [-5 min(ttemp)-dttemp 0;
                       95 max(ttemp)+dttemp range(ttemp)+dttemp]; 
                                          % set lower and upper (rows) limits
                                          % at which to include confidence
                                          % intervals for each parameter
                                          % (columns)
            clear ttemp dttemp
            
        otherwise
            error('Fit type %s not recognized',Summary.angleProfile(i).(s2u).fitTypeName);
    end    
    

    % extract confidence intervals
    if isfield(Summary.angleProfile(i).(s2u),'bint') && ~all(isnan(Summary.angleProfile(i).(s2u).bint(:)))
        fooL = squeeze(Summary.angleProfile(i).(s2u).bint(:,:,1));
        fooU = squeeze(Summary.angleProfile(i).(s2u).bint(:,:,2));
        
        % exclude confidence intervals out of range
        bitIncludeCI_L = fooL >= repmat(limitCI(1,:),size(Summary.angleProfile(i).(s2u).bint,1),1);
        fooL(~bitIncludeCI_L) = NaN;
        bitIncludeCI_U = fooU <= repmat(limitCI(2,:),size(Summary.angleProfile(i).(s2u).bint,1),1);
        fooU(~bitIncludeCI_U) = NaN;
        
        % identify confidence intervals that do not include zero
        bitCISig = fooL .* fooU > 0 & fooU~=0 & fooL~=0;
        
    else
        fooL = [];
        fooU = [];
        bitCISig = [];
    end
    
    nRow = length(param2Combine) + 1; % +1 for GOF
    nParam = size(Summary.angleProfile(i).(s2u).b,2);
    
    % determine if 1:1 parameter:axes
    bit1ParamPerAxes = true;
    for j = 1:length(param2Combine)
        if length(param2Combine{j})>1
            bit1ParamPerAxes = false;
        end
    end
    for j = 1:length(param2Sum)
        if ~isempty(param2Sum{j})
            bit1ParamPerAxes = false;
        end
    end
    
    % time step size:
    dt = Summary.angleProfile(i).(s2u).tMaster(2) - Summary.angleProfile(i).(s2u).tMaster(1);

    for rowN = 1:nRow-1
        
        % set 1 or 2-tailed confidence intervals
        if isfield('lambdaPN',Summary.angleProfile(i).(s2u)) && ...
                bit1ParamPerAxes && ...
                ismember(param2Combine{rowN},Summary.angleProfile(i).(s2u).lambdaPN)
            % 1-tailed
            % note that if one-tailed, there's no reason to plot anything
            % above the median
            prctile2Use = 100*[alpha 0.5]; 
            
        else
            % 2-tailed
            prctile2Use = 100*[alpha/2 1-alpha/2];
        end            
        
        %%% FIGURE 1 -- real data params
        figure(f1);
        subplot(nRow,nSignal, (rowN-1)*nSignal + i,'FontSize',fontSize);
        % note that with certain models, there may be additional dimensions
        % of b
        if size(Summary.angleProfile(i).(s2u).b(:,param2Combine{rowN},:),2) > 1 && size(Summary.angleProfile(i).(s2u).b(:,param2Combine{rowN},:),3) > 1
            error('Cannot handle when B parameters have non-singleton dimensions in the 2nd and 3rd orders of the tensor')
        end
        if ~isempty(fooL)
            errorbar(squeeze(repmat(t',1,length(param2Combine{rowN}),size(Summary.angleProfile(i).(s2u).b,3))),squeeze(Summary.angleProfile(i).(s2u).b(:,param2Combine{rowN},:)),fooL(:,param2Combine{rowN}),fooU(:,param2Combine{rowN}));
        else
            plot(squeeze(repmat(t',1,length(param2Combine{rowN}),size(Summary.angleProfile(i).(s2u).b,3))),squeeze(Summary.angleProfile(i).(s2u).b(:,param2Combine{rowN},:)))
        end
        hold on;
        legText = paramName(param2Combine{rowN});
        % add a summed row if selected
        if ~isempty(param2Sum{rowN})
            plot(t',sum(Summary.angleProfile(i).(s2u).b(:,param2Sum{rowN}),2),'k');
            legText{end+1} = 'Sum';
        end
        
        % legend vs. surrogate
        if ~bit1ParamPerAxes
            legend(legText,'FontSize',fontSize,'box','off');
        else
            % include surrogate data
            if isfield(Summary.angleProfile(i).(s2u),'bS')
                y = median(Summary.angleProfile(i).(s2u).bS(:,param2Combine{rowN},:),3);
                yErrL = prctile(Summary.angleProfile(i).(s2u).bS(:,param2Combine{rowN},:),prctile2Use(1),3);
                yErrU = prctile(Summary.angleProfile(i).(s2u).bS(:,param2Combine{rowN},:),prctile2Use(2),3);
                errorbar(t',y,yErrL-y,yErrU-y,'Color',[0.6 0.6 0.6]);
                clear y yErrU yErrL
            end
        end
        
        % mark significant points
        if ~isempty(bitCISig)
            goo = Summary.angleProfile(i).(s2u).b;
            goo(~bitCISig) = NaN;
            plot(repmat(t',1,length(param2Combine{rowN})),...
                goo(:,param2Combine{rowN}),'x');
            clear goo
        end
        
        % if 2-step plot, add 1-step params
        if isfield(Summary.angleProfile(i).(s2u),'fitMethod') && ...
            strcmp(Summary.angleProfile(i).(s2u).fitMethod,'paramSearch') && ...
            strcmp('_2',s2u(end-1:end)) && isfield(Summary.angleProfile(i),s2u(1:end-2))
            
            plot(squeeze(repmat(t',1,length(param2Combine{rowN}),size(Summary.angleProfile(i).(s2u(1:end-2)).b,3))),...
                squeeze(Summary.angleProfile(i).(s2u(1:end-2)).b(:,param2Combine{rowN},:)),'k--');
            
            legend({'2-step, 1','2-step, 2','1-step'},'box','off');
        end
        
        % labeling
        if rowN == 1
            title(strTitle,'FontSize',fontSize-2);
        end
        if i == 1;
            if bit1ParamPerAxes
                ylabel(paramName{rowN},'FontSize',fontSize);
            else
                ylabel('Parameter value','FontSize',fontSize);
            end
        end
        if rowN == nRow
            xlabel('Time from offer','FontSize',fontSize);
        end
        set(gca,'FontSize',fontSize)
        xlim([Summary.angleProfile(i).(s2u).tMaster(1) Summary.angleProfile(i).(s2u).tMaster(end) + dt]);
    end

    % add GOF plot
    subplot(nRow,nSignal, (rowN)*nSignal + i,'FontSize',fontSize);
    hold on;
    
    
    if isfield(Summary.angleProfile(i).(s2u),'fitMethod') && ...
            strcmp(Summary.angleProfile(i).(s2u).fitMethod,'paramSearch')
        plot(t,Summary.angleProfile(i).(s2u).gof,'k')
        % if 2-step, look for 1 step run and plot error
        if strcmp('_2',s2u(end-1:end)) && isfield(Summary.angleProfile(i),s2u(1:end-2))
            plot(t,Summary.angleProfile(i).(s2u(1:end-2)).gof,'k--');
            legend({'2-step','1-step'},'box','off');
        end
            
    else
        % extract GOF value
        gof(i).sse = NaN(length(Summary.angleProfile(i).(f2u)),1);
        gof(i).sse_max = NaN(length(Summary.angleProfile(i).(f2u)),1);
        gof(i).adjrsquare = NaN(length(Summary.angleProfile(i).(f2u)),1);
        gof(i).rsquare = NaN(length(Summary.angleProfile(i).(f2u)),1);
        for g = 1:length(Summary.angleProfile(i).(f2u))
            gof(i).sse(g) = Summary.angleProfile(i).(f2u)(g).gof.sse;
            gof(i).sse_max(g) = Summary.angleProfile(i).(f2u)(g).gof.sse_max;
            gof(i).adjrsquare(g) = Summary.angleProfile(i).(f2u)(g).gof.adjrsquare;
            gof(i).rsquare(g) = Summary.angleProfile(i).(f2u)(g).gof.rsquare;
        end
        
        plot(t,gof(i).sse,'k')
    %     plot(t,gof(i).adjrsquare,'k')

    %     plot(t,Summary.angleProfile(i).(s2u).gof,'k')

        % optionally include the minimum across all iterations 
        if bitPlotMinFromIter
    %         y = min(Summary.angleProfile(i).(s2u).gof_min,[],2);
            y = gof(i).sse_max;
            plot(t',y,'v','Color','k');
            clear y
        end
    end

    % include surrogate data
    if isfield(Summary.angleProfile(i).(s2u),'gofS')
        y = median(Summary.angleProfile(i).(s2u).gofS,2);
        yErrL = prctile(Summary.angleProfile(i).(s2u).gofS,100*alpha/2,2);
        yErrU = prctile(Summary.angleProfile(i).(s2u).gofS,100*(1-alpha/2),2);
        errorbar(t',y,yErrL-y,yErrU-y,'Color',[0.6 0.6 0.6]);
        clear y yErrL yErrU
        
        % optionally include the MAXIMUM across all iterations of each
        % surrogate
        if bitPlotMinFromIter
            y = max(Summary.angleProfile(i).(s2u).gofS-Summary.angleProfile(i).(s2u).gofS_min,[],2);
            plot(t',y,'v','Color',[0.6 0.6 0.6]);
            clear y
        end
    end
    
    % appearance
    if i == 1
        ylabel('SSE','FontSize',fontSize);
%         ylabel('Adj r-square','FontSize',fontSize);
    end
%     ylim([0 1]);
    xlabel('dRA time (time from offer)','FontSize',fontSize);
    xlim([Summary.angleProfile(i).(s2u).tMaster(1) Summary.angleProfile(i).(s2u).tMaster(end) + dt]);

    %%% FIGURE 2 -- surrogate data params
    figure(f2);
    
    % number of scatters=rows
    nRow = (nParam^2-nParam)/2;
    
    % make look-up table for contents of each row
    pLUT = cell(nParam,nParam);
    for j = 1:size(pLUT,1)
        for k = j:size(pLUT,1)
            pLUT{j,k}= [j k];
        end
    end
    
    posLUT = 0;
    
    % loop through scatters (rows)
    for rowN = 1:nRow
        subplot(nRow,nSignal, (rowN-1)*nSignal + i,'FontSize',fontSize);
        hold on
        
        posLUT = posLUT+1;
        while isempty(pLUT{posLUT}) || length(unique(pLUT{posLUT}))==1
            posLUT = posLUT+1;
        end
        
        % find significant points
        x = Summary.angleProfile(i).(s2u).b(:,pLUT{posLUT}(1));
        y = Summary.angleProfile(i).(s2u).b(:,pLUT{posLUT}(2));
        if ~isempty(bitCISig)
            bitXSig = bitCISig(:,pLUT{posLUT}(1));
            bitYSig = bitCISig(:,pLUT{posLUT}(2));
            bitSig = bitXSig & bitYSig;
        else
            bitSig = false(size(x));
        end
        % datatip labels
        dtLabel = sprintf('t=%0.2f\n',Summary.angleProfile(i).(s2u).tMaster);
        % make cell array
        dtLabel = regexp(dtLabel,'\n','split');
        % eliminate last element, which is empty
        dtLabel = dtLabel(1:end-1)';
        % plot Scatter
        plotScatterN(x,y,'bitSig',bitSig,'statTest','corrPear',...
            'bitEqualAxes',0,'alpha',alpha,'datatipLabel',dtLabel,...
            'pTextFontSize',12);
        axis square

        % labeling
        if rowN == 1
            title(strTitle,'FontSize',fontSize-2);
        end
        xlabel(paramName{pLUT{posLUT}(1)},'FontSize',fontSize);
        ylabel(paramName{pLUT{posLUT}(2)},'FontSize',fontSize);
    
        
    end
    
    % only do if there are surrogate data
    if isfield(Summary.angleProfile(i).(s2u),'bS') && ~isempty(Summary.angleProfile(i).(s2u).bS)
        for pN = 1:nParam
            %%% FIGURE 3 -- surrogate data params
            figure(f3);
            subplot(nParam,nSignal, (pN-1)*nSignal + i,'FontSize',fontSize);
            hold on
            
            % extract data
            d = squeeze(Summary.angleProfile(i).(s2u).bS(:,pN,:))';
            
            % scale and shift
            if length(unique(Summary.angleProfile(i).(s2u).bS_gammaFit_sign(:,pN))) > 1
                error('Sign flip must be constant across dRAs');
            end
            signFlip = Summary.angleProfile(i).(s2u).bS_gammaFit_sign(1,pN);
            
            % converts data into format in which it was originally fit
            d = signFlip*d + (repmat(Summary.angleProfile(i).(s2u).bS_gammaFit_offset(:,pN)',size(d,1),1));
            
            % eliminate outliers
            dz = zscore(d);
            d(abs(dz) > zThresh) = NaN;
            
%             if pN == 1
%                 d = repmat(-max(d,[],1),size(d,1),1)+d;
%             end
            
            % histogram of parameter values across surrogates
            [his,binC] = hist(d,30);
            
            % prepare PDF fits
            % loop through time of dRAs
            for tN = 1:size(Summary.angleProfile(i).(s2u).bS_gammaFit,1)
                
%                 if pN == 1
%                     signFlip = -1;
%                 else
%                     signFlip = 1;
%                 end
                
                y = pdf(Summary.angleProfile(i).(s2u).bS_gammaFit(tN,pN),binC);
                
                % normalized histogram, such that the sum under
                % the histogram = sum under PDF in the current range.
                foo = sum(y)*diff(binC(1:2));
                hisTemp = his(:,tN)/(sum(his(:,tN),1)/foo*diff(binC(1:2)));
                %             his = his./repmat(sum(his,1)/foo*diff(binC(1:2)),size(his,1),1);
                clear foo
                
                % determine the x-coordinate to plot based on the data
                % transformation
                x2Plot = signFlip * (binC - Summary.angleProfile(i).(s2u).bS_gammaFit_offset(tN,pN));
                
                % plot the epirical histogram
                plot(x2Plot,hisTemp,'Color',[0.6 0.6 0.6]);
                %         plot(binC,his','Color',[0.6 0.6 0.6]);
                
                % plot fit
                plot(x2Plot,y,'r');
            end
            
            %%%%
            % PLOT DATA (non-surrogate) values
            yLim = get(gca,'YLim');
            
            % plot real data values vertically with earliest dRA at the top
            b = Summary.angleProfile(i).(s2u).b(:,pN);
            plot(flipud(b),...
                linspace(yLim(1),yLim(2),size(b,1)),'k*');
            
            % overlay those values that are significant
            b(Summary.angleProfile(i).(s2u).(pVal2Use)(:,pN)>alpha) = NaN;
            plot(flipud(b),...
                linspace(yLim(1),yLim(2),size(b,1)),'b*');
            
            %%%%
            
            % labeling
            if pN == 1
                title(strTitle,'FontSize',fontSize-2);
            end
            if i == 1;
                ylabel(sprintf('Norm. Frequency (%s)',paramName{pN}),'FontSize',fontSize);
            end
            if pN == nParam
                xlabel('Parameter value','FontSize',fontSize);
            end
        end
    end
end


%% Summarize Angle Profile Parameters...
% ...for analysis in which we extract angles from surrogate datasets during
% period defined by angle profile from veridical data.
%
% Used to generate Supp Fig 21 in Kimmel et al., 2020.

% set which type of angle fit -- only works for step function fits
% 'main_step' for step function CIS fit (angles that are similar)
% 'main_step_trans' for step function TRANS fit (angles that are anti-similar)
% See fieldnames(Summary.angleProfile) for all options.
s2u = 'main_step'; 

pVal2Use = 'Pe'; %'Pe' for empirical p-value; 'P' for model-based CDF

% Limit number of signal pairs to plot. For instance, the auto correlograms
% are signals 1:3. Leave empty to plot all signal pairs
sig2Plot = [1:3];

% alpha level for significance
% switch animalName
%     case 'Norris'
%         alpha = 10^-2.5;
%         alpha = 10^-2;
%     case 'Kirby'
%         alpha = 10^-3;
%     otherwise
%         error('Animal name %s not recognized',animalName);
% end
% set min Pval for inf points based on value in "ANGLES" section above
minPValStab = 10^-4; % will cause error if value not assigned from above
% override alpha threshold with the maxPVal from "ANGLES" section above. 
alpha = maxPValStab; % will cause error if value not assigned from above
 
% whether to correct for multiple comparisons across times (dRAs)
bitMultComp = false;

fontSize = 16;
markerSize = 7;
surrColor = [0.6 0.6 0.6];

% signalName = {'Benefit','Choice','Expected Reward'};
paramName = {'A','Start','Width'};

%%%%%%%%%%%
% generate signals to plot
if isempty(sig2Plot)
    sig2Plot = 1:length(Summary.angleProfile);
end

nSignal = length(sig2Plot);
nRow = 5;
nCol = nSignal;
axH1 = NaN(1,nSignal);
axH2 = NaN(1,nSignal);
axH3 = NaN(1,nSignal);

%%%%%%%%%%%
figName = 'Angle step-function analysis';
figure('Name',figName);
whitebg(gcf,bgColor);
set(gcf,'Color',bgColor);

% yLim1 = NaN(nSignal,2);
count = 0;
% loop through signals
for i = sig2Plot
    count = count + 1;
    
    % corrected alpha
    if bitMultComp
        % bonforonni for now
        alphaC = alpha / sum(isfinite(Summary.angleProfile(i).(s2u).step_thetaMean));
    else
        alphaC = alpha;
    end
    
    % extract
    t = Summary.angleProfile(i).(s2u).tMaster;
    dt = t(2) - t(1);
    % number of steps
    nSteps = size(Summary.angleProfile(i).(s2u).b,3);
    
    % decide on type of percentile to use (empirical from raw surrogates or
    % as measured from fit to surrogate distribution)
    switch pVal2Use
        case 'P'
            CI2Use = 'step_thetaSDist';
        case 'Pe'
            CI2Use = 'step_thetaS';
    end
    % decide on level of percentile to use depending on alpha:
    CISet = [0.05 0.01 0.001];
    [pos] = find(CISet >= alpha,1,'last');
    switch pos
        case 1
            CI2Use = [CI2Use,'95'];
        case 2
            CI2Use = [CI2Use,'99'];
        case 3
            CI2Use = [CI2Use,'999'];
    end
    %%%% Plot mean theta under step
    subplot(nRow,nCol,count);
    hold on;
    bitSig = Summary.angleProfile(i).(s2u).(['step_theta',pVal2Use]) <= alphaC;
    % generate time matrix for plotting
    t4p = repmat(t',1,nSteps); 
    t4pSig = t4p;
    t4pNoSig = t4p;
    t4pSig(~bitSig) = NaN;
    t4pNoSig(bitSig) = NaN;
    
    % surrogate data mean/median and 99.9th prctile
    h1s = plot(t4p,Summary.angleProfile(i).(s2u).step_thetaSMedian,'o','MarkerSize',markerSize-2);
    h2s = plot(t4p,Summary.angleProfile(i).(s2u).(CI2Use),'v','MarkerSize',markerSize-3);
    % set to surround color for single steps
    if nSteps == 1
        set([h1s h2s],'Color',surrColor);
    else
        for j = 1:length(h1s)
            set(h2s(j),'Color',get(h1s(j),'Color'));
        end
%         alpha([h1 h2],0.5);
        % simply lighten color
%         cfoo = get(h1,'Color');
%         for j = 1:length(cfoo)
%             set(h1(j),'Color',min(cfoo{j}*1.2,1));
%         end
%         
%         cfoo = get(h2,'Color');
%         for j = 1:length(cfoo)
%             set(h2(j),'Color',min(cfoo{j}*1.2,1));
%         end
    end
    
    % all points
    h1 = plot(t4p,Summary.angleProfile(i).(s2u).step_thetaMean,'o-','MarkerSize',markerSize,'LineWidth',1);
    % significant points
    h2 = plot(t4pSig,Summary.angleProfile(i).(s2u).step_thetaMean,'o','MarkerSize',markerSize);
    % set color for single steps
    if nSteps == 1
        set([h1 h2],'Color','k');
        set(h2,'MarkerFaceColor','k');
    else
        for j = 1:length(h1)
            set(h1(j),'Color',get(h1s(j),'Color'));
            set(h2(j),'Color',get(h2s(j),'Color'));
            set(h2(j),'MarkerFaceColor',get(h2(j),'Color'));
        end
    end
    
    % appearance
    set(gca,'FontSize',fontSize,'box','off');
    xlim([t(1)-dt/2 t(end)+dt/2]);
    if count == 1
        ylabel('Stability magnitude (deg)','FontSize',fontSize+4)
    end
    title(Summary.angleProfile(i).(s2u).signalDesc,'FontSize',fontSize);
    foo = minmax(Summary.angleProfile(i).(s2u).step_thetaSMean');
    goo = minmax(Summary.angleProfile(i).(s2u).step_thetaMean');
    soo = [min(foo(1),goo(1)), max(foo(2),goo(2))];
    moo = [-180:10:180];
    pos = [find(moo<soo(1),1,'last'), find(moo>soo(2),1,'first')];
%     yLim1(count,:) = moo(pos);
    clear foo goo soo moo pos
    
    %%%% Plot duration of step
    axH2(count) = subplot(nRow,nCol,count+nCol);
    hold on;
    foo = squeeze(Summary.angleProfile(i).(s2u).b(:,3,:));
    % points (based on magnitude of theta above)
    h1 = plot(t4p,foo,'o-','MarkerSize',markerSize);
    % significant points (based on magnitude of theta above)
    h2 = plot(t4pSig,foo,'o','MarkerSize',markerSize);
    
    % set color for single steps
    if nSteps == 1
        set([h1 h2],'Color','k');
        set(h2,'MarkerFaceColor','k');
    else
        for j = 1:length(h2)
            set([h1(j) h2(j)],'Color',get(h1s(j),'Color'))
            set(h2(j),'MarkerFaceColor',get(h1s(j),'Color'));
        end
    end
    
    % appearance
    set(gca,'FontSize',fontSize,'box','off');
    xlim([t(1)-dt/2 t(end)+dt/2]);
    ylim([0 range([t(1)-dt/2 t(end)+dt/2])]);
    if count == 1
        ylabel('Stability duration (s)','FontSize',fontSize+4)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Plot p-value of amplitude parameter vs. Time
    axH3(count) = subplot(nRow,nCol,count+2*nCol);
    hold on;
    pValTemp = Summary.angleProfile(i).(s2u).(['step_theta',pVal2Use]);
    % find p=0
    bitPZero = pValTemp==0;
    % replace p=0 values with p-val min
    pValTemp(bitPZero) = minPValStab;
    pValTemp = log10(pValTemp);
    % plot all p-values
    plot(t4p,pValTemp,...
        '-','MarkerSize',markerSize,'Color',lineColor);
    % add markers for p>0, non-signficant
    h1 = plot(t4pNoSig(~bitPZero),pValTemp(~bitPZero),...
        'o','MarkerSize',markerSize,'Color',lineColor);
    % add markers for p=0, non-significant
    h2 = plot(t4pNoSig(bitPZero),pValTemp(bitPZero),...
        'v','MarkerSize',markerSize,'Color',lineColor);

    % add markers for p>0, significant
    h1S = plot(t4pSig(~bitPZero),pValTemp(~bitPZero),...
        'o','MarkerSize',markerSize,'Color',lineColor);
    % add markers for p=0, significant
    h2S = plot(t4pSig(bitPZero),pValTemp(bitPZero),...
        'v','MarkerSize',markerSize,'Color',lineColor);

    % set color for single steps
    if nSteps == 1
        set([h1S h2S],'MarkerFaceColor',lineColor);
    else
        for j = 1:length(h2)
            error('We currently do not support summary figure for 2-step fits.')
%             set([h1(j) h2(j)],'Color',get(h1s(j),'Color'))
%             set(h2(j),'MarkerFaceColor',get(h1s(j),'Color'));
        end
    end
    
    % add alpha
    xlim([t(1)-dt/2 t(end)+dt/2]);    
    xLim = get(gca,'XLim');
    line(xLim,log10([alphaC alphaC]),'LineStyle','--','Color',surrColor);
    
    % appearance
    set(gca,'FontSize',fontSize,'box','off');
    if count == 1
        ylabel('log_{10}(P(magnitude))','FontSize',fontSize+4)
    end
    if count == ceil(nSignal/2)
        xlabel('Time from offer (s)','FontSize',fontSize+4)
    end
    % max y-axis should be log10(p=1) --> 0
    yLim = get(gca,'YLim');
    ylim([yLim(1) 0]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%% Plot distribtion of p-value of amplitude parameter
    axH4(count) = subplot(nRow,nCol,count+3*nCol);
    hist(log10(Summary.angleProfile(i).(s2u).(['step_theta',pVal2Use])),10);
    
    % get handle to histogram
    hi = findobj(gca,'-property','FaceColor');
    % Note that findobj retrieves handles in opposite order, so need to
    % reverse order:
    hi = flipud(hi);
    % match face color to line colors above. 
    for j = 1:length(hi)
        set(hi(j),'FaceColor',get(h1s(j),'Color'),...
            'EdgeColor',get(h1s(j),'Color'));
    end
    
    hold on;
    % add alpha
    yLim = get(gca,'YLim');
    line(log10([alphaC alphaC]),yLim,'LineStyle','--');
    
    % appearance
    set(gca,'FontSize',fontSize,'box','off');
    if count == ceil(nCol/2) 
        xlabel('log_{10}(P(step))','FontSize',fontSize+4)
    end
    if count==1
        ylabel('Count','FontSize',fontSize+4)
    end

    %%%% Plot distribution of mean theta under step across surrogates
    axH5(count) = subplot(nRow,nCol,count+4*nCol);
    hold on;
    % this whole process needs to be separate for each step:
    for st = 1:nSteps
        [his,binC] = hist(Summary.angleProfile(i).(s2u).step_thetaSDist(:,:,st)',30);
        dtBin = binC(2) - binC(1);
        for tN = 1:size(Summary.angleProfile(i).(s2u).step_thetaSDistFit,1)
            
            if iscell(Summary.angleProfile(i).(s2u).step_thetaSDistFit(tN,st)) && ...
                    ~isempty(Summary.angleProfile(i).(s2u).step_thetaSDistFit{tN,st})
                % when fit parameters are returned as a fit object
                % (by matlab's fitdist function).
                y = pdf(Summary.angleProfile(i).(s2u).step_thetaSDistFit{tN,st},binC);
                
                % normalized histogram, such that the sum under
                % the histogram = sum under PDF in the current range.
                %             foo = sum(y)*diff(binC(1:2));
            elseif ~iscell(Summary.angleProfile(i).(s2u).step_thetaSDistFit(tN,st)) && ...
                    ~any(isnan(Summary.angleProfile(i).(s2u).step_thetaSDistFit(tN,:,st)))
                % for when fit parameters are returned in a cell array
                % (by GE's sigtest function)
                y = gampdf(binC,Summary.angleProfile(i).(s2u).step_thetaSDistFit(tN,1,st),...
                    Summary.angleProfile(i).(s2u).step_thetaSDistFit(tN,2,st));
                %             foo = sum(y)*diff(binC(1:2));
            else
                y = NaN;
                foo = 1;
            end
            hisTemp = his(:,tN)/(sum(his(:,tN),1)*dtBin);
            %             his = his./repmat(sum(his,1)/foo*diff(binC(1:2)),size(his,1),1);
            clear foo
            
            % plot the epirical histogram
            plot(binC,hisTemp,'Color',surrColor);
            
            % plot fit
            h1 = plot(binC,y,'Color',get(h1s(st),'Color'));
            % set to black if 1 step
            if nSteps == 1
                set(h1,'Color','k')
            end
            
        end
        
        % PLOT DATA (non-surrogate) values
        yLim = get(gca,'YLim');
        
        % plot real data values vertically with earliest dRA at the top
        b = Summary.angleProfile(i).(s2u).step_thetaMean(:,st);
        h1 = plot(flipud(b),...
            linspace(yLim(1),yLim(2),size(b,1)),'o',...
            'Color',get(h1s(st),'Color'));
        
        % fillin those values that are significant
        b(~bitSig(:,st),st) = NaN;
        h2 = plot(flipud(b),...
            linspace(yLim(1),yLim(2),size(b,1)),'o',...
            'Color',get(h1s(st),'Color'),...
            'MarkerFaceColor',get(h1s(st),'Color'));
        
        % set to black if 1 step
        if nSteps == 1
            set([h1 h2],'Color','k');
            set(h2,'MarkerFaceColor','k');
        end
        
        clear b
    end
    
    % appearance
    set(gca,'FontSize',fontSize,'box','off');
    if count == 1
        ylabel('Count','FontSize',fontSize+4)
    end
    if count == ceil(nSignal/2)
        xlabel('Stability similarity of surrogates (deg)','FontSize',fontSize+4)
    end
   
end

% annotate figure
annotateFig(figName);

% set limits of axes that have angle magnitude
% set(axH1,'YLim',[min(yLim1(:,1)) max(yLim1(:,2))]);
% set(axH5,'XLim',[min(yLim1(:,1)) max(yLim1(:,2))]);

%% COMPUTE magnitude of dRA projected into low-D space
% TODO: this in theory could be moved to the analysis script file. It could
% also probably be generalized fairly easily. However, because it is so
% spefific and really only applies with the plotting code in the following
% block, we'll leave it here for now.

normFN = 'normdRAs';

% matrix of sRAs
% sRA = NaN(length(Summary.(otdr).sRA.sRA1),length(condLabel));
% for i = 1:length(condLabel)
%     sRA(:,i) = Summary.(otdr).sRA.(['sRA',num2str(i)]);
% end
sRA = Summary.(otdr).sRA.RA;
nSRA = size(sRA,2);

clear lowD_val lowD_null lowD_P

% loop through conditions
for i = 1:length(condLabel)
    
%     switch i
%         case 1
%             foo = [normFN,'1'];
%         case 2
%             foo = [normFN,'2'];
%         case 3
%             foo = [normFN,'3'];
%         otherwise
%             error('number of conditions exceeds max');
%     end
    
    % extract normlized dRA
    dRA = Summary.(tdr).dRAs(:,:,i)';
    [rT] = size(dRA,2); % number of regression times
    % extract dRA magnitude vs. time
    dRA_mag = Summary.(tdr).(normFN)(:,i)';
    % rescale normalized dRA to a non-normalized dRA:
    dRA_orig = bsxfun(@times,dRA,dRA_mag);
    clear dRA_mag
    
    % compute vector overlap with low D space
    lowD_val(i) = vectorOverlap(sRA,dRA_orig);

    % extract output as separate vars
    dRA_lowD_prop = lowD_val(i).dRA_lowD_prop;
    dRA_lowD_prop_bySRA = lowD_val(i).dRA_lowD_prop_bySRA;
    dRA_lowD_mag = lowD_val(i).dRA_lowD_mag;
    dRA_lowD_mag_bySRA = lowD_val(i).dRA_lowD_mag_bySRA;
    dRA_lowD_magDiff = lowD_val(i).dRA_lowD_magDiff;
    dRA_lowD_magDiff_bySRA = lowD_val(i).dRA_lowD_magDiff_bySRA;

%     % compute dRA magnitude, which should be 1 as the dRAs are normalized
%     % to unit vectors
%     dRA_mag_unit = sqrt(sum(dRA .^ 2,2));
%     if any(abs(dRA_mag_unit - 1) > eps*100)
%         error('dRA magnitude should be 1');
%     end
%     % project high-D dRA into low-D space
%     dRA_lowD = dRA * sRA;
%     
%     % find the proportion of the dRA captured by each of the 3 sRAs. We
%     % normalize this magnitude by the magnitude of the high-D dRA, even
%     % though this is redudant since the high-D sRA has magnitude 1.
%     dRA_lowD_prop_bySRA = bsxfun(@rdivide,abs(dRA_lowD),dRA_mag_unit);
%     
%     % compute proportion of dRA captured by lowD space. Since the sRAs are
%     % orthogonal, we find the proportion of dRA captured by taking the
%     % magnitude of the projected vector (i.e., sqrt of sum of squared
%     % components, summed across sRAs). We then normalize this magnitude by
%     % the magnitude of the high-D dRA, even though this is redudant since
%     % the high-D dRA has magnitude 1.
%     dRA_lowD_prop = sqrt(sum(dRA_lowD .^ 2,2)) ./ dRA_mag_unit;
%     
%     % check
%     if any(bsxfun(@minus,sqrt(sum(dRA.^2,2)),sqrt(sum(dRA_lowD.^2,2))) < 0)
%         error('projection onto dRA cannot be more than projection magnitude into 3D space')
%     end
%     
%     % Backsolve for the actual magnitude of the projected dRA by
%     % scaling the high-D magnitude by the proportion captured
%     dRA_lowD_mag = dRA_mag .* dRA_lowD_prop;
%     % Do the same for the magnitude of each component 
%     dRA_lowD_mag_bySRA = bsxfun(@times, dRA_lowD_prop_bySRA, dRA_mag);
% 
    
    %%%%%%%%%%%%%
    % NULL MODEL of random vector or space capturing linear represenations
    % (dRAs) of signlas of interest.
    dataTensor = Summary.oTDRSummary.dataTensor;
    [T, numNeurons, C] = size(dataTensor);
    bigA = reshape(permute(dataTensor,[1 3 2]), [], numNeurons);
    numSamples = 10000;
%     if isempty(sRA_rand)
        % Generate random samples
        [out] = sampleRandSubspaces(nSRA, cov(bigA), ...
            'vectorOverlap', numSamples, dRA_orig);
%     else
%         % compute variance explained by random samples
%         vecOL_rand = vectorOverlap(sRA_rand, dRA_orig);
%     end
%     
    
    % Initialize
    lowD_null(i).dRA_lowD_prop = NaN(numSamples,rT,1);
    lowD_null(i).dRA_lowD_prop_bySRA = NaN(numSamples,rT,nSRA);
    lowD_null(i).dRA_lowD_mag = NaN(numSamples,rT,1);
    lowD_null(i).dRA_lowD_mag_bySRA = NaN(numSamples,rT,nSRA);
    lowD_null(i).dRA_lowD_magDiff = NaN(numSamples,rT,1);
    lowD_null(i).dRA_lowD_magDiff_bySRA = NaN(numSamples,rT,nSRA);

    % Collect random samples
    for s = 1:numSamples
        lowD_null(i).dRA_lowD_prop(s,:,:) = out{s}.dRA_lowD_prop;
        lowD_null(i).dRA_lowD_prop_bySRA(s,:,:) = out{s}.dRA_lowD_prop_bySRA;
        lowD_null(i).dRA_lowD_mag(s,:,:) = out{s}.dRA_lowD_mag;
        lowD_null(i).dRA_lowD_mag_bySRA(s,:,:) = out{s}.dRA_lowD_mag_bySRA;
        lowD_null(i).dRA_lowD_magDiff(s,:,:) = out{s}.dRA_lowD_magDiff;
        lowD_null(i).dRA_lowD_magDiff_bySRA(s,:,:) = out{s}.dRA_lowD_magDiff_bySRA;
    end
    
    % compute p-values
    lowD_P(i).dRA_lowD_prop = sigTest(dRA_lowD_prop',lowD_null(i).dRA_lowD_prop,'lower','empirical')';
    lowD_P(i).dRA_lowD_mag = sigTest(dRA_lowD_mag',lowD_null(i).dRA_lowD_mag,'upper','empirical')';
    lowD_P(i).dRA_lowD_magDiff = sigTest(dRA_lowD_magDiff',lowD_null(i).dRA_lowD_magDiff,'upper','empirical')';
    lowD_P(i).dRA_lowD_prop_bySRA = NaN(size(dRA_lowD_prop_bySRA));
    lowD_P(i).dRA_lowD_mag_bySRA = NaN(size(dRA_lowD_mag_bySRA));
    lowD_P(i).dRA_lowD_magDiff_bySRA = NaN(size(dRA_lowD_magDiff_bySRA));
%     err_gauss = NaN(size(dRA_lowD_prop_bySRA));
%     err_gamma = NaN(size(dRA_lowD_prop_bySRA));
    % loop through each sRA
    for j = 1:nSRA
        % Right now we compare the measured values for sRA i to the null
        % distribution for null_RA 1 given dRA i. This is because the set
        % of null_RAs is ordered from most to least explanatory. Therefore,
        % we should compare sRA i to null_RA 1. However, the question
        % remains as to which null_RA to compare sRA k given dRA i (where k
        % ~= i). Should sRA k be compared to null_RA 2, null_RA 3?
        if i==1 && j==1
            warning(['When computing the difference in variance-related metrics ',... 
                'between an OFF-TARGET sRA-dRA pair, we have not determined ',...
                'the appropriate null distribution. This does not affect the ',...
                'main analysis, but will cause an error if attempting to plot ',...
                'p-values for off-target differences in the next code block']);
        end
        
        if i == j
            % for on-target sRA, compare to null sRA 1 (i.e., null vector
            % with highest alignment to dRA)
            foo = 1;
        else
            foo = j;
            continue % do not compute off-target p-value (see above wanring)
        end
        [lowD_P(i).dRA_lowD_prop_bySRA(:,j)] = sigTest(dRA_lowD_prop_bySRA(:,j)',lowD_null(i).dRA_lowD_prop_bySRA(:,:,foo),'lower','empirical');
        [lowD_P(i).dRA_lowD_mag_bySRA(:,j)] = sigTest(dRA_lowD_mag_bySRA(:,j)',lowD_null(i).dRA_lowD_mag_bySRA(:,:,foo),'upper','empirical');
        [lowD_P(i).dRA_lowD_magDiff_bySRA(:,j)] = sigTest(dRA_lowD_magDiff_bySRA(:,j)',lowD_null(i).dRA_lowD_magDiff_bySRA(:,:,foo),'upper','empirical');
    end
    %%%%%%%%%%%%%
    
end

%% PLOT magnitude of dRA projected into low-D space
%
% This is a test of how much signal related to the task-relevant variables
% is available for (linear) readout at each time bin by the dRAs. 
%
% To recreate Supp Fig 15 in Kimmel et al., 2020, use the following
% arguments:
%   Panels (a) and (d) (corresponding to rows 1 and 3 in the resulting
%   figure, respectively):  
%       bitPlotMagDiff = false;
%       bitMatchSRA2DRA = false;
%   Panels (b), (c), (e) and (f) (corresponding to rows 1, 2, 3 and 4 in
%   the resulting figure, respectively):
%       bitPlotMagDiff = true;
%       bitMatchSRA2DRA = true;
%   (Keep all other arguments as below)
       
% logical on whether to plot the magnitude of the dRA projected
% into the low-D space (==true) or the proportion of the dRA captured by
% the low-D space (==false).
bitPlotAbsMag = true;

% When plotting the magnitude of the low-D captured dRA (i.e.,
% bitPlotAbsMag == TRUE), logical on whether to plot the absolute magnitude
% (==FALSE) or the difference between high-D and captured magnitude
% (==TRUE).
bitPlotMagDiff = false;

% logical on whether, when plotting proportions, to plot on right-hand
% y-axis the p-value of proportion (TRUE) or magnitude of dRA in high-D
% space (FALSE). Only applies when bitPlotAbsMag == false
bitPlotPVal = true;

% logical on whether to plot p-value on separate axes (==TRUE) vs. on same
% axes but referenced to right-hand y-axis (==FALSE). P-value is not
% plotted when bitPlotPVal == FALSE (and thus bitPlotPValSepAx is
% irrelevant).
bitPlotPValSepAx = true;

% When plotting absolute magnitude of projection (bitPlotAbsMag == TRUE),
% logical on whether to plot the p-value that the difference in magnitude
% between high-D dRA and dRA as projected into low-D space is greater than
% chance (==TRUE) or instead that the absolute magnitude of projected
% vector (not difference) is greater than chance (==FALSE). The former
% assesses whether significant magnitude was missed by the low-D space,
% while the latter assesses whether the low-D space captures a significance
% magnitude. NOTE: the chance that the difference of magnitude is greater
% than chance difference is equivlent to the chance that the proportion of
% magnitude captured is less than the chance proportion.
bitMagPValDiff = true;

% logical on whether to limit sRA-specific plots to the on-target sRA/dRA
% pair (==TRUE). E.g., plot the magnitude/proportion of the benefit dRA
% captured by the Benefit sRA (==TRUE), or instead captured by each of the
% sRAs (==FALSE)
bitMatchSRA2DRA = true;

% define color for high-D magnitude
lcHiD = [0.5 0.5 0.5];

%%%%%%%%%%%%%%%%%%%%%

% set logicals based on higher-level logical
if ~bitPlotPVal
    % if not plotting p-vals, no need to specify separate axes
    bitPlotPValSepAx = false;
end
if ~bitPlotAbsMag
    % if not plotting magnitude of projections, then we cannot plot
    % difference of magntidue
    bitPlotMagDiff = false;
end

% checks
if bitPlotMagDiff && bitMatchSRA2DRA
    % The p-values for dRA_k projected onto sRA_i when i ~= k is not
    % well-defined. It's unclear which null distribution should be used for
    % these off-target sRAs (see above code block).
    error(['We do not currently support p-values for the ',...
        'difference between OFF-TARGET dRA-sRA pairs because ',...
        'the null distribution is not well defined (see warning ',...
        'in above code block). Restrict your analysis to either ',...
        'the absolute magnitude of the dRA (set bitPlotMagDiff = ',...
        'TRUE) or difference between on-target pairs (set ',...
        'bitMatchSRA2DRA = TRUE).'])
end

% get colors
lc = colorCategorical(length(condLabel));

% make figure
figure;
whitebg(gcf,bgColor)
figName = ['dRA magnitude in low-D space - ',analName];
set(gcf,'Color',bgColor,'Name',figName);

% get times
x = Summary.(tdr).Times.regressTimes';
dt = x(2) - x(1);

% determine number of rows
if bitPlotPVal && bitPlotPValSepAx
    nRow = 4;
else
    nRow = 2; % 1 row for projection of dRA on to each sRA; 1 row for projection of dRA into 3D space
end
nCol=length(condLabel);
    
clear ax axH h h1 h2 foo goo 

%%%%%% LOOP THROUGH SIGNALS FOR PLOTTING
for i = 1:length(condLabel)

    % extract output as separate vars
    dRA_lowD_prop = lowD_val(i).dRA_lowD_prop;
    dRA_lowD_prop_bySRA = lowD_val(i).dRA_lowD_prop_bySRA;
    dRA_lowD_mag = lowD_val(i).dRA_lowD_mag;
    dRA_lowD_mag_bySRA = lowD_val(i).dRA_lowD_mag_bySRA;
    dRA_lowD_magDiff = lowD_val(i).dRA_lowD_magDiff;
    dRA_lowD_magDiff_bySRA = lowD_val(i).dRA_lowD_magDiff_bySRA;

    % extract dRA magnitude vs. time
    dRA_mag = Summary.(tdr).(normFN)(:,i)';

    %%%%%%%%%%%%%
    % ROW 1 -- projection of dRA onto each sRA
    % subplot
    axH = subplot(nRow,nCol,i);
    hold on
    % title
    title(sprintf('%s dRA',condLabel{i}),'FontSize',18);
    
    % determine whether to plot difference or proportion of dRA capture
    if bitPlotAbsMag
        if bitPlotMagDiff
            % magnitude of high-D dRA captured by each sRA
            y = bsxfun(@minus,dRA_mag',dRA_lowD_mag_bySRA);
            yLabelStr = 'Magnitude of high-D - sRA-captured dRA (a.u)';
        else
            % magnitude of high-D dRA captured by each sRA
            y = dRA_lowD_mag_bySRA;
            yLabelStr = 'Magnitude of dRA on sRA (a.u.)';
        end
        
        if bitMagPValDiff
            y2 = log10(lowD_P(i).dRA_lowD_magDiff_bySRA); % p-value
            y2LabelStr = 'log_{10}1 - P(magnitude diff > chance)';
        else
            y2 = log10(lowD_P(i).dRA_lowD_mag_bySRA); % p-value
            y2LabelStr = 'log_{10}1 - P(magnitude > chance)';
        end
    else
        % proportion of dRA captured: projection onto each sRA, normalized
        % by dRA magnitude (which should be 1)
        y = dRA_lowD_prop_bySRA;
        yLabelStr = 'Proportion of dRA captured by sRA';
        ylim([0 1]);
        y2 = log10(lowD_P(i).dRA_lowD_prop_bySRA); % p-value
        y2LabelStr = 'log_{10}1 - P(proportion < chance)';
    end
    
    % limit plotting to when the sRA signal matches the dRA signal (i.e.,
    % on-target sRA). 
    if bitMatchSRA2DRA 
        y = y(:,i);
        y2 = y2(:,i);
    end
    
    % plotting p-value on separate axes
    if bitPlotPValSepAx
        % plot magnitude/proportion
        h = plot(x,y);
        
        % if plotting difference of magnitude, include null distribution
        % percentiles, but only for the on-target sRA
        if bitPlotMagDiff
            % if matching sRA to dRA, only plot percentiles for on-target
            % signals
            if bitMatchSRA2DRA
                foo = 1;
                goo = i;
            else
                % Note: this section will never be reached. See error
                % message above.
                
                % plot percentiles for all signals. Reorder null
                % distributions such that the first null distribution is
                % plotted in the color for the current signal (i).
                foo = 1:length(condLabel);
                goo = circshift(foo,-i+1);                
            end
            
            for j = foo
                hPrc = plot(x,...
                    prctile(lowD_null(i).dRA_lowD_magDiff_bySRA(:,:,j),[50 99])',...
                    'Color',lc(goo(j),:));
                set(hPrc(1),'LineStyle','--');
                set(hPrc(2),'LineStyle',':');
            end
        end
        
        % add additional row for p-value
        axH(2) = subplot(nRow,nCol,i+nCol);
        hold on
        
        % plot p-value on second axes:
        h2 = plot(x,y2,'LineWidth',2);
        
    elseif bitPlotPVal
        % plotting p-value on same axes, referencing right-hand y-axis
        [axH,h,h2] = plotyy(x,y,x,y2);
        % set line style for p-values
        set(h2,'LineStyle','--');
        % make sure tick marks adjust
        set(axH(2),'YTickMode','auto');
    else % no p-values
        % plot magnitude/proportion
        h = plot(x,y);
        
    end
    
    % link main and p-value axes
    linkaxes(axH,'x');
    
    % set line width
    set(h,'LineWidth',2);
    
    % if p-values
    if bitPlotPVal
        % invert p-value y-axis
        set(axH(2),'YDir','reverse');
    end
    
    % set each line color
    clear foo
    if bitMatchSRA2DRA
        % gather handles
        if exist('h2','var') && isobject(h2)
            foo = [h,h2];
        else
            foo = h;
        end
        set(foo,'Color',lc(i,:))
        clear foo
        % make on-target sRA thick
        set(h,'LineWidth',3);
    else
        for j=1:length(h)
            % gather handles
            if exist('h2','var') && isobject(h2)
                foo = [h(j),h2(j)];
            else
                foo = h(j);
            end
            set(foo,'Color',lc(j,:))
            clear foo
            % make on-target sRA thick
            if i==j
                set(h(j),'LineWidth',3);
            end
        end
    end
    
    axes(axH(1)); % make main axes focus
    
    % add high-D dRA magnitude
    if bitPlotAbsMag && ~bitPlotMagDiff
        hold on;
        plot(x,dRA_mag,'Color',lcHiD,'LineWidth',2);
    end
    
    % appearance stuff
    yLim = get(gca,'YLim');
    set(axH,'box','off','FontSize',16);
    set(axH,'XLim',[min(x), max(x)]);
    set(axH,'YColor',lineColor); % for top row, all y-axes are same color
    % line at x = 0
    line([0 0],yLim,'Color',lineColor);
    
    if i == 2 
        if bitPlotPValSepAx
            foo = axH(2);
        else
            foo = axH(1);
        end
        xlabel(foo,'Time from offer (s)','FontSize',18);
        clear foo
    end
    if i == 1
        ylabel(yLabelStr,'FontSize',18)
    end
    
    % p-value axes appearance
    if bitPlotPVal 
        if (i == 1 && bitPlotPValSepAx) || (i == 3 && ~bitPlotPValSepAx)
            ylabel(axH(2),y2LabelStr,'FontSize',18)
        end
    
        % set limits of p-value axis if plotting magnitude (where point is to
        % find small p-values indicating larger than chance magnitudes). do not
        % set limits if plotting proportion (where point is to find large
        % p-values indicating comparable to chance proportions)
        yLim2 = get(axH(2),'YLim');
        if bitPlotAbsMag && ~bitMagPValDiff
            set(axH(2),'YLim',[min(yLim2(1),-4) 0]);
        else
            set(axH(2),'YLim',[min(yLim2(1),-2) 0]);
        end
        
        % plot horizontal line at p = log10(0.01) = -2
        line(axH(2),[min(x), max(x)],[-2 -2],'LineStyle','--','Color',lineColor);
        
        % plot vertical line at x = 0
        yLim2 = get(axH(2),'YLim');
        line(axH(2),[0 0],yLim2,'Color',lineColor);
        
        % return focus to second axes
        axes(axH(2));
        
    end
    
    clear axH h h2 y y2 y3 yLim2
    
    %%%%%%%%%%%%%%%%%%
    % ROW 2 -- projection of dRA into 3D space
    % subplot
    % determine axes offset
    if bitPlotPValSepAx
        offset = nCol*2;
    else
        offset = nCol;
    end
    axH = subplot(nRow,nCol,i+offset);
    hold on;
    
    % determine whether to plot difference or proportion of dRA capture
    if bitPlotAbsMag
        if bitPlotMagDiff
            % difference between magnitude of high-D dRA and the portion
            % captured by low-D space 
            y = dRA_mag' - dRA_lowD_mag;
            yLabelStr = 'Magnitude of high-D - low-D captured dRA (a.u.)';
        else
            % magnitude of high-D dRA captured by low-D space
            y = dRA_lowD_mag;
            yLabelStr = 'Magnitude of dRA in low-D space (a.u.)';
        end
        if bitMagPValDiff
            y2 = log10(lowD_P(i).dRA_lowD_magDiff); % p-value
            y2LabelStr = 'log_{10}1 - P(magnitude diff > chance)';
        else
            y2 = log10(lowD_P(i).dRA_lowD_mag); % p-value
            y2LabelStr = 'log_{10}1 - P(magnitude > chance)';
        end
    else
        % set y limit for proportion
        ylim([0 1]);
        
        % proportion of dRA captured by low-D space
        y = dRA_lowD_prop;
        yLabelStr = 'Proportion of dRA captured by low-D space';
        y2 = log10(lowD_P(i).dRA_lowD_prop); % p-value
        y2LabelStr = 'log_{10}1 - P(proportion < chance)';
        
        % absolute magnitude info
        y3 = dRA_mag; % dRA magnitude
        y3LabelStr = 'dRA Magnitude (a.u.)';
        
    end
    
        
    if bitPlotPVal && ~bitPlotPValSepAx
        % plotting p-value on same axes, referencing right-hand y-axis
        [axH,h,h2] = plotyy(x,y,x,y2);
        % set line style for p-values
        set(h2,'LineStyle','--');
        % right-hand y-axis ticks
        set(axH(2),'YTickMode','auto');
        
    else % no p-values on same axis
        if ~bitPlotAbsMag
            % if plotting proportion, include absolute magnitude on
            % right-hand y-axis
            [axH_mag,h,hMag] = plotyy(x,y,x,y3);
            % collect axes handle
            axH(1) = axH_mag(1);
            % link axes
            linkaxes(axH_mag,'x');
            
            % set appearance of absolute magnitude line and right-hand axes
            set(hMag,'LineWidth',2,'Color',lcHiD);
            set(axH_mag(2),'YColor',lcHiD,'box','off','FontSize',16);
            
            % add ylabel
            if i==3
                ylabel(axH_mag(2),y3LabelStr,'FontSize',18);
            end
        else
            % plot magnitude
            h = plot(x,y);
            
            % if plotting difference of magnitude, include null
            % distribution percentiles
            if bitPlotMagDiff
                hPrc = plot(x,prctile(lowD_null(i).dRA_lowD_magDiff,[50 99])',...
                    'Color',lcHiD);
                set(hPrc(1),'LineStyle','--');
                set(hPrc(2),'LineStyle',':');
                
            else
                % add high-D dRA magnitude only if not plotting difference
                hold on;
                plot(x,dRA_mag,'Color',lcHiD,'LineWidth',2);
            end
        end 
        
        % plot pVals on separate axis
        if bitPlotPVal 
            
            % add additional row for p-value
            axH(2) = subplot(nRow,nCol,i+offset+nCol);
            hold on
            
            % plot p-value on second axes:
            h2 = plot(x,y2,'LineWidth',2);
            
        else
            % dummy handle
            h2 = [];
        end
    end
    
    % link main and p-value axes
    linkaxes(axH,'x');
    
    % set line width
    set(h,'LineWidth',2);
    
    % if p-values
    if bitPlotPVal
        % invert p-value y-axis
        set(axH(2),'YDir','reverse');
    end
    
    % set line color
    set([h h2],'Color',lineColor)
    
%     axes(axH(1)); % make main axes focus
    
    % appearance stuff
    yLim = get(axH(1),'YLim');
    set(axH,'box','off','FontSize',16);
    set(axH,'XLim',[min(x), max(x)]);
    set(axH,'YColor',lineColor);
    
    % line at x = 0
    line(axH(1),[0 0],yLim,'Color',lineColor);
    
    if i == 2 
        if bitPlotPValSepAx
            foo = axH(2);
        else
            foo = axH(1);
        end
        xlabel(foo,'Time from offer (s)','FontSize',18);
        clear foo
    end
    if i == 1
        ylabel(axH(1),yLabelStr,'FontSize',18)
    end
    
    % p-value axes appearance
    if bitPlotPVal 
        if (i == 1 && bitPlotPValSepAx) || (i == 3 && ~bitPlotPValSepAx)
            ylabel(axH(2),y2LabelStr,'FontSize',18)
        end
    
        % set limits of p-value axis if plotting magnitude (where point is to
        % find small p-values indicating larger than chance magnitudes). do not
        % set limits if plotting proportion (where point is to find large
        % p-values indicating comparable to chance proportions)
        yLim2 = get(axH(2),'YLim');
        if bitPlotAbsMag && ~bitMagPValDiff
            set(axH(2),'YLim',[min(yLim2(1),-4) 0]);
        else
            set(axH(2),'YLim',[min(yLim2(1),-2) 0]);
        end
        
        % plot horizontal line at p = log10(0.01) = -2
        line(axH(2),[min(x), max(x)],[-2 -2],'LineStyle','--','Color',lineColor);
        
        % plot vertical line at x = 0
        yLim2 = get(axH(2),'YLim');
        line(axH(2),[0 0],yLim2,'Color',lineColor);
        
        % return focus to second axes
%         axes(axH(2));
        
    end
    
    clear axH h h2 y y2 y3
end

% annotation
[figTitleHand] = annotateFig(figName,mfilename);