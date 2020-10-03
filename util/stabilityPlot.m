function [fH,fH2] = stabilityPlot(angleProfile,varargin)
% TODO -- update docstring
% 
% [fH,fH2] = stabilityPlot(angleProfile,varargin)
%
% Plots TWO figures:
% (1) "Topo" maps of dRA angles, with each dRA's set of angles plotted as
%   a separate trace.
% (2) Average angle profile aligned to time of reference dRA, including
%   1st, 2nd and 3rd order polynomial fits to the profile.
%
% Accepts ANGLEPROFILE, which is a struct compiled in CB_gamal_data_plot
% from the output of CB_stability. It has S elements for the S signals.
% Each element contains two structs, .main and .fitO, which are the output
% vars from CB_stability.
%
% Returns figure handles FH
%
% Daniel Kimmel, 06 Feb 2016

%% default params

% set which class of angle fit
% 'decay' for exponential decay fits
% 'step' for step function fits
% 'step_trans' for step function TRANS fit (angles that are anti-similar)
fitClass2u = 'step'; 

% logical on whether to plot auto-angle (within signal) or cross-angles
% (between signal)
bitAuto = true;

% logical on whether to plot fits from 2-step fit functions
bit2Step = false;

%%% plotting stuff
% logical whether to horizontally align each row of angles such that the
% time of the reference dRA is time = 0 (TRUE). When FALSE, horizontal
% position refers to the absolute time of the comparison dRA (and time = 0
% refers to the task event to which the original PTSHs were aligned).
bitAlignRows = 1; 

% set vertical offset (in degrees) between rows
offset = -30;

% hide the peak angle
bitHidePeak = 0;

% logical on whether to plot error shading as simple lines for output to
% illustrator (applies to average profile only).
bit4Illustrator = false;

% fontSize
fontSize = 16;
textColor = 'k';
lineColor = 'k';
bgColor = 'w';
transparency = 0.4;

%% collect optionally provided parameters
warnopts(assignopts(who, varargin));

%% prepare vars

switch fitClass2u
    case 'decay'
        s2u = 'main_decay';
        f2u = 'fitO_decay';
    case 'step'
        s2u = 'main_step';
        f2u = 'fitO_step';
    case 'step_trans'
        s2u = 'main_step_trans';
        f2u = 'fitO_step_trans';
    otherwise
        error('Fit class %s not recognized',fitClass2u)
end

if strcmp('step',fitClass2u(1:4)) && bit2Step
    s2u = [s2u,'_2'];
end

%% determine which sets of angles to plot

% initialize
pos2Plot = [];
signalSet = [];

% loop through analyses in angleProfile
for i = 1:length(angleProfile)
    if (bitAuto && length(unique(angleProfile(i).(s2u).signalN)) == 1) || ...
            (~bitAuto && length(unique(angleProfile(i).(s2u).signalN)) > 1)
        % include:
        pos2Plot(end+1) = i;
        signalSet(end+1,:) = angleProfile(i).(s2u).signalN;
    end
end

if bitAuto
    nRow = 1;
    nCol = length(pos2Plot);
else
    nRow = 2;
    nCol = ceil(length(pos2Plot)/2);
end

%% plot raw angle and fits vs. time for each dRA (e.g., "tree roots")

% figure for profile per reference dRA
figName = 'Angle profile';
fH = figure('Name',figName);
whitebg(gcf,bgColor);
set(gcf,'Color',bgColor);

axN = 0;
xLimAll = NaN(length(pos2Plot),2);
yLimAll = NaN(length(pos2Plot),2);

for i = pos2Plot

    axN = axN + 1;
    
    % extract vars
    out = angleProfile(i).(s2u);
    if isfield(angleProfile(i),f2u)
        fo = angleProfile(i).(f2u);
    else
        fo = [];
    end
    t = angleProfile(i).(s2u).tMaster;
    dt = t(2) - t(1);
    theta = angleProfile(i).(s2u).thetaMaster;
    
%     switch i
%         case 1
%             strTitle = 'Benefit';
%         case 2
%             strTitle = 'Choice';
%         case 3
%             strTitle = 'Expected Reward';
%     end

    strTitle = angleProfile(i).(s2u).signalDesc;

    % determine theta offset used for fits, which is stored in upper bound
    % of fits
    if strcmp(fitClass2u,'step') || strcmp(fitClass2u,'step_trans')
        thetaOffset = angleProfile(i).(s2u).aBound(:,2);
    else
        thetaOffset = 90;
    end
    
    % determine scalar
    if isfield(angleProfile(i).(s2u),'angleScale') && ~isempty(angleProfile(i).(s2u).angleScale)
        aS = angleProfile(i).(s2u).angleScale;
    else
        aS = 1;
    end
    
    axH(axN) = subplot(nRow,nCol,axN);
    hold on

    % time shift
    if bitAlignRows 
        shiftT = t';
    else
        shiftT = zeros(length(t),1);
    end

    % shift vertically
    shiftVec = ([0:size(out.theta,1)-1]' * offset);
    theta4Fit = out.theta + repmat(shiftVec,1,size(out.theta,2));
    theta4All = aS*(repmat(thetaOffset,1,size(theta,2)) - theta) + repmat(shiftVec,1,size(theta,2));

    % plot all theta
%     plot((repmat(t,size(theta4All,1),1)-repmat(shiftT,1,size(theta4All,2)))',theta4All','Color',[0.7 0.7 0.7]);
    hold on;
    % plot theta used for fits
    plot((out.t - repmat(shiftT,1,size(out.t,2)))',theta4Fit','Color',lineColor);
        
    %%% plot FITS
    if isfield(angleProfile(i).(s2u),'fitMethod') && ...
            strcmp('paramSearch',angleProfile(i).(s2u).fitMethod)
        for j = 1:size(out.theta,1)
            stepFnBest = zeros(length(t),1);
            % accommodate multiple steps
            for k = 1:size(angleProfile(i).(s2u).b,3)
                foo = zeros(length(t),1);
                startPos = find(t==angleProfile(i).(s2u).b(j,2,k));
                width = angleProfile(i).(s2u).b(j,3,k)/dt;
                if round(width) - width > eps
                    error('Expect integer width')
                else
                    width = round(width);
                end
                foo(startPos:startPos+width-1) = angleProfile(i).(s2u).b(j,1,k);
                stepFnBest = stepFnBest + foo;
                clear foo
            end
            bitExclude = isnan(out.theta(j,:));
            % time shift
            if bitAlignRows 
                shiftT_fit = t(j);
            else
                shiftT_fit = 0;
            end            
            if ~all(bitExclude)
                stairs(t(~bitExclude) - shiftT_fit, stepFnBest(~bitExclude) + shiftVec(j), 'r');
            end
        end
        
    else
        for j = 1:length(fo)
            % skip plotting fit if all coefficients are set to 0 (this is a
            % sign that the fit was eliminated)
            if all(isnan(out.b(j,:)))
                continue
            end
            % position of times used:
            pos2Use = isfinite(out.t(j,:));
            % extract times used
            tT = out.t(j,pos2Use);
            % find time range used:
            tR = range(tT);
            % build times for eval:
            tE = linspace(0,tR,length(tT));
            
            t2Plot = tT;
            if bitAlignRows
                if isnan(out.t(j,j))
                    % when skipping peak value
                    t2Plot = tT - out.tMaster(j);
                else
                    t2Plot = tT - tT(1);
                end
            end
            % plot against true theta times, offset vertically
            stairs(t2Plot',fo(j).f(tE') + shiftVec(j),'Color','r');
        end
    end
    

    if i == ceil(length(angleProfile)/2)
        if bitAlignRows
            xlabel('Time from time of dRA (s)','FontSize',fontSize);
        else
            xlabel('Time from offer (s)','FontSize',fontSize);
        end
    end
    
    if i==1
        ylabel('90 - Theta (deg)','FontSize',fontSize);
    end
    
    % appearance
    set(gca,'FontSize',fontSize);
    % eliminate y labels
    set(gca,'YTick',[]);
    % store limits
    xLimAll(axN,:) = get(gca,'XLim');
    yLimAll(axN,:) = get(gca,'YLim');
    
    % remove axes:
%     set(gca,'XColor',get(gcf,'Color'));
%     set(gca,'YColor',get(gcf,'Color'));
    set(gca,'Visible','off');
    
    % Title axes. Don't use TITLE() function because will be removed with
    % axes
%     title(strTitle,'FontSize',fontSize);
    h = text(0.75,1,strTitle,'Units','normalized','FontSize',fontSize,...
        'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'Color',textColor);
    set(h,'Clipping','off');
    
    
    % make scale bars
    if axN==1
        xText = t(1) - shiftT(end)-2*dt;
        yText = shiftVec(end)+offset;
        h = line([xText xText+1],[yText yText],'Color',textColor); %horizontal line
        set(h,'Clipping','off');
        h = line([xText xText],[yText yText+90],'Color',textColor); %vertical line
        set(h,'Clipping','off');
        h = text(xText,yText+90,'90 deg','HorizontalAlignment','right',...
            'VerticalAlignment','bottom','Color',textColor,'FontSize',fontSize);
        set(h,'Clipping','off');
        h = text(xText+1,yText,'1 s','VerticalAlignment','top',...
            'Color',textColor,'FontSize',fontSize);
        set(h,'Clipping','off');
    end
    
    % label time of reference dRA
%     
%     set(gca,'XTick',[0],'XTickLabel',0)
%     xlabel('Time from reference dRA (s)','FontSize',fontSize);
    
    % force limits so that all subplots have same scale
%     xlim([xText t(end)+2*dt]);
%     ylim([shiftVec(end) 100]);
end

% apply axes limits to all subplots
xlim(axH,[min(xLimAll(:,1)) max(xLimAll(:,2))]);
ylim(axH,[min(yLimAll(:,1)) max(yLimAll(:,2))]);

% annotate original figure
annotateFig(figName,mfilename);

%% SECONDARY FIGURE (AVERAGE)

% figure for average profile
figName2 = 'Angle profile average';
fH2 = figure('Name',figName2);
whitebg(gcf,bgColor);
set(gcf,'Color',bgColor);

axN = 0;

for i = pos2Plot

    axN = axN + 1;
    
    % extract vars
    out = angleProfile(i).(s2u);
    t = angleProfile(i).(s2u).tMaster;
    dt = t(2) - t(1);
    
%     switch i
%         case 1
%             strTitle = 'Benefit';
%         case 2
%             strTitle = 'Choice';
%         case 3
%             strTitle = 'Expected Reward';
%     end

    strTitle = angleProfile(i).(s2u).signalDesc;

    %%% BUILD matrix of theta aligned to reference time
    thetaAlign2Ref = NaN(size(out.theta ,1),size(out.theta,1)*2-1);
    % loop through each row of theta, placing in matrix
    for q = 1:size(out.theta,1)
        thetaAlign2Ref(q,size(out.theta,1)-q+1:2*size(out.theta,1)-q) = ...
            out.theta(q,:);
    end
    % build time vector
    tAlign2Ref = [out.tMaster, out.tMaster(2:end) + range(out.tMaster(1:end))];
%     % note that first row(s) of out.t may not contain data, so find first
%     % row:
%     posTRow1 = find(any(isfinite(out.t),2),1,'first');
%     tAlign2Ref = [out.t(end,1:end-1), out.tMaster(end), out.t(1,2:end) + range(out.t(1,:))+dt];
    % center time vector
    tAlign2Ref = tAlign2Ref - out.tMaster(end);
%     tAlign2Ref = tAlign2Ref - out.t(1,end);
    % bitmask for positive (dt >= 0) portion
    bitPos = tAlign2Ref >0;
    % bitmask for reference dRA(t) s.t., t >=0
    bitRAPos = out.tMaster > 0;
    
    % compute mean, std -- limit to positive portion
    thetaAlign2Ref_mean = nanmean(thetaAlign2Ref(bitRAPos,bitPos));
    thetaAlign2Ref_std = nanstd(thetaAlign2Ref(bitRAPos,bitPos));
    thetaAlign2Ref_n = sum(isfinite(thetaAlign2Ref(bitRAPos,bitPos)));
    
    % fit line to positive portion
    foo = thetaAlign2Ref(bitRAPos,bitPos);
    goo = repmat(tAlign2Ref(bitPos),size(foo,1),1);
    foo = foo(:);
    goo = goo(:);
    bitElim = isnan(foo) | isnan(goo);
    [coeff1,stats1] = polyfit(goo(~bitElim),foo(~bitElim),1);
    t4Fit = unique(goo(~bitElim));
    
%     coeff = regress(foo(:),[goo(:) ones(size(foo(:)))]);
    
    % fit 2nd order poly to positive portion
    [coeff2,stats2] = polyfit(goo(~bitElim),foo(~bitElim),2);

    % fit 3rd order poly to positive portion
    [coeff3,stats3] = polyfit(goo(~bitElim),foo(~bitElim),3);

    clear foo goo
    
    % COMPARE 1ST AND 2ND ORDER POLY
    % compare first and 2nd order polys with f-test ratio see
    % http://www.curvefit.com/2_models__1_dataset.htm also see
    % http://www.graphpad.com/support/faqid/1765/ compute ratio of percent
    % change in SS to percent change in DF:
    fRatio = ...
        ((stats1.normr^2 - stats2.normr^2) / ...
        stats2.normr^2) / ...
        ((stats1.df - stats2.df) / ...
        stats2.df);
    % compute p value using f distribution where the f ratio is the
    % test statistic and the numerator as DF1 - DF2 degrees of
    % freedom and the denomenator has DF2 degrees of freedom.
    pFRatio = ...
        1 - fcdf(fRatio,...
        stats1.df - stats2.df,...
        stats2.df);

    % COMPARE 1ST AND 3RD ORDER POLY
    fRatio13 = ((stats1.normr^2 - stats3.normr^2) / ...
        stats3.normr^2) / ((stats1.df - stats3.df) / stats3.df);
    pFRatio13 = 1 - fcdf(fRatio13,stats1.df - stats3.df,stats3.df);
    
    %%% PLOT average theta as function of time from reference theta on
    %%% separate figure
    axH2(axN) = subplot(nRow,nCol,axN);
    hold on
    if any(bitPos)
        [hL,hF] = errorArea(tAlign2Ref(bitPos),thetaAlign2Ref_mean,...
            thetaAlign2Ref_std ./ sqrt(thetaAlign2Ref_n),...
            lineColor,'transparency',transparency,...
            'bit4Illustrator',bit4Illustrator);
        set(hL,'LineWidth',2);
    end
    
    % add linear fit 
    plot(t4Fit,polyval(coeff1,t4Fit),'--','Color',lineColor,...
        'lineWidth',2);
    
    % add 2nd order fit 
    plot(t4Fit,polyval(coeff2,t4Fit),'-.','Color',lineColor,...
        'lineWidth',2);
    
    % add 3rd order fit 
    plot(t4Fit,polyval(coeff3,t4Fit),'.','Color',lineColor,...
        'lineWidth',2);
    
    % write p that 1st order poly is as good or better than 2nd order poly
    text(1,1,sprintf('F-test 1st %s 2nd order poly, p = %0.2g',char(8807),pFRatio),...
        'Units','Normalized','HorizontalAlignment','right','FontSize',12)
    % write p that 1st order poly is as good or better than 3rd order poly
    text(1,0.9,sprintf('F-test 1st %s 3rd order poly, p = %0.2g',char(8807),pFRatio13),...
        'Units','Normalized','HorizontalAlignment','right','FontSize',12)
    
    % figure plots as symmetric about x = 0, but really just focus on 1
    % half
    xlim([0 max(t4Fit)]);
    
    % limit to y>=0 (can be <0 because of linear fit)
    yLim = get(gca,'YLim');
    ylim([0 yLim(2)]);

    % appearance
    set(gca,'FontSize',fontSize);
    
    % Title axes
    title(strTitle);
    
    % y label
    if i==1
        ylabel('Mean similarity, 90 - theta \pms.e.m. (deg)')
    end
    % x label
    if i==pos2Plot(ceil(length(pos2Plot)/2))
        xlabel('Time from reference dRA (s)');
    end
  
end

    % plot on same figure as single dRA traces:
%     errorArea(axH2,tAlign2Ref - shiftT(end),nanmean(thetaAlign2Ref) + shiftVec(end) + offset*2,...
%         nanstd(thetaAlign2Ref) ./ sqrt(sum(isfinite(thetaAlign2Ref))),...
%         lineColor,'transparency',transparency,...
%         'bit4Illustrator',bit4Illustrator);


% annotate 
annotateFig(figName2,mfilename);

% link axes
linkaxes(axH2);

% %% Plot "topo" map of angles -- OLD WAY
%     
% % % try some basic fitting
% % for j = 1:size(theta,1)
% %     
% %     % fit left side:
% %     expMu(i,j,1) = expfit(theta(j,1:j));
% %     
% %     % fit right side:
% %     expMu(i,j,2) = expfit(theta(j,j:end));
% % end
% 
% % flip the theta vertically (so that angles decay with time)
% theta = 90-theta;
% 
% % hide the peak
% if bitHidePeak
%     foo = eye(size(theta,1));
%     foo(foo==1) = NaN;
%     theta = theta + foo;
%     clear foo
% end
% 
% % shift vertically
% foo = [1:size(theta,1)]' * offset;
% dataT = theta + repmat(foo,1,size(theta,2));
% clear foo
% 
% % shift horizontally
% tT = repmat(t',1,size(theta,2));
% foo = repmat(fliplr([0,t(1:end-1)]),size(theta,1),1);
% tT = tT + foo;
% clear foo
% 
% plot(tT,dataT','Color','k');
% 

