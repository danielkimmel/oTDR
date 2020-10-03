function oTDRFigures(Summary, SummaryPP)

all_times = Summary.Times.all_times;
sRA1 = Summary.sRA.RA1;
sRA2 = Summary.sRA.RA2;
sRA3 = Summary.sRA.RA3;

BenefitValues = Summary.Predictors(:,1);
ChoiceValues = Summary.Predictors(:,2);
clr = [26 173 150
        213 124 53
        123 136 192]./255;
%% Projections

genProjFigures(Summary.sRA.proj1, all_times, BenefitValues, ChoiceValues, 'Benefit Activity')
if isfield(Summary.sRA, 'proj1Single')
   plot(all_times, Summary.sRA.proj1Single,'--', 'linewidth', 2, 'color', [1.0000    0.2118    0.1373])
end
genProjFigures(Summary.sRA.proj2, all_times, BenefitValues, ChoiceValues, 'Choice Activity')
genProjFigures(Summary.sRA.proj3, all_times, BenefitValues, ChoiceValues, 'Expected Reward Activity')
if isfield(Summary, 'sRA_PP') 
genProjFigures(Summary.sRA_PP.proj1, all_times, BenefitValues, ChoiceValues, 'Benefit Activity')
genProjFigures(Summary.sRA_PP.proj2, all_times, BenefitValues, ChoiceValues, 'Choice Activity')
genProjFigures(Summary.sRA_PP.proj3, all_times, BenefitValues, ChoiceValues, 'Expected Reward Activity') 
end
%% Variance
tmin = all_times(1);
tmax = all_times(end);
figure
yTicks = [0 round(max(Summary.sRA.varAnalysis.VE_RA(:))*1.01)];
hold on
ylabel('variance explained (%)','FontSize',16)
set(gca,'ytick',yTicks)
set(gca,'ylim',[yTicks(1)-0.5 yTicks(end)])
xlabel('Time (ms)','FontSize',16)
set(gca,'xtick',[tmin 0 tmax])
set(gca,'xlim',[tmin tmax])
set(gca,'FontSize',16)
h2 = plot(all_times, Summary.sRA.varAnalysis.VE_RA(:,1), 'color', clr(1,:),'linewidth',2);
h3 = plot(all_times, Summary.sRA.varAnalysis.VE_RA(:,2), 'color', clr(2,:),'linewidth',2);
h4 = plot(all_times, Summary.sRA.varAnalysis.VE_RA(:,3), 'color', clr(3,:),'linewidth',2);
if exist('SummaryPP', 'var') 
    plot(SummaryPP.Times.all_times, SummaryPP.sRA.varAnalysis.VE_RA(:,1), 'color', clr(1,:),'linewidth',2,'lineStyle', '-.');
    plot(SummaryPP.Times.all_times, SummaryPP.sRA.varAnalysis.VE_RA(:,2), 'color', clr(2,:),'linewidth',2,'lineStyle', '-.');
    plot(SummaryPP.Times.all_times, SummaryPP.sRA.varAnalysis.VE_RA(:,3), 'color', clr(3,:),'linewidth',2,'lineStyle', '-.');
end
plot( [0 0], [yTicks(1) yTicks(end)],'k-.', 'LineWidth', 2)
legend([h2(1) h3(1) h4(1)],[{['Benefit']},{['Choice']},{['Expected Reward']}])


tmin = all_times(1);
tmax = all_times(end);
figure
hold on;
h0 = fill([all_times(1) all_times(end)  all_times(end) all_times(1)],[1 1 0.05 0.05], [0.5 0.5 0.5]);
h2 = plot(all_times, Summary.sRA.varAnalysis.pV_RA(:,1), 'color',clr(1,:), 'LineWidth', 2); 
h3 = plot(all_times, Summary.sRA.varAnalysis.pV_RA(:,2), 'color',clr(2,:), 'LineWidth', 2); 
h4 = plot(all_times, Summary.sRA.varAnalysis.pV_RA(:,3), 'color',clr(3,:), 'LineWidth', 2); 
if exist('SummaryPP', 'var') 
    plot(SummaryPP.Times.all_times, SummaryPP.sRA.varAnalysis.pV_RA(:,1), 'color',clr(1,:), 'LineWidth', 2,'lineStyle', '-.');
    plot(SummaryPP.Times.all_times, SummaryPP.sRA.varAnalysis.pV_RA(:,2), 'color',clr(2,:), 'LineWidth', 2,'lineStyle', '-.');
    plot(SummaryPP.Times.all_times, SummaryPP.sRA.varAnalysis.pV_RA(:,3), 'color',clr(3,:), 'LineWidth', 2,'lineStyle', '-.');
end
plot( [0 0], [yTicks(1) yTicks(end)],'k-.', 'LineWidth', 2)
xlabel('Time (ms)','FontSize',16)
set(gca,'xtick',[tmin 0 tmax])
set(gca,'xlim',[tmin tmax])
ylabel('p-value Var. w.r.t. chance (a.u.)','FontSize',16)
set(gca,'FontSize',16)
legend([h2(1) h3(1) h4(1) h0(1)],[{['Benefit']},{['Choice']},{['Expected Reward']}, {'Chance'}])
set(gca,'yscale','log')
set(gca,'YDir','reverse');
yTicks = get(gca,'yTick');
set(gca,'ylim',[min(Summary.sRA.varAnalysis.pV_RA(:)) 1]);
plot( [0 0], [yTicks(1) yTicks(end)],'k-.', 'LineWidth', 2)

%% RSV and ISV in one fig

sclev = 10*ceil(max(Summary.sRA.varAnalysis.RSVE_RA(:))./10);
tmin = all_times(1)-0.05;
tmax = all_times(end)+0.05;
yTicks = [0:10:sclev];
figure
hold on;
hb1 = plot(all_times, Summary.sRA.varAnalysis.RSVE_RA(:,1,1) , 'color',clr(1,:), 'LineWidth', 2); 
hb2 = plot(all_times, Summary.sRA.varAnalysis.ISVE_RA(:,1,1) ,':', 'color',clr(1,:), 'LineWidth', 2); 
hc1 = plot(all_times, Summary.sRA.varAnalysis.RSVE_RA(:,2,2) , 'color',clr(2,:), 'LineWidth', 2); 
hc2 = plot(all_times, Summary.sRA.varAnalysis.ISVE_RA(:,2,2) ,':', 'color',clr(2,:), 'LineWidth', 2); 
he1 = plot(all_times, Summary.sRA.varAnalysis.RSVE_RA(:,3,3) , 'color',clr(3,:), 'LineWidth', 2); 
he2 = plot(all_times, Summary.sRA.varAnalysis.ISVE_RA(:,3,3) ,':', 'color',clr(3,:), 'LineWidth', 2); 

plot( [0 0], [yTicks(1) yTicks(end)],'k-', 'LineWidth', 2)
set(gca,'ytick',yTicks)
set(gca,'ylim',[yTicks(1) yTicks(end)])
xlabel('Time (s)','FontSize',16)
set(gca,'xtick',[tmin 0 tmax])
set(gca,'xlim',[tmin tmax])
ylabel('RSV/ISV Percent of total variance')
set(gca,'FontSize',16)
% legend([h2(1) h3(1) h4(1) h1(1) h0(1)],[{'Benefit'},{'Choice'},{'Expected Reward'}, {'PCA projections'}, {'Chance'}])
legend([hb1(1) hc1(1) he1(1)],[{'Benefit'}, {'Choice'}, {'Expected Reward'}])
legend boxoff


sclep = floor(min(log10(Summary.sRA.varAnalysis.pRSV_RA(:))));
sclep = floor((sclep/2))*2;
tmin = all_times(1)-0.05;
tmax = all_times(end)+0.05;
yTicks = [10.^[max(sclep, -10):0]];
figure
hold on;
hb1 = plot(all_times, Summary.sRA.varAnalysis.pRSV_RA(:,1,1) , 'color',clr(1,:), 'LineWidth', 2); 
hb2 = plot(all_times, Summary.sRA.varAnalysis.pISV_RA(:,1,1) ,':', 'color',clr(1,:), 'LineWidth', 2); 
hc1 = plot(all_times, Summary.sRA.varAnalysis.pRSV_RA(:,2,2) , 'color',clr(2,:), 'LineWidth', 2); 
hc2 = plot(all_times, Summary.sRA.varAnalysis.pISV_RA(:,2,2) ,':', 'color',clr(2,:), 'LineWidth', 2); 
he1 = plot(all_times, Summary.sRA.varAnalysis.pRSV_RA(:,3,3) , 'color',clr(3,:), 'LineWidth', 2); 
he2 = plot(all_times, Summary.sRA.varAnalysis.pISV_RA(:,3,3) ,':', 'color',clr(3,:), 'LineWidth', 2); 
plot(all_times, 0.05*ones(size(all_times)),'k')
plot( [0 0], [yTicks(1) yTicks(end)],'k-', 'LineWidth', 2)
set(gca,'ytick',yTicks)
set(gca,'ylim',[yTicks(1) yTicks(end)])
xlabel('Time (s)','FontSize',16)
set(gca,'xtick',[tmin 0 tmax])
set(gca,'xlim',[tmin tmax])
ylabel('p-value of RSV/ISV')
set(gca,'FontSize',16)
% legend([h2(1) h3(1) h4(1) h1(1) h0(1)],[{'Benefit'},{'Choice'},{'Expected Reward'}, {'PCA projections'}, {'Chance'}])
legend([hb1(1) hc1(1) he1(1)],[{'Benefit'},{'Choice'}, {'Expected Reward'}])
legend boxoff
set(gca,'yscale','log')
set(gca,'YDir','reverse');


%%

%% between signal correlation (mixedness)

figure;
subplot(1,3,1)
scatter(Summary.sRA.RA1Star,Summary.sRA.RA2Star)
xlabel('Benefit'); ylabel('Choice');
[r,p] = corr(Summary.sRA.RA1Star,Summary.sRA.RA2Star);
title(sprintf('r = %0.2g, p = %0.2g',r,p));
axis equal
axis square

subplot(1,3,2)
scatter(Summary.sRA.RA1Star,Summary.sRA.RA3Star)
xlabel('Benefit'); ylabel('Expected Reward');
[r,p] = corr(Summary.sRA.RA1Star,Summary.sRA.RA3Star);
title(sprintf('r = %0.2g, p = %0.2g',r,p));
axis equal
axis square

subplot(1,3,3)
scatter(Summary.sRA.RA2Star,Summary.sRA.RA3Star)
xlabel('Choice'); ylabel('Expected Reward');
[r,p] = corr(Summary.sRA.RA2Star,Summary.sRA.RA3Star);
title(sprintf('r = %0.2g, p = %0.2g',r,p));
axis equal
axis square

end
function genProjFigures(A_projVOI,  all_times, BenefitValues, ChoiceValues, ActivityType)


A_projVOI_0 = A_projVOI(:, ChoiceValues<1);
A_projVOI_1 = A_projVOI(:, ChoiceValues>=1);
BenefitValues_0 = BenefitValues(ChoiceValues<1);
[BenefitValues_0, ix] = sort(BenefitValues_0,'ascend');
A_projVOI_0 = A_projVOI_0(:,ix);
BenefitValues_1 = BenefitValues(ChoiceValues>=1);
[BenefitValues_1, ix] = sort(BenefitValues_1,'ascend');
A_projVOI_1 = A_projVOI_1(:,ix);

clr = prism(length(unique(sort(BenefitValues, 'ascend')')));
uniqueBenefit = unique(sort(BenefitValues, 'ascend')');
%%
CB_offerColor= [1.0000    0.2118    0.1373
    1.0000    0.5059         0
    1.0000    0.6235    0.0353
    1.0000    0.7647    0.3020
    0.3216    0.8588    1.0000
    0.0745    0.7098    0.9529
    0.4000    0.6000    1.0000
         0    0.2745    0.4980];
     CB_offerColor = flipud(CB_offerColor([1 3 5 7 8], :));

%% Benefit
yTicks = [min(A_projVOI(:)) max(A_projVOI(:))];
figure;
hold on
plot( [0 0], [yTicks(1) yTicks(end)],'k-', 'LineWidth', 2)
j = 1;
Benefit = [];
k = [];
for i = uniqueBenefit
    if ~isempty(A_projVOI_0(:,BenefitValues_0==i))
        %plot(all_times, A_projVOI_0(:,BenefitValues_0==i), 'lineWidth', 1, 'color',clr(j,:))
        plot(all_times, A_projVOI_0(:,BenefitValues_0==i),':','lineWidth', 2, 'color',CB_offerColor(j, :))
    end
    if ~isempty(A_projVOI_1(:,BenefitValues_1==i))
        %plot(all_times, A_projVOI_1(:,BenefitValues_1==i), 'lineWidth', 2, 'color',clr(j,:))
        plot(all_times, A_projVOI_1(:,BenefitValues_1==i), 'lineWidth', 2, 'color',CB_offerColor(j, :))
    end
    %h = findobj('Color',clr(j,:));
    h = findobj('Color',CB_offerColor(j, :));
    k(j) = h(1);
    Benefit{j} = ['Benefit = ' num2str(uniqueBenefit(j))];
    j = j+1;

end
set(gca, 'xtick', [all_times(1)-0.05 0 all_times(end)+0.05])
xlim([all_times(1)-0.05 all_times(end)+0.05])
ylim(yTicks)

xlabel('time (s)')
ylabel(ActivityType)
legend(k(end:-1:1), Benefit(end:-1:1))

set(gca, 'fontsize',16)
legend boxoff
% %% condition separation
% if ~isempty(strfind(lower(ActivityType), 'benefit'))
%     condSep_t = condSeparation(A_projVOI(:), BenefitValues);
% elseif ~isempty(strfind(lower(ActivityType), 'choice'))
%         condSep_t = condSeparation(A_projVOI(:), ChoiceValues);
% elseif ~isempty(strfind(lower(ActivityType), 'reward'))
%         condSep_t = condSeparation(A_projVOI(:), BenefitValues.*ChoiceValues);
% end
% yTicks = [min(condSep_t(:)) max(condSep_t(:))];
% figure;
% hold on
% plot(all_times, condSep_t, 'k-', 'linewidth',2)
% plot( [0 0], [yTicks(1) yTicks(end)],'y-.', 'LineWidth', 3)
% xlim([all_times(1) all_times(end)])
% ylim(yTicks)
% xlabel('time (s)')
% ylabel([ActivityType ' Relevant Var.'])
% [mxSep, ixmx]= max(condSep_t);
% [mnSep, ixmn]= min(condSep_t);
% plot(all_times(ixmx), mxSep, 'r*');
% plot(all_times(ixmn), mnSep, 'b*');

end