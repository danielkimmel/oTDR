%% output:
%    Summary.
%       TDR.
%           DoVE: difference of variance analyses defined as
%           DoVE(ti, tj) = |VE_ti(tj)-VE_tj(tj)| , where VE_ti(tj) is the
%           variance explained of the TDR regression vector (dRAs) defined 
%           at time ti and the variance explained is evaluated at time tj.
%           pVal: pvalue of the DoVE measure (upper tail test).
%       oTDR.
%           DoVE: difference of variance analyses defined as
%           DoVE(ti, tj) = |VE_*(ti)-VE_*(tj)| , where VE_*(ti) is the
%           variance explained of the oTDR regression vector (sRAs) and the 
%           variance explained is evaluated at time ti.
%           pVal: pvalue of the DoVE measure (upper tail test).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Summary] = DOVE(bigA, randVects, regressVects, Vref, all_times, Predictor)
trialL = length(all_times);
numConds = size(bigA,1)./trialL;
numNeus = size(bigA,2);
[~, regressVects] = normVects(regressVects);
[~, Vref] = normVects(Vref);

varNeu_t = nan(trialL, numNeus);
for n = 1:numNeus
   varNeu_t(:,n) = var(reshape(bigA(:,n), trialL, numConds),[],2);
end
totalVar_t = sum(varNeu_t,2);

bigA_regressProj = bsxfun(@minus,bigA, mean(bigA))*regressVects;
bigA_randProj = bsxfun(@minus,bigA, mean(bigA))*randVects;
bigA_refProj = bsxfun(@minus,bigA, mean(bigA))*Vref;

varERegressVect_t = nan(length(all_times),trialL);
for i = 1:length(all_times)
    varERegressVect_t(i,:) = var(reshape(bigA_regressProj(:,i), trialL, numConds),[],2)./totalVar_t.*100;
end

varERand_t = nan(size(randVects,2),trialL);
for i = 1:size(randVects,2)
    varERand_t(i, :) = var(reshape(bigA_randProj(:,i), trialL, numConds),[],2)./totalVar_t.*100;
end

varErefV_t = var(reshape(bigA_refProj, trialL, numConds),[],2)./totalVar_t.*100;
%%
for i = 1:length(all_times)
    DoVE1(i, :) = abs(varERegressVect_t(i, :) - diag(varERegressVect_t(:, :))');
end


for i = 1:size(randVects,2)
    VE1 = varERand_t(randi(size(randVects,2)),:);
    VE2 = varERand_t(randi(size(randVects,2)),:);
    DoVE1_rand(i, :) = abs(VE1-VE2);
end

for i = 1:length(all_times)
    DoVE2(i, :) = abs(varErefV_t-varErefV_t(i));
    DoVE2_rand(:, :,i) = abs(bsxfun(@minus, varERand_t, varERand_t(:,i)));
end


for i = 1:length(all_times)
   for j = 1:length(all_times) 
       pVal1(i,j) = sum(DoVE1(i,j)<= DoVE1_rand(i, :))/length(DoVE1_rand(i, :));
   end
end


for i = 1:length(all_times)
   for j = 1:length(all_times) 
       DovE2randij = DoVE2_rand(:,j, i);
       pVal2(i,j) = sum(DoVE2(i,j)<= DovE2randij(:) )/length(DovE2randij(:));
       if sum(abs(DoVE2_rand(:,j, i)- DoVE2_rand(:,i, j)))>eps
          %keyboard; 
       end
   end
end

Summary.TDR.DoVE = DoVE1;
Summary.TDR.pVal = pVal1;
Summary.oTDR.DoVE = DoVE2;
Summary.oTDR.pVal = pVal2;

% % % % figure;
% % % % heatmap(DoVE1, all_times, all_times, [], 'TickAngle', 45,...
% % % %         'ShowAllTicks', true, 'TickFontSize', 16,'Colormap', jet,'Colorbar',true, 'MinColorValue', 0, 'MaxColorValue', 100,'UseLogColormap', false);
% % % % ylabel('regressVect (dRA_t)')
% % % % xlabel('time (s)')
% % % % title('DoVE (dRA)')
% % % % 
% % % % 
% % % % 
% % % % figure;
% % % % heatmap(pVal1, all_times, all_times, [], 'TickAngle', 45,...
% % % %         'ShowAllTicks', true, 'TickFontSize', 16,'Colormap', flipud(jet),'Colorbar',true, 'MinColorValue', 0, 'MaxColorValue', 0.05,'UseLogColormap', false);
% % % % ylabel('regressVect (dRA_t)')
% % % % xlabel('time (s)')
% % % % title('p-value DoVE (dRA)')
% % % % 
% % % % 
% % % % figure;
% % % % heatmap(DoVE2, all_times, all_times, [], 'TickAngle', 45,...
% % % %         'ShowAllTicks', true, 'TickFontSize', 16,'Colormap', jet,'Colorbar',true, 'MinColorValue', 0, 'MaxColorValue', 100,'UseLogColormap', false);
% % % % ylabel('time (s)')
% % % % xlabel('time (s)')
% % % % title('DoVE (sRA)')
% % % % 
% % % % 
% % % % figure;
% % % % heatmap(pVal2, all_times, all_times, [], 'TickAngle', 45,...
% % % %         'ShowAllTicks', true, 'TickFontSize', 16,'Colormap', flipud(jet),'Colorbar',true, 'MinColorValue', 0, 'MaxColorValue', 0.05,'UseLogColormap', false);
% % % % ylabel('time (s)')
% % % % xlabel('time (s)')
% % % % title('p-value DoVE (sRA)')
% % % % 


%%
Summary = struct;
VEref_t = var(reshape(bigA_refProj, trialL, numConds),[],2)./totalVar_t.*100;
RSVEref_t = condSeparation(bigA_refProj, Predictor)./totalVar_t.*100;
DoVE1 = [];
DoRSVE1 = [];
for i = 1:length(all_times)
    bigA_regressProj = bsxfun(@minus,bigA, mean(bigA))*regressVects(:,i);
    VE = var(reshape(bigA_regressProj, trialL, numConds),[],2)./totalVar_t.*100;
    RSVE = condSeparation(bigA_regressProj, Predictor)./totalVar_t.*100;
    DoVE1(i,1) = abs(VEref_t(i)-VE(i));
    DoRSVE1(i,1) = abs(RSVEref_t(i)-RSVE(i));
end


% DoVE1 = abs(varErefV_t - diag(varERegressVect_t(:, :)));

DoVE1_rand = [];
DoRSVE1_rand = [];
for i = 1:size(randVects,2)
    bigA_randProj1 = bsxfun(@minus,bigA, mean(bigA))*randVects(:, randi(size(randVects,2)));
    bigA_randProj2 = bsxfun(@minus,bigA, mean(bigA))*randVects(:, randi(size(randVects,2)));
    
    VE1 = var(reshape(bigA_randProj1, trialL, numConds),[],2)./totalVar_t.*100;
    VE2 = var(reshape(bigA_randProj2, trialL, numConds),[],2)./totalVar_t.*100;
    DoVE1_rand(:, i) = abs(VE1-VE2);
    
    [RSVE1] = condSeparation(bigA_randProj1, Predictor)./totalVar_t.*100;
    [RSVE2] = condSeparation(bigA_randProj2, Predictor)./totalVar_t.*100;
    
    DoRSVE1_rand(:, i) = abs(RSVE1-RSVE2);
end

sigThreshold = norminv(0.95);

meanDoVE_randProj = mean(DoVE1_rand,2);
stdDoVE_randProj = var(DoVE1_rand, 0 ,2).^0.5;
CIBoundaryDoVE = meanDoVE_randProj+sigThreshold*stdDoVE_randProj;

tmin = all_times(1);
tmax = all_times(end);
figure
yTicks = [0 10*ceil(max(DoVE1)/10)];
hold on
h0 = plot(all_times, meanDoVE_randProj + sigThreshold*stdDoVE_randProj, 'color',[0.5 0.5 0.5], 'LineWidth', 2);
area(all_times',meanDoVE_randProj + sigThreshold*stdDoVE_randProj, 'Facecolor',[0.5 0.5 0.5]);
ylabel('DoVE (%)','FontSize',16)
set(gca,'ytick',yTicks)
set(gca,'ylim',[yTicks(1)-0.5 yTicks(end)])
xlabel('Time (ms)','FontSize',16)
set(gca,'xtick',[tmin 0 tmax])
set(gca,'xlim',[tmin tmax])
set(gca,'FontSize',16)
h2 = plot(all_times, DoVE1, 'color', 'k','linewidth',2);
plot( [0 0], [yTicks(1) yTicks(end)],'k-.', 'LineWidth', 2)


zDoVE1 = bsxfun(@times, bsxfun(@minus, DoVE1,meanDoVE_randProj),1./stdDoVE_randProj);

tmin = all_times(1);
tmax = all_times(end);
yTicks = [0 sigThreshold round(max([zDoVE1(:)])+1)];
figure
hold on;
h0 = fill([all_times(1) all_times(end)  all_times(end) all_times(1)],[-sigThreshold -sigThreshold sigThreshold sigThreshold], [0.5 0.5 0.5]);
h2 = plot(all_times, zDoVE1 , 'color','k', 'LineWidth', 2); 
plot( [0 0], [yTicks(1) yTicks(end)],'k-.', 'LineWidth', 2)
set(gca,'ytick',yTicks)
set(gca,'ylim',[yTicks(1) yTicks(end)])
xlabel('Time (ms)','FontSize',16)
set(gca,'xtick',[tmin 0 tmax])
set(gca,'xlim',[tmin tmax])
ylabel('z-score DoVE (%)','FontSize',16)
set(gca,'FontSize',16)


%%%%%%%%%%%%%%%
meanDoRSVE_randProj = mean(DoRSVE1_rand,2);
stdDoRSVE_randProj = var(DoRSVE1_rand, 0 ,2).^0.5;
CIBoundaryDoRSVE = meanDoRSVE_randProj+sigThreshold*stdDoRSVE_randProj;

tmin = all_times(1);
tmax = all_times(end);
figure
yTicks = [0 10*ceil(max(DoRSVE1)/10)];
hold on
h0 = plot(all_times, meanDoRSVE_randProj + sigThreshold*stdDoRSVE_randProj, 'color',[0.5 0.5 0.5], 'LineWidth', 2);
area(all_times',meanDoRSVE_randProj + sigThreshold*stdDoRSVE_randProj, 'Facecolor',[0.5 0.5 0.5]);
ylabel('DoRSVE (%)','FontSize',16)
set(gca,'ytick',yTicks)
set(gca,'ylim',[yTicks(1)-0.5 yTicks(end)])
xlabel('Time (ms)','FontSize',16)
set(gca,'xtick',[tmin 0 tmax])
set(gca,'xlim',[tmin tmax])
set(gca,'FontSize',16)
h2 = plot(all_times, DoRSVE1, 'color', 'k','linewidth',2);
plot( [0 0], [yTicks(1) yTicks(end)],'k-.', 'LineWidth', 2)


zDoRSVE1 = bsxfun(@times, bsxfun(@minus, DoRSVE1,meanDoRSVE_randProj),1./stdDoRSVE_randProj);

tmin = all_times(1);
tmax = all_times(end);
yTicks = [0 sigThreshold round(max([zDoRSVE1(:)])+1)];
figure
hold on;
h0 = fill([all_times(1) all_times(end)  all_times(end) all_times(1)],[-sigThreshold -sigThreshold sigThreshold sigThreshold], [0.5 0.5 0.5]);
h2 = plot(all_times, zDoRSVE1 , 'color','k', 'LineWidth', 2); 
plot( [0 0], [yTicks(1) yTicks(end)],'k-.', 'LineWidth', 2)
set(gca,'ytick',yTicks)
set(gca,'ylim',[yTicks(1) yTicks(end)])
xlabel('Time (ms)','FontSize',16)
set(gca,'xtick',[tmin 0 tmax])
set(gca,'xlim',[tmin tmax])
ylabel('z-score DoRSVE (%)','FontSize',16)
set(gca,'FontSize',16)

Summary.DoVE = DoVE1;
Summary.DoVE_rand = DoVE1_rand;
Summary.CIBoundaryDoVE = CIBoundaryDoVE;

Summary.DoRSVE = DoRSVE1;
Summary.DoRSVE_rand = DoRSVE1_rand;
Summary.CIBoundaryDoRSVE = CIBoundaryDoRSVE;

end