function Summary = TDRstabilityAnalysis(dataTensor, dRAs, codedParams, regressTimes, loocvFlg, numPCs, trCountmtx, varargin)
% TODO: update docstring 
% TODO: code is hardcoded to expect 3 task variables
% (for instance, it defines Matlab variables like angle11, angle12,
% angle13, etc). We need to update this function to accommodate an
% arbitrary number of task variables and consequently and arbitrary number
% of pairwise angles. As an example of what the solution may look like,
% check-out TDRvarAnalysis(), which solves this problem for the variance
% analysis (where we had P regression axes, each of which had to be
% evaluated for variance explained WRT each of the P task variables). So
% for the angle analysis, we could do something similar, generating a
% tensor of pairwise angles that was T x T x P x P for T time points and P
% task variables. However, there would be many redundant elements. Another
% approach, would be to store several T x T matrices (either stacked into a
% tensor, or stored in separate elements of a struct), where each matrix
% would define the angels between task variables P_x and P_y. Given P task
% variables, there would be P-choose-2 such matrices. We?d then keep a
% separate key linking a given matrix to task vars P_x and P_y. This
% arrangement would be more akin to how we currently do the angle analysis,
% only extensible.
%% default param values
figFlg = false;
numSamples = 1000;

%% collect optionally provided parameters
warnopts(assignopts(who, varargin));

%%
[T, N, C] = size(dataTensor);
nRA = size(dRAs,3);

%% optimize maximum entropy distribution
[targetSigma_T, targetSigma_N, targetSigma_C, M] = extractFeatures(dataTensor);
covConst{1} = targetSigma_T;
covConst{2} = targetSigma_N;
covConst{3} = eye(C)*trace(targetSigma_C)/C;
[~, mxEntropySummary] = sampleMxEntropyTensor(covConst, 1);
%% measure the surrogates angle threshold
    
surrdRA1 = nan(T, N, numSamples);
surrdRA2 = nan(T, N, numSamples);
surrdRA3 = nan(T, N, numSamples);

surrAngle11 = nan(T, T, numSamples);
surrAngle12 = nan(T, T, numSamples);
surrAngle13 = nan(T, T, numSamples);
surrAngle22 = nan(T, T, numSamples);
surrAngle23 = nan(T, T, numSamples);
surrAngle33 = nan(T, T, numSamples);

% strange new parfor behavior (likely result of MacOS or MATLAB update).
% The parfor loop requires that surrdRA4 be definied even if nRA <= 4. So
% we define it here and then will clear it after the loop
surrdRA4 = nan(T, N, numSamples);

% only define these vars if nRA > 4
if nRA > 4    
    surrAngle14 = nan(T, T, numSamples);
    surrAngle24 = nan(T, T, numSamples);
    surrAngle34 = nan(T, T, numSamples);
    surrAngle44 = nan(T, T, numSamples);
end

fprintf('GENERATING %d SURROGATE DATASETS FOR STABILITY ANALYSIS...\n',numSamples);
parfor i = 1:numSamples
    tic
    [surrTensor] = sampleMxEntropyTensor(covConst, 1, mxEntropySummary.bigCovEigValues); 
    surrTensor = reshape(surrTensor, T, N, C)+M.TN; 
    
    % We do not subtract mean from surrogates (even though the real data
    % are mean-subtracted) because the surrogates are based on
    % mean-subtracted data. This this feature is already built into the
    % statistics of the surrogates, it does not make sense to apply it
    % again.
%     surrTensor = surrTensor-repmat(mean(surrTensor,3), 1, 1, C); % subtract mean across conditions
    warning('OFF','TDR:ProjectOntoNonorthogAxes') % turn off warning about non-orthogonal axes
    [surrdRAs] = runTDR(surrTensor, numPCs, [codedParams, ones(C,1)], trCountmtx, loocvFlg); % include bias term
%     [surrdRAs] = runTDR(surrTensor, numPCs, codedParams, trCountmtx, loocvFlg); % exclude bias term
    warning('ON','TDR:ProjectOntoNonorthogAxes') % turn warning back on
     
    surrdRA1(:,:,i) = surrdRAs(:, :, 1);
    surrdRA2(:,:,i) = surrdRAs(:, :, 2);
    surrdRA3(:,:,i) = surrdRAs(:, :, 3);
    if nRA > 4
        surrdRA4(:,:,i) = surrdRAs(:, :, 4);
    end
    
    [surrAngle11(:,:,i)] = getAllAngles(surrdRA1(:,:,i)', surrdRA1(:,:,i)');
    [surrAngle12(:,:,i)] = getAllAngles(surrdRA1(:,:,i)', surrdRA2(:,:,i)');
    [surrAngle13(:,:,i)] = getAllAngles(surrdRA1(:,:,i)', surrdRA3(:,:,i)');
    [surrAngle22(:,:,i)] = getAllAngles(surrdRA2(:,:,i)', surrdRA2(:,:,i)');
    [surrAngle23(:,:,i)] = getAllAngles(surrdRA2(:,:,i)', surrdRA3(:,:,i)');
    [surrAngle33(:,:,i)] = getAllAngles(surrdRA3(:,:,i)', surrdRA3(:,:,i)');
    
    if nRA > 4
        [surrAngle14(:,:,i)] = getAllAngles(surrdRA1(:,:,i)', surrdRA4(:,:,i)');
        [surrAngle24(:,:,i)] = getAllAngles(surrdRA2(:,:,i)', surrdRA4(:,:,i)');
        [surrAngle34(:,:,i)] = getAllAngles(surrdRA3(:,:,i)', surrdRA4(:,:,i)');
        [surrAngle44(:,:,i)] = getAllAngles(surrdRA4(:,:,i)', surrdRA4(:,:,i)');
    end    
    toc

end

% clean up surrdRA4 if only needed to make parfor happy (see above)
if nRA <= 4
    clear surrdRA4
end

%% measure the significance of data
[angle11] = getAllAngles(dRAs(:,:,1)', dRAs(:,:,1)');
[angle12] = getAllAngles(dRAs(:,:,1)', dRAs(:,:,2)');
[angle13] = getAllAngles(dRAs(:,:,1)', dRAs(:,:,3)');
[angle22] = getAllAngles(dRAs(:,:,2)', dRAs(:,:,2)');
[angle23] = getAllAngles(dRAs(:,:,2)', dRAs(:,:,3)');
[angle33] = getAllAngles(dRAs(:,:,3)', dRAs(:,:,3)');

[pVal11] = sigAngles(angle11, surrAngle11);
[pVal12] = sigAngles(angle12, surrAngle12);
[pVal13] = sigAngles(angle13, surrAngle13);
[pVal22] = sigAngles(angle22, surrAngle22);
[pVal23] = sigAngles(angle23, surrAngle23);
[pVal33] = sigAngles(angle33, surrAngle33);

if nRA > 4
    [angle14] = getAllAngles(dRAs(:,:,1)', dRAs(:,:,4)');
    [angle24] = getAllAngles(dRAs(:,:,2)', dRAs(:,:,4)');
    [angle34] = getAllAngles(dRAs(:,:,3)', dRAs(:,:,4)');
    [angle44] = getAllAngles(dRAs(:,:,4)', dRAs(:,:,4)');

    [pVal14] = sigAngles(angle14, surrAngle14);
    [pVal24] = sigAngles(angle24, surrAngle24);
    [pVal34] = sigAngles(angle34, surrAngle34);
    [pVal44] = sigAngles(angle44, surrAngle44);
end

Summary.angle11 = angle11;
Summary.angle12 = angle12;
Summary.angle13 = angle13;
Summary.angle22 = angle22;
Summary.angle23 = angle23;
Summary.angle33 = angle33;

Summary.pVal11 = pVal11;
Summary.pVal12 = pVal12;
Summary.pVal13 = pVal13;
Summary.pVal22 = pVal22;
Summary.pVal23 = pVal23;
Summary.pVal33 = pVal33;
Summary.surrAngle11 = surrAngle11;
Summary.surrAngle12 = surrAngle12;
Summary.surrAngle13 = surrAngle13;
Summary.surrAngle22 = surrAngle22;
Summary.surrAngle23 = surrAngle23;
Summary.surrAngle33 = surrAngle33;

Summary.surrdRA1 = surrdRA1;
Summary.surrdRA2 = surrdRA2;
Summary.surrdRA3 = surrdRA3;

if nRA > 4
    Summary.angle14 = angle14;
    Summary.angle24 = angle24;
    Summary.angle34 = angle34;
    Summary.angle44 = angle44;
    
    Summary.pVal14 = pVal14;
    Summary.pVal24 = pVal24;
    Summary.pVal34 = pVal34;
    Summary.pVal44 = pVal44;

    Summary.surrAngle14 = surrAngle14;
    Summary.surrAngle24 = surrAngle24;
    Summary.surrAngle34 = surrAngle34;
    Summary.surrAngle44 = surrAngle44;
    
    Summary.surrdRA4 = surrdRA4;
end


Summary.regressTimes = regressTimes;
if figFlg;plotStabilityFigures(Summary);end
%%
end


function [angleij] = getAllAngles(dRAsi, dRAsj)
% angleij = real(acos(abs(dRAsi'*dRAsj))*180/pi); %% angle in abs sense 0-90
angleij = real(acos((dRAsi'*dRAsj))*180/pi); %% angle 0-180
end




function [pValsij] = sigAngles(anglesij, surrAnglesij)
pValsij = nan(size(anglesij));
L = length(anglesij);
for i = 1:L
    for j = setdiff(1:L, i)
%         pValsij(i,j) = sigGamma(anglesij(i,j), reshape(surrAnglesij(i,j,:), [], 1), 'lower');
        % DK: replaced below line with more general function for empirical
        % significance testing. Also, apply GOF error threshold (returning
        % NaNs for poor fits), but do not throw and error for poor fits,
        % since pValues here are a rough estimate. We typically re-compute
        % pValues later based on angles that have been reflected about
        % 90deg.
%         pValsij(i,j) = sigGamma(anglesij(i,j), reshape(surrAnglesij(i,j,:), [], 1), 'lower');
        pValsij(i,j) = sigTest(anglesij(i,j),reshape(surrAnglesij(i,j,:), [], 1),...
            'lower','gamma','bitIgnorePoorFit',true);

    end
end
end


function plotStabilityFigures(Summary)
SummaryFields = fieldnames(Summary);
for i = 1:length(SummaryFields)
   eval([SummaryFields{i} '= Summary.' SummaryFields{i} ';']); 
end
angle11= angle11+diag(nan(length(regressTimes),1));
angle22= angle22+diag(nan(length(regressTimes),1));
angle33= angle33+diag(nan(length(regressTimes),1));
scleAngle = [10*floor(min([angle11(:);angle22(:);angle33(:)])./10) 90];
figure;
hp = pcolor(regressTimes, regressTimes, angle11);
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
caxis(scleAngle)
mp = colormap(jet);
colormap(flipud(mp))
colorbar
set(gca,'FontSize', 22)
xlabel('Benefit time s')
ylabel('Benefit time s')
title('\theta')
axis square
set(gca, 'fontsize',16)

figure;
hp = pcolor(regressTimes, regressTimes, angle12);
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
caxis([0 90])
mp = colormap(jet);
colormap(flipud(mp))
colorbar
set(gca,'FontSize', 22)
ylabel('Benefit time s')
xlabel('Choice time s')
title('\theta')
axis square
set(gca, 'fontsize',16)



figure;
hp = pcolor(regressTimes, regressTimes, angle13);
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
caxis([0 90])
mp = colormap(jet);
colormap(flipud(mp))
colorbar
set(gca,'FontSize', 22)
ylabel('Benefit time s')
xlabel('Expected Reward time s')
title('\theta')
axis square
set(gca, 'fontsize',16)


figure;
hp = pcolor(regressTimes, regressTimes, angle22);
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
caxis(scleAngle)
mp = colormap(jet);
colormap(flipud(mp))
colorbar
set(gca,'FontSize', 22)
ylabel('Choice time s')
xlabel('Choice time s')
title('\theta')
axis square
set(gca, 'fontsize',16)


figure;
hp = pcolor(regressTimes, regressTimes, angle23);
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
caxis([0 90])
mp = colormap(jet);
colormap(flipud(mp))
colorbar
set(gca,'FontSize', 23)
ylabel('Choice time s')
xlabel('Expected Reward time s')
title('\theta')
axis square
set(gca, 'fontsize',16)





figure;
hp = pcolor(regressTimes, regressTimes, angle33);
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
caxis(scleAngle)
mp = colormap(jet);
colormap(flipud(mp))
colorbar
set(gca,'FontSize', 22)
ylabel('Expected Reward time s')
xlabel('Expected Reward time s')
title('\theta')
axis square
set(gca, 'fontsize',16)
%%
if regressTimes(end)>6
    scleP = [-log10(0.05) 5 ];

else
    scleP = [-log10(0.05) 4 ];
end
figure;
hp = pcolor(regressTimes, regressTimes, -log10(pVal11));
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
hold on
caxis(scleP)
cmp = colormap('jet');
% colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
colormap(cmp)
hcb = colorbar;
xlabel('Benefit time s')
ylabel('Benefit time s')
title('\theta pVal')
axis square
set(gca, 'fontsize',16)
set(hcb,'YTick',scleP)
set(hcb,'YTickLabel',10.^-scleP)


figure;
hp = pcolor(regressTimes, regressTimes, -log10(pVal12));
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
hold on
caxis(scleP)
cmp = colormap('jet');
% colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
colormap(cmp)
hcb = colorbar;
ylabel('Benefit time s')
xlabel('Choice time s')
title('\theta pVal')
axis square
set(gca, 'fontsize',16)
set(hcb,'YTick',scleP)
set(hcb,'YTickLabel',10.^-scleP)



figure;
hp = pcolor(regressTimes, regressTimes, -log10(pVal13));
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
hold on
caxis(scleP)
cmp = colormap('jet');
% colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
colormap(cmp)
hcb = colorbar;
ylabel('Benefit time s')
xlabel('Expected Reward time s')
title('\theta pVal')
axis square
set(gca, 'fontsize',16)
set(hcb,'YTick',scleP)
set(hcb,'YTickLabel',10.^-scleP)


figure;
hp = pcolor(regressTimes, regressTimes, -log10(pVal22));
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
hold on
caxis(scleP)
cmp = colormap('jet');
% colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
colormap(cmp)
hcb = colorbar;
ylabel('Choice time s')
xlabel('Choice time s')
title('\theta pVal')
axis square
set(gca, 'fontsize',16)
set(hcb,'YTick',scleP)
set(hcb,'YTickLabel',10.^-scleP)


figure;
hp = pcolor(regressTimes, regressTimes, -log10(pVal23));
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
hold on
caxis(scleP)
cmp = colormap('jet');
% colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
colormap(cmp)
hcb = colorbar;
ylabel('Choice time s')
xlabel('Expected Reward time s')
title('\theta pVal')
axis square
set(gca, 'fontsize',16)
set(hcb,'YTick',scleP)
set(hcb,'YTickLabel',10.^-scleP)




figure;
hp = pcolor(regressTimes, regressTimes, -log10(pVal33));
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
hold on
caxis(scleP)
cmp = colormap('jet');
% colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
colormap(cmp)
hcb = colorbar;
ylabel('Expected Reward time s')
xlabel('Expected Reward time s')
title('\theta pVal')
axis square
set(gca, 'fontsize',16)
set(hcb,'YTick',scleP)
set(hcb,'YTickLabel',10.^-scleP)

%%

% % % % % % % % % % % % % figure;
% % % % % % % % % % % % % % hp = pcolor(regressTimes, regressTimes, pVal11);
% % % % % % % % % % % % % hp = pcolor(regressTimes, regressTimes, -log10(pVal11));
% % % % % % % % % % % % % set(gca,'Ydir','reverse')
% % % % % % % % % % % % % set(hp,'edgecolor','none')
% % % % % % % % % % % % % % shading interp
% % % % % % % % % % % % % hold on
% % % % % % % % % % % % % % caxis([0 1])
% % % % % % % % % % % % % caxis([-log10(0.05) 10])
% % % % % % % % % % % % % cmp = colormap('jet');
% % % % % % % % % % % % % % colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
% % % % % % % % % % % % % colormap(cmp);
% % % % % % % % % % % % % hcb = colorbar;
% % % % % % % % % % % % % xlabel('Benefit time s')
% % % % % % % % % % % % % ylabel('Benefit time s')
% % % % % % % % % % % % % title('\theta pVal')
% % % % % % % % % % % % % axis square
% % % % % % % % % % % % % set(gca, 'fontsize',16)
% % % % % % % % % % % % % set(hcb,'YTick',[-log10(0.05) 10])
% % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % figure;
% % % % % % % % % % % % % % hp = pcolor(regressTimes, regressTimes, pVal12);
% % % % % % % % % % % % % hp = pcolor(regressTimes, regressTimes, -log10(pVal12));
% % % % % % % % % % % % % set(gca,'Ydir','reverse')
% % % % % % % % % % % % % set(hp,'edgecolor','none')
% % % % % % % % % % % % % % shading interp
% % % % % % % % % % % % % hold on
% % % % % % % % % % % % % % caxis([0 1])
% % % % % % % % % % % % % caxis([-log10(0.05) 10])
% % % % % % % % % % % % % cmp = colormap('jet');
% % % % % % % % % % % % % % colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
% % % % % % % % % % % % % colormap(cmp);
% % % % % % % % % % % % % hcb = colorbar;
% % % % % % % % % % % % % ylabel('Benefit time s')
% % % % % % % % % % % % % xlabel('Choice time s')
% % % % % % % % % % % % % title('\theta pVal')
% % % % % % % % % % % % % axis square
% % % % % % % % % % % % % set(gca, 'fontsize',16)
% % % % % % % % % % % % % set(hcb,'YTick',[-log10(0.05) 10])
% % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % figure;
% % % % % % % % % % % % % % hp = pcolor(regressTimes, regressTimes, pVal13);
% % % % % % % % % % % % % hp = pcolor(regressTimes, regressTimes, -log10(pVal13));
% % % % % % % % % % % % % set(gca,'Ydir','reverse')
% % % % % % % % % % % % % set(hp,'edgecolor','none')
% % % % % % % % % % % % % % shading interp
% % % % % % % % % % % % % hold on
% % % % % % % % % % % % % % caxis([0 1])
% % % % % % % % % % % % % caxis([-log10(0.05) 10])
% % % % % % % % % % % % % cmp = colormap('jet');
% % % % % % % % % % % % % % colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
% % % % % % % % % % % % % colormap(cmp);
% % % % % % % % % % % % % hcb = colorbar;
% % % % % % % % % % % % % ylabel('Benefit time s')
% % % % % % % % % % % % % xlabel('Expected Reward time s')
% % % % % % % % % % % % % title('\theta pVal')
% % % % % % % % % % % % % axis square
% % % % % % % % % % % % % set(gca, 'fontsize',16)
% % % % % % % % % % % % % set(hcb,'YTick',[-log10(0.05) 10])
% % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % figure;
% % % % % % % % % % % % % % hp = pcolor(regressTimes, regressTimes, pVal22);
% % % % % % % % % % % % % hp = pcolor(regressTimes, regressTimes, -log10(pVal22));
% % % % % % % % % % % % % set(gca,'Ydir','reverse')
% % % % % % % % % % % % % set(hp,'edgecolor','none')
% % % % % % % % % % % % % % shading interp
% % % % % % % % % % % % % hold on
% % % % % % % % % % % % % % caxis([0 1])
% % % % % % % % % % % % % caxis([-log10(0.05) 10])
% % % % % % % % % % % % % cmp = colormap('jet');
% % % % % % % % % % % % % % colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
% % % % % % % % % % % % % colormap(cmp);
% % % % % % % % % % % % % hcb = colorbar;
% % % % % % % % % % % % % ylabel('Choice time s')
% % % % % % % % % % % % % xlabel('Choice time s')
% % % % % % % % % % % % % title('\theta pVal')
% % % % % % % % % % % % % axis square
% % % % % % % % % % % % % set(gca, 'fontsize',16)
% % % % % % % % % % % % % set(hcb,'YTick',[-log10(0.05) 10])
% % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % figure;
% % % % % % % % % % % % % % hp = pcolor(regressTimes, regressTimes, pVal23);
% % % % % % % % % % % % % hp = pcolor(regressTimes, regressTimes, -log10(pVal23));
% % % % % % % % % % % % % set(gca,'Ydir','reverse')
% % % % % % % % % % % % % set(hp,'edgecolor','none')
% % % % % % % % % % % % % % shading interp
% % % % % % % % % % % % % hold on
% % % % % % % % % % % % % % caxis([0 1])
% % % % % % % % % % % % % caxis([-log10(0.05) 10])
% % % % % % % % % % % % % cmp = colormap('jet');
% % % % % % % % % % % % % % colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
% % % % % % % % % % % % % colormap(cmp);
% % % % % % % % % % % % % hcb = colorbar;
% % % % % % % % % % % % % ylabel('Choice time s')
% % % % % % % % % % % % % xlabel('Expected Reward time s')
% % % % % % % % % % % % % title('\theta pVal')
% % % % % % % % % % % % % axis square
% % % % % % % % % % % % % set(gca, 'fontsize',16)
% % % % % % % % % % % % % set(hcb,'YTick',[-log10(0.05) 10])
% % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % figure;
% % % % % % % % % % % % % % hp = pcolor(regressTimes, regressTimes, pVal33);
% % % % % % % % % % % % % hp = pcolor(regressTimes, regressTimes, -log10(pVal33));
% % % % % % % % % % % % % set(gca,'Ydir','reverse')
% % % % % % % % % % % % % set(hp,'edgecolor','none')
% % % % % % % % % % % % % hold on
% % % % % % % % % % % % % % caxis([0 1])
% % % % % % % % % % % % % caxis([-log10(0.05) 10])
% % % % % % % % % % % % % cmp = colormap('jet');
% % % % % % % % % % % % % % colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
% % % % % % % % % % % % % colormap(cmp);
% % % % % % % % % % % % % hcb = colorbar;
% % % % % % % % % % % % % ylabel('Expected Reward time s')
% % % % % % % % % % % % % xlabel('Expected Reward time s')
% % % % % % % % % % % % % title('\theta pVal')
% % % % % % % % % % % % % axis square
% % % % % % % % % % % % % set(gca, 'fontsize',16)
% % % % % % % % % % % % % set(hcb,'YTick',[-log10(0.05) 10])
%%
% % % surrTauPrct11 = prctile(surrTau11, [0.5 0.95], 2);
% % % surrTauPrct22 = prctile(surrTau22, [0.5 0.95], 2);
% % % surrTauPrct33 = prctile(surrTau33, [0.5 0.95], 2);
% % % clr = colorCategorical(3);
% % % scle = mean(diff(regressTimes)); % convert from bins to sec
% % % figure
% % % hold on
% % % area(regressTimes(:), surrTauPrct11(:,2)*scle, 'facecolor', 0.9*[1 1 1], 'linewidth', 2)
% % % plot(regressTimes, tau11*scle, 'color', clr(1,:))
% % % xlim([0 ceil(regressTimes(end))])
% % % xlabel('time (s)')
% % % ylabel('benefit \tau (s)')
% % % set(gca, 'fontsize', 16)
% % % 
% % % figure
% % % hold on
% % % area(regressTimes(:), surrTauPrct22(:,2)*scle, 'facecolor', 0.9*[1 1 1], 'linewidth', 2)
% % % plot(regressTimes, tau22*scle, 'color', clr(2,:))
% % % xlim([0 ceil(regressTimes(end))])
% % % xlabel('time (s)')
% % % ylabel('choice \tau (s)')
% % % set(gca, 'fontsize', 16)
% % % 
% % % figure
% % % hold on
% % % area(regressTimes(:), surrTauPrct33(:,2)*scle, 'facecolor', 0.9*[1 1 1], 'linewidth', 2)
% % % plot(regressTimes, tau33*scle, 'color', clr(3,:))
% % % xlim([0 ceil(regressTimes(end))])
% % % xlabel('time (s)')
% % % ylabel('expected-reward \tau (s)')
% % % set(gca, 'fontsize', 16)
% % % 
% % % 
% % % hf(1) = figure;
% % % set(hf(1), 'color', [1 1 1]);
% % % set(hf(1),'position',[100 100 400 600])
% % % hold on
% % % h11 = bar(1, mean(tau11)*scle);
% % % set(h11(1),'facecolor',[26 173 150]/255,'barwidth',0.6,'edgecolor','none')
% % % errorbar(1, mean(tau11)*scle, 2*stdEr(tau11*scle), 'ko', 'markerfacecolor', 'k')
% % % h22 = bar(2, mean(tau22)*scle);
% % % set(h22(1),'facecolor',[213 124 53]/255,'barwidth',0.6,'edgecolor','none')
% % % errorbar(2, mean(tau22)*scle, 2*stdEr(tau22*scle), 'ko', 'markerfacecolor', 'k')
% % % h33 = bar(3, mean(tau33)*scle);
% % % set(h33(1),'facecolor',[123 136 192]/255,'barwidth',0.6,'edgecolor','none')
% % % errorbar(3, mean(tau33)*scle, 2*stdEr(tau33*scle), 'ko', 'markerfacecolor', 'k')
% % % xlim([0.5 3.5])
% % % set(gca,'xtick',[1:3], 'xtickLabel',{'Benefit', 'Choice', 'EReward'})
% % % set(gca,'FontSize',16)
% % % ylabel('\tau (s)')

end
