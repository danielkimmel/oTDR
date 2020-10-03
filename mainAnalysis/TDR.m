% TODO: update docstring -- currently outdated and should be ignored by
% public users. 
%
% TDR Input:
% 
% Data
%
% codedParams -- C x P matrix where each row is a unique condition and each
%       column is a task variable (i.e. regression predictor) for which a
%       dRA will be computed. Number of conditions C must match number of
%       elements in Data (see above). (NOTE: Predictor scaling to [0, 1]
%       occurs within function.)
% 
% numPCs -- number of PCs into which to project the neural data prior to
% finding dRAs (functions as a denoising steps). Leave empty [] to include
% all PCs, i.e., no denoising.
% 
% regressBins
% 
% varAnalysisFlg
% 
% stabilityAnalysisFlg
% 
% Optional parameters (can be passed as name,value pairs):
% loocvFlg = true; % logical on whether to use ridge regression in solving for betas
% numSamples4Surrogates = 1000;
% bitPlot -- logical on whether to plot output figures. Default = False.
%
% TDR Output:
% -Summary
%   .dRA: tensor of size TxNx(K+1) (T: time points, N: number of neurons, K number of coded params)
%       .dRA(:,:,1): Benefit normalized regression vectors at all times.
%       .dRA(:,:,2): Choice normalized regression vectors at all times.
%       .dRA(:,:,3): Expected reward normalized regression vectors at all times.
%       .dRA(:,:,4): Baseline normalized regression vectors at all times.
%   .Times:
%       .all_times: times w r t offer of the data
%   .varAnalysis: variance significance analysis
%       .VChance_t: variance across conditions at all times of projections on chance vectors of size (all_times x num chance vectors)
%       .RSVChance1: relevant benefit variance across conditions at all times of projections on chance vectors of size (all_times x num chance vectors)
%       .RSVChance2: relevant choice variance across conditions at all times of projections on chance vectors of size (all_times x num chance vectors)
%       .RSVChance3: relevant expected reward variance across conditions at all times of projections on chance vectors of size (all_times x num chance vectors)
%       .VChance_t: variance across conditions at all times of projections on chance vectors of size (all_times x num chance vectors)
%       .ISVChance1: irrelevant benefit variance across conditions at all times of projections on chance vectors of size (all_times x num chance vectors)
%       .ISVChance2: irrelevant choice variance across conditions at all times of projections on chance vectors of size (all_times x num chance vectors)
%       .ISVChance3: irrelevant expected reward variance across conditions at all times of projections on chance vectors of size (all_times x num chance vectors)
%       .totalVar_t: total variance of all neurons at each time point
%       .totalVarMedD_t: total variance of medD space at each time point
%       .V_dRA1: variance across conditions at all times of projections on dRA1 vectors of size (all_times x all_times num vectors)
%       .V_dRA2: variance across conditions at all times of projections on dRA2 vectors of size (all_times x all_times num vectors)
%       .V_dRA3: variance across conditions at all times of projections on dRA3 vectors of size (all_times x all_times num vectors)
%       .VE_dRA1: percent variance explained across conditions at all times of projections on dRA1 vectors of size (all_times x all_times num vectors)
%       .VE_dRA2: percent variance explained across conditions at all times of projections on dRA2 vectors of size (all_times x all_times num vectors)
%       .VE_dRA3: percent variance explained across conditions at all times of projections on dRA3 vectors of size (all_times x all_times num vectors)
%       .RSV_dRA1: relevant benefit variance across conditions at all times of projections on dRA1 vectors of size (all_times x all_times num vectors)
%       .RSV_dRA2: relevant choice variance across conditions at all times of projections on dRA2 vectors of size (all_times x all_times num vectors)
%       .RSV_dRA3:  relevant expected reward variance across conditions at all times of projections on dRA3 vectors of size (all_times x all_times num vectors)
%       .RSVE_dRA1: percent relevant benefit variance across conditions at all times of projections on dRA1 vectors of size (all_times x all_times num vectors)
%       .RSVE_dRA2: percent relevant choice variance across conditions at all times of projections on dRA2 vectors of size (all_times x all_times num vectors)
%       .RSVE_dRA3: percent relevant expected reward variance across conditions at all times of projections on dRA3 vectors of size (all_times x all_times num vectors)
%       .ISV_dRA1: irrelevant benefit variance across conditions at all times of projections on dRA1 vectors of size (all_times x all_times num vectors)
%       .ISV_dRA2: irrelevant choice variance across conditions at all times of projections on dRA2 vectors of size (all_times x all_times num vectors)
%       .ISV_dRA3:  irrelevant expected reward variance across conditions at all times of projections on dRA3 vectors of size (all_times x all_times num vectors)
%       .ISVE_dRA1: percent irrelevant benefit variance across conditions at all times of projections on dRA1 vectors of size (all_times x all_times num vectors)
%       .ISVE_dRA2: percent irrelevant choice variance across conditions at all times of projections on dRA2 vectors of size (all_times x all_times num vectors)
%       .ISVE_dRA3: percent irrelevant expected reward variance across conditions at all times of projections on dRA3 vectors of size (all_times x all_times num vectors)
%       .zV_dRA1: z-scored deviation from chance of variance across conditions at all times of projections on dRA1 vectors of size (all_times x all_times num vectors)
%       .zV_dRA2: z-scored deviation from chance of variance across conditions at all times of projections on dRA2 vectors of size (all_times x all_times num vectors)
%       .zV_dRA3: z-scored deviation from chance of variance across conditions at all times of projections on dRA3 vectors of size (all_times x all_times num vectors)
%       .zRSV_dRA1: z-scored deviation from chance of relevant benefit variance across conditions at all times of projections on dRA1 vectors of size (all_times x all_times num vectors)
%       .zRSV_dRA2: z-scored deviation from chance of relevant choice variance across conditions at all times of projections on dRA2 vectors of size (all_times x all_times num vectors)
%       .zRSV_dRA3:  z-scored deviation from chance of relevant expected reward variance across conditions at all times of projections on dRA3 vectors of size (all_times x all_times num vectors)
%       .zISV_dRA1: z-scored deviation from chance of irrelevant benefit variance across conditions at all times of projections on dRA1 vectors of size (all_times x all_times num vectors)
%       .zISV_dRA2: z-scored deviation from chance of irrelevant choice variance across conditions at all times of projections on dRA2 vectors of size (all_times x all_times num vectors)
%       .zISV_dRA3:  z-scored deviation from chance of irrelevant expected reward variance across conditions at all times of projections on dRA3 vectors of size (all_times x all_times num vectors)
%   .angleAnalysis: angle between vectors analysis
%       .thetaChance: angles null distribution between any two vector living in the data space
%       .regressVectTimes: regression vectors times
%       .angle11: angle between regression vectors (in med D space) of Benefit and Benefit
%       .angle12: angle between regression vectors (in med D space) of Benefit and Choice
%       .angle13: angle between regression vectors (in med D space) ofBenefit and Expected Reward
%       .angle22: angle between regression vectors (in med D space) of Choice and Choice
%       .angle23: angle between regression vectors (in med D space) of Choice and Expected Reward
%       .angle33: angle between regression vectors (in med D space) of Expected Reward and Expected Reward
%       .pVal11: p value of the significance test (ie prob that the angle is higher than expected by chance) angle between regression vectors (in med D space) of Benefit and Benefit
%       .pVal12: p value of the significance test (ie prob that the angle is higher than expected by chance) angle between regression vectors (in med D space) of Benefit and Choice
%       .pVal13: p value of the significance test (ie prob that the angle is higher than expected by chance) angle between regression vectors (in med D space) ofBenefit and Expected Reward
%       .pVal22: p value of the significance test (ie prob that the angle is higher than expected by chance) angle between regression vectors (in med D space) of Choice and Choice
%       .pVal23: p value of the significance test (ie prob that the angle is higher than expected by chance) angle between regression vectors (in med D space) of Choice and Expected Reward
%       .pVal33: p value of the significance test (ie prob that the angle is higher than expected by chance) angle between regression vectors (in med D space) of Expected Reward and Expected Reward


function [Summary] = TDR(Data, codedParams, numPCs, regressBins, varAnalysisFlg, stabilityAnalysisFlg, varargin)

%% Default param values
% The following parameters can be passed as name,value pairs 

loocvFlg = true; % logical on whether to use ridge regression in solving for betas
numSamples4Surrogates = 1000;

% bitPlot -- logical on whether to plot output figures
bitPlot = false;

%% collect optionally provided parameters
warnopts(assignopts(who, varargin));

%% preprocessing 

% check that number of conditions defined in Data and codedParams match:
if length(Data) ~= size(codedParams,1)
    error('Number of conditions in Data and codedParams must match')
end

% extract number of trials per condition
trCountmtx = vertcat(Data.trialRept)';

% scale and center codedParams to range [0, 1]
codedParams = bsxfun(@minus, codedParams, min(codedParams));
codedParams = bsxfun(@times, codedParams, 1./range(codedParams));

% existing time bins
all_times = Data(1).times;

% extract data dimensions
C = length(Data);
T = length(all_times);
N = size(Data(1).A, 2);

% if number of PCs not provided, assume all PCs
if isempty(numPCs)
    numPCs = N;
end

% preprocess data (mean subtraction, etc)
processedData = tensor2Data(preprocess4TDR(data2Tensor(Data)), rmfield(Data, 'A'));


%% Determine Time bins
% in this section, we collapse sets of N timebins as given by regressBins.
% originally, Gamal had written this to combining bins from the beginning
% of the sequence, so for a series of bins: -1.5 -1.0 -0.5 0.5 1.0 1.5,
% the code would combine for N=2 into: [-1.5 -1], [-0.5 0.5], [1.0 1.5]
% However, we want to treat t=0 as a special case, such that the bins
% respect this boundary, and thus we build out the sets of bins outward
% from t=0. Under this system, the bins would be: 
% [-1.5], [-1 0.5], [0.5 1.0], [1.5]
% 
% In some cases, there will not be an integer multiple of bins to combine
% and thus we will have to eliminate some bins from the beginning or end of
% the timeseries for the purposes of computing the dRAs. Note that we never
% eliminate these excess time bins from the original dataset, thus the
% analyses below that depend on projecting the original data onto the dRAs
% (e.g., projection and varAnalysis) include these original time bins.
%
% Note that we never remove any 

% convert input DATA into tensor (T x N x C)
dataTensor_allTimes = permute(reshape(vertcat(processedData.A)', N, T, C), [2 1 3]);

% determine time bins prior to zero
times_pre = all_times(all_times < 0); 
if regressBins > length(times_pre)
    error('Requested number of bins to combine for regression (%d) exceeds number of times bins prior to t=0 (%d)',regressBins,length(times_pre))
end
% build logical vector that eliminates times from beginning of vector that exceed an integer multiple
% of the number of bins to combine.
times_bitIncl = ismember(all_times,times_pre(end-(floor(length(times_pre)/regressBins)*regressBins-1):end));

% determine time bins on or after zero
times_post = all_times(all_times >= 0); 
if regressBins > length(times_post)
    error('Requested number of bins to combine for regression (%d) exceeds number of times bins at or after t=0 (%d)',regressBins,length(times_post))
end
% extend logical vector to eliminates times from end of vector that exceed an integer multiple
% of the number of bins to combine.
times_bitIncl = times_bitIncl | ismember(all_times,times_post(1:floor(length(times_post)/regressBins)*regressBins));

% eliminate excess times from dataTensor4TDR 
dataTensor4TDR = dataTensor_allTimes(times_bitIncl,:,:);
dataTensor4TDR = squeeze(mean(reshape(dataTensor4TDR, regressBins, sum(times_bitIncl)/regressBins, N, C), 1));

% build vector of times corresponding to new combined bins
regressTimes = all_times(times_bitIncl);
regressTimes = round(mean(reshape(regressTimes, regressBins, sum(times_bitIncl)/regressBins),1), 2);

% GAMAL's OLD WAY -- which allows bins to span t=0
% dataTensor4TDR = dataTensor_allTimes(1:regressBins*floor(T/regressBins), : , :);
% dataTensor4TDR = squeeze(mean(reshape(dataTensor4TDR, regressBins, floor(T/regressBins), N, C), 1));
% regressTimes = all_times(1:floor(T/regressBins)*regressBins);
% regressTimes = round(mean(reshape(regressTimes, regressBins, floor(T/regressBins)),1), 2);


%% Run TDR
% DK -- include constant, aka, bias, term in model:
[dRAs, normdRAs, Summary.runTDRSummary] = runTDR(dataTensor4TDR, numPCs, [codedParams, ones(C,1)], trCountmtx, loocvFlg);
% [dRAs, normdRAs, Summary.runTDRSummary] = runTDR(dataTensor4TDR, numPCs, codedParams, trCountmtx, loocvFlg);

Summary.dRAs = dRAs;
Summary.normdRAs = normdRAs;

%% project data from time t onto dRAs at time t

nRA = size(dRAs,3); % number of regression axes per time
proj = NaN(T,C,nRA);

% find corresponding regression time -- note that any excess time bins
% before of after the series of regression bins (see above) are assigned to
% the nearest regression bin.
posRT = interp1(regressTimes,[1:length(regressTimes)],all_times,'nearest','extrap');

% loop through data times
for t = 1:T
    % project data from t_th time bin onto rt_th dRA, where the t_th time
    % bin corresponds to the rt_th regression time bin.
    proj(t,:,:) = projData(dataTensor_allTimes(t,:,:), squeeze(dRAs(posRT(t),:,:)));
end

% store
Summary.proj = proj;
clear proj

%% time-varying PC
% By Daniel Kimmel, 2016 Dec 09 -- here we compute the first PC in each
% regression bin as a measure of the dimension which would capture the
% maximal variance at each time bin.

% initialize (times x neurons)
dPC = NaN(size(dataTensor4TDR,1),size(dataTensor4TDR,2));

% loop through regression time bins
for t = 1:size(dataTensor4TDR,1)
    % compute PCA across conditions (rows) and neurons (columns)
    [foo] = pca(squeeze(dataTensor4TDR(t,:,:))');
    % store just first PC
    dPC(t,:) = foo(:,1);
end

% save output
Summary.dPC = dPC;

%% variance analysis
%% results here are reshaped to be with respect to times then with respect to components
if varAnalysisFlg
    Summary.varAnalysis = getvarAnalysisFields(TDRvarAnalysis(reshape(permute(dRAs(:,:,1:3), [2 1 3]), N, []), data2Tensor(processedData), codedParams),...
        all_times, regressTimes, 3, 3, bitPlot);
end

%% variance analysis for time-varying PC (DK, 2016 Dec 09)
%% results here are reshaped to be with respect to times then with respect to components
if varAnalysisFlg
    Summary.varAnalysis_dPC = getvarAnalysisFields(TDRvarAnalysis(reshape(permute(dPC, [2 1]), N, []), data2Tensor(processedData), codedParams),...
        all_times, regressTimes, 1, 3, false);
end

%% stability analysis
if stabilityAnalysisFlg
Summary.angleAnalysis = TDRstabilityAnalysis(dataTensor4TDR, dRAs, codedParams, regressTimes, loocvFlg, numPCs, trCountmtx,...
    'numSamples',numSamples4Surrogates);
end
%%
Summary.loocvFlg = loocvFlg;
Summary.Times.all_times = round(all_times,2);
Summary.Times.regressTimes = regressTimes;
Summary.Times.regressBins = regressBins;

Summary.numPCs = numPCs;

end

function Summary = getvarAnalysisFields(varAnalysis, all_times, regressTimes, nDim, nSignal, bitPlot)
T = length(all_times);
rT = length(regressTimes);

V_RA = reshape(varAnalysis.V_RA, T, rT, nDim);
RSV_RA = reshape(varAnalysis.RSV_RA, T, rT, nDim, nSignal);	
ISV_RA = reshape(varAnalysis.ISV_RA, T, rT, nDim, nSignal);

VE_RA = reshape(varAnalysis.VE_RA, T, rT, nDim);
RSVE_RA = reshape(varAnalysis.RSVE_RA, T, rT, nDim, nSignal);	
ISVE_RA = reshape(varAnalysis.ISVE_RA, T, rT, nDim, nSignal);

pV_RA = reshape(varAnalysis.pV_RA, T, rT, nDim);
pRSV_RA = reshape(varAnalysis.pRSV_RA, T, rT, nDim, nSignal);	
pISV_RA = reshape(varAnalysis.pISV_RA, T, rT, nDim, nSignal);

%% save results to summary
Summary.V_RA = V_RA;
Summary.RSV_RA = RSV_RA;
Summary.ISV_RA = ISV_RA;
Summary.VE_RA = VE_RA;
Summary.RSVE_RA = RSVE_RA;
Summary.ISVE_RA = ISVE_RA;
Summary.pV_RA = pV_RA;
Summary.pRSV_RA = pRSV_RA;
Summary.pISV_RA = pISV_RA;

if ~bitPlot
    return
end 
   
%% Plotting

%%%%%%%%%%%%%% Variance explained
minColorBar = 0;
maxColorBar = max(VE_RA(:));
figure;
hp = pcolor(all_times, regressTimes, VE_RA(:,:,1)');
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
hold on
plot([0 0], [all_times(end) all_times(1)],'k-.', 'linewidth',2)
caxis([minColorBar maxColorBar])
colormap(jet)
colorbar
set(gca,'FontSize', 10)
ylabel('Benefit RA at time (s)')
xlabel('time (s)')
title('Variance explained')


figure;
hp = pcolor(all_times, regressTimes, VE_RA(:,:,2)');
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
hold on
plot([0 0], [all_times(end) all_times(1)],'k-.', 'linewidth',2)
caxis([minColorBar maxColorBar])
colormap(jet)
colorbar
set(gca,'FontSize', 10)
ylabel('Choice RA at time (s)')
xlabel('time (s)')
title('Variance explained')

figure;
hp = pcolor(all_times, regressTimes, VE_RA(:,:,3)');
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
hold on
plot([0 0], [all_times(end) all_times(1)],'k-.', 'linewidth',2)
caxis([minColorBar maxColorBar])
colormap(jet)
colorbar
set(gca,'FontSize', 10)
ylabel('Expected Reward RA at time (s)')
xlabel('time (s)')
title('Variance explained')

%%
    %%%%%%%%%%%%%% relevant signal var Variance explained

minColorBar = 0;
maxColorBar = max(RSVE_RA(:));
figure;
hp = pcolor(all_times, regressTimes, RSVE_RA(:,:,1,1)');
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
hold on
plot([0 0], [all_times(end) all_times(1)],'k-.', 'linewidth',2)
caxis([minColorBar maxColorBar])
colormap(jet)
colorbar
set(gca,'FontSize', 10)
ylabel('Benefit RA at time (s)')
xlabel('time (s)')
title('Relevant signal variance explained')

figure;
hp = pcolor(all_times, regressTimes, RSVE_RA(:,:,2,2)');
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
hold on
plot([0 0], [all_times(end) all_times(1)],'k-.', 'linewidth',2)
caxis([minColorBar maxColorBar])
colormap(jet)
colorbar
set(gca,'FontSize', 10)
ylabel('Choice RA at time (s)')
xlabel('time (s)')
title('Relevant signal variance explained')

figure;
hp = pcolor(all_times, regressTimes, RSVE_RA(:,:,3,3)');
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
hold on
plot([0 0], [all_times(end) all_times(1)],'k-.', 'linewidth',2)
caxis([minColorBar maxColorBar])
colormap(jet)
colorbar
set(gca,'FontSize', 10)
ylabel('Expected Reward RA at time (s)')
xlabel('time (s)')
title('Relevant signal variance explained')



%%
    %%%%%%%%%%%%%%%%% p value variance
minColorBar = 0;
maxColorBar = 1;
figure;
hp = pcolor(all_times, regressTimes, pV_RA(:,:,1)');
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
hold on
plot([0 0], [all_times(end) all_times(1)],'k-.', 'linewidth',2)
caxis([minColorBar maxColorBar])
cmp = colormap('jet');
colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
colorbar
% shading interp
set(gca,'FontSize', 10)
ylabel('Benefit RA at time (s)')
xlabel('time (s)')
title(' P value variance')

figure;
hp = pcolor(all_times, regressTimes, pV_RA(:,:,2)');
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
hold on
plot([0 0], [all_times(end) all_times(1)],'k-.', 'linewidth',2)
caxis([minColorBar maxColorBar])
cmp = colormap('jet');
colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
colorbar
% shading interp
set(gca,'FontSize', 10)
ylabel('Choice RA at time (s)')
xlabel('time (s)')
title(' P value variance')

figure;
hp = pcolor(all_times, regressTimes, pV_RA(:,:,3)');
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
hold on
plot([0 0], [all_times(end) all_times(1)],'k-.', 'linewidth',2)
caxis([minColorBar maxColorBar])
cmp = colormap('jet');
colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
colorbar
% shading interp
set(gca,'FontSize', 10)
ylabel('Expected Reward RA at time (s)')
xlabel('time (s)')
title(' P value variance')

    %%%%%%%%%%%%%%%%% p value relevant variance
minColorBar = 0;
maxColorBar = 1;
figure;
hp = pcolor(all_times, regressTimes, pRSV_RA(:,:,1,1)');
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
hold on
plot([0 0], [all_times(end) all_times(1)],'k-.', 'linewidth',2)
caxis([minColorBar maxColorBar])
cmp = colormap('jet');
colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
colorbar
% shading interp
set(gca,'FontSize', 10)
ylabel('Benefit RA at time (s)')
xlabel('time (s)')
title(' P value relevant signal variance')

figure;
hp = pcolor(all_times, regressTimes, pRSV_RA(:,:,2,2)');
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
hold on
plot([0 0], [all_times(end) all_times(1)],'k-.', 'linewidth',2)
caxis([minColorBar maxColorBar])
cmp = colormap('jet');
colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
colorbar
% shading interp
set(gca,'FontSize', 10)
ylabel('Choice RA at time (s)')
xlabel('time (s)')
title(' P value relevant signal variance')

figure;
hp = pcolor(all_times, regressTimes, pRSV_RA(:,:,3,3)');
set(gca,'Ydir','reverse')
set(hp,'edgecolor','none')
% shading interp
hold on
plot([0 0], [all_times(end) all_times(1)],'k-.', 'linewidth',2)
caxis([minColorBar maxColorBar])
cmp = colormap('jet');
colormap(resamplecmap(flipud(cmp), size(cmp,1), logspace(0,log10(size(cmp,1)),size(cmp,1))));
colorbar
% shading interp
set(gca,'FontSize', 10)
ylabel('Expected Reward RA at time (s)')
xlabel('time (s)')
title(' P value relevant signal variance')
end



