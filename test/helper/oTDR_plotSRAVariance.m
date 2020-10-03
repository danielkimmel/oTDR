function [figH] = oTDR_plotSRAVariance(Summary,varargin)
% TODO(D): update docstring. 
% TODO: Note this function is specific to CB task
% because of the names of the summary structs (e.g., oTDRSummary,
% oTDRSummary_PP, etc) for which it depends. However, it would probably be
% a relatively straightforward task to revise the function to work more
% generally. Consider this for a future release.
%
% [figH] = oTDR_plotSRAVariance(Summary,[name,value])
% 
% Plots variance explained by sRAs (regression axes) from oTDR 
%
% Accept SUMMARY struct with substructs
%   .oTDRSummary -- from oTDR.m
%   .TDRSummary -- from TDR.m
%
% Accepts optional inputs as name,value pairs (see code).
%
% Returns FIGH -- handle to figure.
%
% Daniel Kimmel, excepted into separate function December 16, 2016

%% default input parameter values

% set which temporal epoch to use by passing integer value:
%   1: present-trial epoch projected onto present-trial sRAs
%   2: previous-trial epoch projected onto present-trial sRAs
%   3: previous-trial epoch projected onto previous-trial sRAs
%   4: present-trial epoch coded WRT PRESENT-trial conditions projected
%       onto previous-trial sRAs.
%   5: present-trial epoch coded WRT PREVIOUS-trialconditions projected
%       onto previous-trial sRAs. 
eN = 1;

% set variance (V), relevant signal variance (RSV), anti-relevant signal
% varaince (ARSV; when r < 0) irrelevant signal variance (ISV). Set to 'R2'
% to plot the R2_tk field from TDRSummary.runTDRSummary, which measures the
% contribution of the signal-specific term to the total estimated firing
% rate. Set to 'normdRA' to plot the norm (magnitude) of the non-unit
% vector dRA(t). This is an estimate of the magnitude of signal at time t.
% Include multiple entries to plot multiple signals (except for R2 and
% normdRA, which can only take 1 signal)
var2Plot =  {'V'}; % {'R2'}; {'RSV','ISV','ARSV'}; {'normdRA'}
 
% set type of variance to plot by passing integer value:
%   1: absolute variance 
%   2: variance explained expressed as percentage relative to total variance
%   3: (DEFUNCT)
%   4: p-value of variance explained based on comparing veridical value to
%       gamma distribution fit to empirical null distribution
varType = 2;

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
bitMatchSignal2SRA = false;

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
% median and 95% timing of key task events, respectively. Note this is
% task-specific and depends on the input struct (Summary) having field
% .oTDRSummary(or oTDRSummary_PP).Times.timeIntervals, which itself is a
% struct returned by getNstructFromSpks(). Fur future release, if one is
% generalizing the current plotting function, timeIntervals could be made
% into a separate input arg and reformatted into a more general
% task-agnostic format.
bitTimeOfEvent = true;

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

% default names of oTDR and TDR substructs
otdr = 'oTDRSummary';
tdr = 'TDRSummary';
otdrSerial = 'oTDRSummary_serial';

% set analysis (tdr vs otdr). To work for TDR, one has to add code to
% select which TDR vector to plot.
fn = 'oTDRSummary';

% set which total variance to normalize by: totalVar_t or totalVarMedD_t
totVarFN = 'totalVar_t';

% chance level
alpha = 0.05;

% appearance
bgColor = 'w';
lineColor = 'k';
textColor = 'k';

% logical on whether to use horizontal line (instead of fill patch) for
% p-value treshold. When true, overrides bit4EPS.
bitLine4PThresh = false;

% logical on whether to change patches to transparent given matlab's
% strange behavior when exporting filled patches. Irrelevant when
% bitLine4PThresh == TRUE.
bit4EPS = false; 

% labels
condLabel = {''};

% name of study (needed for selecting colors). Defaults to 'CB'. Set to []
% for generic study.
studyName = 'CB';

%% collect optionally provided parameters
warnopts(assignopts(who, varargin));

%% checks

%%% force certain parameters
% When using common condition PCs, force the bitUsePC flag
if bitUsePC_CC
    bitUsePC = true;
end

if varType == 3
    error('We no longer compute z-scored variance. All normalization is done by fitting gamma distributions and computing p-values from them.')
end

if any(ismember(var2Plot,{'R2','normdRA'})) && length(var2Plot) > 1
    error('Can only plot single variable at a time when plotting R2 or normdRA')
end 

if any(ismember(var2Plot,{'R2','normdRA'})) && ~bitMatchSignal2SRA
    warning('When plotting R2 or normdRA, bitMatchSignal2SRA=true is implied. Will set bitMatchSignal2SRA=true.')
    bitMatchSignal2SRA = true;
end

if all(ismember(var2Plot,{'V'})) && ~bitMatchSignal2SRA
    warning('Plotting variance related to separate signals (bitMatchSignal2SRA==false) is not well-defined for varType==''V''. Setting to bitMatchSignal2SRA=true');
    bitMatchSignal2SRA = true;
end

if bitSumDim && bitMatchSignal2SRA && ~bitUsePC && any(ismember(var2Plot,{'RSV','ISV','ARSV'}))
    warning('When summing across dimensions (bitSumDim=true), RSV/ISV for each signal is plotted based on the summed variance (implying bitMatchSignal2SRA=false), thus setting bitMatchSignal2SRA=true is ignored.');
end

if bitUsePC && bitMatchSignal2SRA && any(ismember(var2Plot,{'RSV','ISV','ARSV'}))
    warning('When plotting RSV or ISV for PCs, all signals for each dimension are plotted (bitMatchSignal2SRA=false), so setting bitMatchSignal2SRA=true is ignored')
end

if bitSumDim && bitSepAx && any(ismember(var2Plot,'V'))
    error('Summing across dimensions onto separate axes for variance metric is not defined.');
elseif bitSumDim && bitSepAx && ~bitUsePC
    error('Summing across dimensions onto separate axes for RSV/ISV is possible for sRAs, but treats axes as separate SIGNALS not separate DIMENSIONS, which is the expected behavior. Thus we do not permit this.')
elseif bitSumDim && bitSepAx && bitUsePC
    error('Summing across dimensions onto separate axes for RSV/ISV is not possible for PCs.')
end

% determine if plotting each signal on a separate axis:
if bitSumDim
    % if summing across dimensions, do use separate axes
%     bitSepAx = false;
    
    % do not plot max variance:
    bitMaxVar = false;

elseif ~bitMatchSignal2SRA
    % if plotting RSV/ISV for all signals on each sRA, use separate
    % axes (unless summing across dimensions, which is already handled by
    % the above conditional)
    bitSepAx = true;

elseif bitUsePC && any(ismember(var2Plot,{'RSV','ISV','ARSV'}))
    % when plotting RSV or ISV for PCs, use separate axes
    bitSepAx = true;
elseif bitUsePC
    bitSepAx = false;
end


if bitRSVBoundDiff
    
    if bitMatchSignal2SRA
        warning('Cannot plot difference of RSV and RSV bound (bitRSVBoundDiff==true) when matching signal to sRA (bitMatchSignal2SRA==true)')
        bitRSVBoundDiff = false;
    end
    
    if bitSumDim
        warning('Cannot plot difference of RSV and RSV bound (bitRSVBoundDiff==true) when summing dimensions (bitSumDim==true)')
        bitRSVBoundDiff = false;
    end
    
    % force bitRSVBound = false
    bitRSVBound = false;
end    

% determine type of max variance to plot and set flags accordingly
bitMaxVar = false;
bitUseTDR4MaxVar = false;
bitUseSerialSRA4MaxVar = false;
bitUsePC4MaxVar = false;
switch maxVar2Plot
    case 'dPC'
        bitUsePC4MaxVar = true;
    case 'dRA'
        if ~isfield(Summary,tdr)
            % check that fields are available
            warning('Requested dRA for max var, but field %s not available. Using dPC instead',tdr)
            bitUsePC4MaxVar = true;
        elseif ~all(ismember(var2Plot,{'RSV','ISV','ARSV'}))
            % can only use dRA when plotting RSV, ISV
            warning('Requested dRA for max var, but can only use for RSV and ISV. Using dPC instead')
            bitUsePC4MaxVar = true;
        else
            bitUseTDR4MaxVar = true;
        end
    case 'serialSRA'
        if ~isfield(Summary,otdrSerial)
            % check that fields are available
            warning('Requested serial sRA for max var, but field %s not available. Using dPC instead',tdr)
            bitUsePC4MaxVar = true;
        else
            bitUseSerialSRA4MaxVar = true;
        end
    case 'none'
        % bitMaxVar = false;
        % do nothing
    otherwise
        error('maxVar2Plot %s not recognized',maxVar2Plot);
end

if any([bitUseTDR4MaxVar,bitUseSerialSRA4MaxVar,bitUsePC4MaxVar])
    bitMaxVar = true;
end

% if ~bitMaxVar && bitRSVBound
%     warning('Can only plot RSV bounds when plotting max variance for the on-target signal (see maxVar2Plot)');
% end


figH = figure;
whitebg(gcf,bgColor)
set(gcf,'Color',bgColor,'Name',['Variance Explained with chance level']);


% set varAnalysis field name
if eN == 1 || eN == 2
    % present offer vectors (usually fn = "oTDRSummary").
    % no need to change fn
    
    if eN == 1
        % present trial data
        fnVA = 'varAnalysis';
    elseif eN == 2
        % previous trial data
        fnVA = 'varAnalysis_PP';
    end
    
    % prefix for sRA labels
    sRAPrefix = '';
    
elseif eN == 3 || eN == 4 || eN == 5
    % previous offer vectors (usually fn = "oTDRSummary_PP").
    fn = [fn,'_PP'];
    
    if eN == 3
        % previous trial data (note different field name within previous
        % offer summary struct)
        fnVA = 'varAnalysis';
    elseif eN == 4
        % present trial data coded WRT previous offer conditions
        fnVA = 'varAnalysis_present';
    elseif eN == 5
        % present trial data coded WRT previous offer conditions
        fnVA = 'varAnalysis_presRespPrevCond';
    end
%     fnVA = 'varAnalysisPP';
%     fnSRA = 'sRA_PP';

    % add prefix for sRA labels. This refers to which variable the sRAs
    % were discovered WRT. 
    if ~isempty(condLabel) 
        condLabel = cellfun(@(c) ['Prev ',c],condLabel,'UniformOutput',0);
    end
else
    error('epoch number not recognized')
end

% depending on temporal epoch, set x-axis label and determine intervals to
% plot
if eN == 1 || eN == 4 || eN == 5
    
    % set name of alignment to use for event intervals
    alignName = 'Offer';
    
    if eN == 5
        % TDR in present trial but based on previous trial conditions
        tdr = [tdr,'_presRespPrevCond'];
    end
    
elseif eN == 2 || eN == 3
   
    % set name of alignment to use for event intervals
    alignName = 'Fixation';
    
%     if eN == 3
        % TDR in previous trial and based on previous trial conditions.
        % This works for comparing to both present-trial and previous-trial
        % sRAs, since the regressors are the same in the previous trial.
        tdr = [tdr,'_PP'];
%     else
%         % currently do not have TDR for previous trial data on present
%         % trial conditions
%         tdr = [];
%         
%         if ~strcmp(maxVar2Plot,'none')
%             error('Cannot use a TDR-based var for maxVar2Plot when plotting previous trial data projected on present trial sRAs')
%         end
%     end    
end

% set xlabel (which depends on alignment)
xLabel = sprintf('Time from %s (s)',lower(alignName));;

% set type of regression axis (sRA vs PC) to use
if bitUsePC
    if bitUsePC_CC
        fnRA = 'PCs_CC';
    else
        fnRA = 'PCs';
    end
else
    fnRA = 'sRA';
end    

% get colors -- 
if strcmp(bgColor,'w')
    % new colors
    sRAColor = colorCategorical(length(condLabel),studyName);
    PCColor = [0 0 0];
else
    sRAColor = colorCategorical(length(condLabel));
    PCColor = [1 1 1];
end
% repeat colors if signal is repeated
% sRAColor = sRAColor(Summary.(otdr).(fnRA).param4RA,:);
    
% SPECIAL CASE, when plotting R2 term (proportional to correlation
% between projection and predictor) or normdRA term
if any(ismember(var2Plot,{'R2','normdRA'}))
    x = Summary.(tdr).Times.regressTimes;
    
    % for certain other parameters
    varType = 2;
    bitUsePC = false;
else
    switch eN
        case 1 % present trial data from present offer summary struct
            x = Summary.(fn).Times.all_times';
        case 2 % previous trial data from present offer summary struct
            x = Summary.(fn).Times.all_times_PP';
        case 3 % previous trial data from previous offer summary struct
            x = Summary.(fn).Times.all_times';
        case {4,5} % present trial data from previous offer summary struct
            x = Summary.(fn).Times.all_times_present';
        otherwise
            error('epoch %d not recognized',eN)
    end
end

% timebins
dt = x(2) - x(1);

% number of RAs -- set by number of sRAs
% nRA = length(condLabel);
nRA = size(Summary.(otdr).sRA.RA,2);
    

% compute chance level for Variance (or for RSV or ISV, not both, when
% plotting PCs), and only when plotting absolute var or var explained, and
% not when plotting R2
chY = [];
if ~any(ismember(var2Plot,{'R2','normdRA'})) && ...
        (~all(ismember({'RSV','ISV','ARSV'},var2Plot)) || (bitSumDim && varType==4)) ... % makes exception for bitSumDim==1 & varType==4 because we need to extract chance distribution so as to compute p(RSV/ISV)  
        && ~(((bitUsePC && ~bitSumDim) || (~bitMatchSignal2SRA && ~bitSumDim)) && any(ismember(var2Plot,{'RSV','ISV','ARSV'}))) ... % makes exception for bitSumDim==1 because we need to extract chance distribution so as to compute p(RSV/ISV) 
        && (ismember(varType,[1 2]) || (varType == 4 && bitSumDim)) % This gets gamma fits even when varType == 4 but summing variance across dimensions. That way we can compute p-values below
%         (~all(ismember({'RSV','ISV','ARSV'},var2Plot)) || (bitUsePC && length(var2Plot) == 1)) SINCE RSV FOR EACH PC NOW ARE PLOTTED TOGETHER, WE CAN'T PLOT CHANCE FOR EACH RSV...

    if ~isfield(Summary.(otdr).(fnRA).varAnalysis,'varRand_t')
        warning('Expected field Summary.%s.%s.varAnalysis.varRand_t was not present and thus no chance level is plotted.',otdr,fnRA);
    else
        %     if bitSigma
        %         % 2 sigma
        %         chY = (mean(Summary.(otdr).(fnRA).varAnalysis.varRand_t,2) ...
        %             + chanceLevel*sqrt(var(Summary.(otdr).(fnRA).varAnalysis.varRand_t,[],2)));
        %     else
        %         % 95th percentile:
        %         chY = prctile(Summary.(otdr).(fnRA).varAnalysis.varRand_t,100*(1-alpha),2);
        %     end
        
        %%%%%%%%%%%%%%%%%%
        % extract random vectors
        varRand = Summary.(otdr).(fnRA).varAnalysis.varRand_t;
        RSVRand = Summary.(otdr).(fnRA).varAnalysis.RSVRand_t;
        ISVRand = Summary.(otdr).(fnRA).varAnalysis.ISVRand_t;
        
        % normalize random vectors if plotting percent variance explained
        if varType == 2
            varRand = 100 * varRand ./ ...
                repmat(Summary.(otdr).(fnRA).varAnalysis.(totVarFN),1,...
                size(varRand,2));
            RSVRand = 100 * RSVRand ./ ...
                repmat(Summary.(otdr).(fnRA).varAnalysis.(totVarFN),1,...
                size(RSVRand,2),size(RSVRand,3));
            ISVRand = 100 * ISVRand ./ ...
                repmat(Summary.(otdr).(fnRA).varAnalysis.(totVarFN),1,...
                size(ISVRand,2),size(ISVRand,3));
        end
        
        % fit gamma distributions to chance distributions
        nTime = size(varRand,1);
        nSignal = size(Summary.(otdr).(fnRA).varAnalysis.RSVRand_t,3);
        gamma_V = NaN(nTime,2);
        gamma_RSV = NaN(nTime,nSignal,2);
        gamma_ISV = NaN(nTime,nSignal,2);
        ch_V = NaN(nTime,1);
        ch_RSV = NaN(nTime,nSignal);
        ch_ISV = NaN(nTime,nSignal);
        % loop through time points
        for i = 1:nTime
            gamma_V(i,:) = gamfit(varRand(i,:));
            ch_V(i) = gaminv(1-alpha,gamma_V(i,1),gamma_V(i,2));
            % loop through signals
            for j = 1:nSignal
                gamma_RSV(i,j,:) = gamfit(squeeze(RSVRand(i,:,j)));
                ch_RSV(i,j) = gaminv(1-alpha,gamma_RSV(i,j,1),gamma_RSV(i,j,2));
                gamma_ISV(i,j,:) = gamfit(squeeze(ISVRand(i,:,j)));
                ch_ISV(i,j) = gaminv(1-alpha,gamma_ISV(i,j,1),gamma_ISV(i,j,2));
            end
        end
        
        if all(strcmp(var2Plot,'RSV'))% && bitUsePC
            chY = ch_RSV;
        elseif all(strcmp(var2Plot,'ISV')) %&& bitUsePC
            chY = ch_ISV;            
        elseif all(strcmp(var2Plot,'V'))
            chY = ch_V;
        end
    end
end

% instantiate
legText = {};
h = [];
axH = [];
axH2 = [] ; % used for plotyy
yLimAll = [];


%%%%%%
% Loop through different types of variance requested
for v = 1:length(var2Plot)

    % Label for type of variance. Since there can be multiple types, we leave
    % the label generic fornow:
    varLabel = 'Variance';
    % switch var2Plot{v}
    %     case 'V'
    %         varLabel = 'Variance';
    %     case 'RSV'
    %         varLabel = 'Relevant Variance';
    %     case 'ISV'
    %         varLabel = 'Irrelevant Variance';
    % end
    
    
    % check that only plot relevant variance in z-score form since the
    % chance distribution is only defined for absolute variance
    % if ~strcmp('V',var2Plot{v}) && varType~=3
    %     error('When plotting %s, you may only report variance in terms of deviation from chance (z-score)',var2Plot{v});
    % end
    
    % define variance field to plot
    % and determine how to set chance level
    
    varFN = var2Plot{v};
    varMaxFN = var2Plot{v};
    
            
    % Set fields to extract based on type of variance requested
    if varType == 1 || (varType == 4 && bitSumDim == 1)
        % absolute variance, OR SPECIAL CASE when plotting probability of
        % variance summed across dimensions, we have to compute that
        % probability de novo, since it is not returned by TDR/oTDR
        % functions
        yLabel = ['Absolute ',varLabel];
    elseif varType == 2 
        % percent variance
        varFN = [varFN,'E'];
        varMaxFN = [varMaxFN,'E'];
        
        yLabel = [varLabel,' explained (%)'];
        
        % noramlize chance level if plotting % variance explained and scale
        % by 100
        %             chY = 100 * chY ./ Summary.(otdr).(fnRA).varAnalysis.(totVarFN);
        
        %         case 3 % z-scored
        %             varFN = ['z',varFN];
        %             varMaxFN = ['z',varMaxFN];
        %
        %             yLabel = [varLabel,': Deviation from chance (z-score)'];
        %
        %             % replace chance level with Y = 2
        % %             chY = ones(size(x)) * chanceLevel;
        %             %         chY(:) = chanceLevel;
        
    elseif varType == 4
        % p-value
        %         if ~ismember(var2Plot{v},{'RSV','ISV','ARSV'})
        %             error('Can only plot p-value for RSV and ISV')
        %         end
        
        varFN = ['p',varFN];
        varMaxFN = ['p',varMaxFN];
        
    end
    
    % separately set chance level and yLabel when varType == 4. We do not
    % do this above because of special case when varType==4 && bitSumDim,
    % in which case we want to set field names as though computing absolute
    % variance
    if varType == 4
        % replace chance level with Y = 1 - alpha.
        chY = ones(size(x)) * 0.05;
        %         chY = ones(size(x)) * 10.^-log10(0.05);
        %         chY(:) = 10.^-log10(0.05);
        
        %         yLabel = [varLabel,': 10^{-log(p-value)}'];
        yLabel = sprintf('P(%s)',varLabel);
    end

    
    % varFN = [varFN,'_sRA'];
    varFN = [varFN,'_RA'];
    
    % add "partial" suffix if requested
    if bitRSVPartial && ~bitMatchSignal2SRA ...
                && ismember(var2Plot(v),{'RSV','ISV','ARSV'}) && ~bitSumDim  
        % save the non-partial name
        varFNOrig = varFN;
        % and add suffix
        varFN = [varFN,'_partial'];
        
        % flag
        bitRSVPartialUsed = true;
        
    else
        varFNOrig = [];
        % flag
        bitRSVPartialUsed = false;
    end
    
    % SPECIAL CASE, when plotting R2 term (proportional to correlation
    % between projection and predictor)
    if strcmp('R2',var2Plot{v})
        fn = 'TDRSummary';
        fnRA = 'runTDRSummary';
        varFN = 'R2_tk';
        varMaxFN = [];
        fnVA = [];
        yLabel = 'R2: contribution of signal';
    % SPECIAL CASE, when plotting R2 term (proportional to correlation
    % between projection and predictor)
    elseif strcmp('normdRA',var2Plot{v})
        fn = 'TDRSummary';
        fnRA = [];
        varFN = 'normdRAs';
        varMaxFN = [];
        fnVA = [];
        yLabel = 'norm dRA: signal magnitude';
    end
    
    
    hold on;
    

    % set how many times to loop through signals. When we are summing
    % across dimensions and plotting variance of sRA/PC, there will only be 1
    % line and we should loop = 1. When summing but plotting RSV/ISV,
    % there will be 3 lines, and requires loop = 3, one for each signal. 
    if bitSumDim && ...
            (ismember(var2Plot(v),{'V'})) % || ...
%             (ismember(var2Plot(v),{'RSV','ISV','ARSV'}) && ~bitSepAx))
%             && ~(bitUsePC && ismember(var2Plot(v),{'RSV','ISV','ARSV'}))
        nLoopSignal = 1;
    else
        nLoopSignal = nRA;
    end
    
    % initialize logical on which signals (columns) are on or off-target
    % (rows)
    bitTarget = false(nLoopSignal,length(condLabel));
    
    % loop through conditions
    clear foo y
    
    for i = 1:nLoopSignal
        
        % if plotting on separate axes
        if bitSepAx
            % plot each signal on a separate subplot
            if v==1 % first variance type to plot
                axH(i) = subplot(1,nRA,i);
                if bitUsePC
                    if bitUsePC_CC
                        title(sprintf('CC PC_%d',i),'FontSize',16)
                    else
                        title(sprintf('PC_%d',i),'FontSize',16)
                    end
                else
                    title(condLabel{Summary.(otdr).(fnRA).param4RA(i)},'FontSize',16)
                end
                hold on;
                
                % when we makes new axes, reset the line count for that
                % axes
                nLine(i) = 0; % number of lines PER AXIS

            else
                axes(axH(i));
            end
            axN = i;
        else
            % make axes in first signal on first type of var
            if i == 1 && v == 1
                axN = 1;

                % when we makes new axes, reset the line count. in the case
                % of a single axis, this only happens once
                nLine = 0; % number of lines PER AXIS
                
            end
        end
        
        y=[];
        
        %%% DEFUNCT -- used to plot PCs separately from sRAs. Now code is
        %%% consolidated below
%         % set set of axes (dimension in neural space) to plot
%         if bitUsePC && bitSepAx
%             
%             % here we're going to treat PCs a little differently, instead
%             % of looping through each PC for a given signal (thus grouping
%             % dimensions per signal), we will loop through each signal for
%             % a given PC (thus grouping signals per dimension), which is
%             % akin to how we plot RSV of each signal per sRA. This means we
%             % will hijack the signal loop (i) as the PC loop, and here loop
%             % through signals (s). This relies on there being the same
%             % number of PCs as signals in the var analysis
%             if size(Summary.(fn).(fnRA).(fnVA).RSV_RA,2) ~= size(Summary.(fn).(fnRA).(fnVA).RSV_RA,3)
%                 error('Number of PCs used for var analysis must match number of signals')
%             end
%             
%             if ismember(var2Plot{v},{'RSV','ISV','ARSV'})
%                 
%                 % loop through each signal for PC=i
%                 for s = 1:3
%                     % when current type of variance is RSV/ISV
%                     if ndims(Summary.(fn).(fnRA).(fnVA).(varFN)) == 3
%                         y(:,s) = Summary.(fn).(fnRA).(fnVA).(varFN)(:,i,s);
%                     else
%                         error('number of dimensions (%d) nor as expected',ndims(Summary.(fn).(fnRA).(fnVA).(varFN)))
%                     end
%                 end
%             else
%                 % this really shouldn't happen, because if summing across
%                 % dimensions and plotting 'V' (i.e., not separately by
%                 % signal) then an error should have occurred above
%                 error('Summing across dimensions for variance "V" is not defined.')
% 
% %                 % plot single metric for PC = i.
% %                 if ndims(Summary.(fn).(fnRA).(fnVA).(varFN)) == 2
% %                     y(:,p) = Summary.(fn).(fnRA).(fnVA).(varFN)(:,i);
% %                     
% %                 else
% %                     error('number of dimensions (%d) nor as expected',ndims(Summary.(fn).(fnRA).(fnVA).(varFN)))
% %                 end
%                  
%             end
%             
%             
%             % when summing across dimensions, sum across PCs
%             if bitSumDim
%                 y = sum(y,2);
%             end
%             
%         % plotting non-PCs
%         else

        % if summing across dimensions, extract all dimensions in one
        % pass (this assuming we are not looping through each
        % dimensions)
        if bitSumDim
            ind = 1:nRA;
        else
            ind = i;
        end

        if ((~bitMatchSignal2SRA && ~strcmp(var2Plot{v},'V')) ...
                || (bitUsePC && ismember(var2Plot(v),{'RSV','ISV','ARSV'}))) && ~bitSumDim
            % if plotting RSV/ISV for each signal regardless of sRA,
            % extract all signals and plot them simultaneously for each
            % sRA (1 sRA per loop). This is also the default behavior
            % for plotting PCs: we always plot the RSV/ISV for all
            % signals for each PC. The exception to this is when we are
            % summing across dimensions (PCs or sRAs) and plotting
            % RSV/ISV. In this case, we always plot RSV/ISV for all
            % signals, but of course there is only a single "dimension"
            % (really the sum across dimensions). In this case, instead
            % of using the loop to access different dimensions, we use
            % the loop to access different signals and plot just one
            % per loop.
            % Note, in this case ind2 refers to signals (condLabel), not
            % sRAs. Thus we set it according to the number of signals
            ind2 = 1:length(condLabel);
        else
            % if plotting RSV/ISV for the signal specific to the sRA,
            % or when plotting RSV/ISV summed across dimensions
            % (plotting 1 signal per loop).
            % Note, in this case ind2 refers to signals (condLabel), not
            % sRAs. Thus we set it according to the signal number
            if bitUsePC
                % in case of PCs, the assignment between PC number and
                % signal is arbitrary
                ind2 = min(i,length(condLabel));
            else
                ind2 = Summary.(fn).(fnRA).param4RA(i);
            end
        end

        % SPECIAL CASE, when plotting R2 term (proportional to correlation
        % between projection and predictor)
        if strcmp('R2',var2Plot{v})
            y = Summary.(fn).(fnRA).(varFN)(:,ind);
            % normalize by max
            y = y / max(y);
            %                 % normalize by z-score
            %                 y = zscore(y);
            % SPECIAL CASE, when plotting normdRA term
        elseif strcmp('normdRA',var2Plot{v})
            y = Summary.(fn).(varFN)(:,ind);
            % normalize by max
            y = y / max(y);
            %                 % normalize by z-score
            %                 y = zscore(y);
        else
            nDim = ndims(Summary.(fn).(fnRA).(fnVA).(varFN));
            if nDim == 2
                y = Summary.(fn).(fnRA).(fnVA).(varFN)(:,ind);
            elseif nDim >= 3
                
                % decide which are the on-target signals (used for full
                % correlation or for reference) and which are the
                % off-target signals (used for partial correlation or as
                % difference with bounds)
                if bitUsePC
                    % special case for PCs...
                    % in the case of PCs, we will define the on-target
                    % predictor as that with the highest RSV for the given PC
                    
                    if bitRSVPartialUsed
                        % get field name prior to adding partial suffix
                        soo = varFNOrig;
                    else
                        soo = varFN;
                    end
                    % collect full RSV/ISV info
                    if varType == 4
                        % when plotting p-value, still select axes based on
                        % max RSV/ISV
                        goo = strrep(soo,'p','');
                        temp = Summary.(fn).(fnRA).(fnVA).(goo)(:,ind,ind2);
                    else
                        temp = Summary.(fn).(fnRA).(fnVA).(soo)(:,ind,ind2);
                    end
                    clear goo soo
                    
                    % find axes with max RSV/ISV value
                    [~,foo] = max(mean(temp));

                    bitTarget(i,foo) = true;
                    clear foo temp
                elseif length(ind2) == size(bitTarget,2) 
                    % use full correlation when signal matches axes, and
                    % use partial correlation otherwise
                    bitTarget(i,:) = ind2 == Summary.(fn).(fnRA).param4RA(ind);
                elseif length(ind2) == 1
                    bitTarget(i,ind2) = true;
                end
                
                if sum(bitTarget(i,:)) ~= 1
                    error('only expect 1 signal to match sRA')
                end
                
                if nDim == 3
                    y = Summary.(fn).(fnRA).(fnVA).(varFN)(:,ind,ind2);
                elseif nDim == 4
                    % for RSV/ISV partial

                    % intialize
                    y = squeeze(NaN(size(Summary.(fn).(fnRA).(fnVA).(varFNOrig)(:,ind,ind2))));


                    % target signal
                    y(:,bitTarget(i,:)) = Summary.(fn).(fnRA).(fnVA).(varFNOrig)(:,ind,bitTarget(i,:));
                    % off-target signals --
                    pos = Summary.(fn).(fnRA).(fnVA).signalN4Partial(bitTarget(i,:)',:);
                    pos = pos == ind2(~bitTarget(i,:));
                    y(:,~bitTarget(i,:)) = Summary.(fn).(fnRA).(fnVA).(varFN)(:,ind,bitTarget(i,:),pos);
                end
            else
                error('number of dimensions (%d) nor as expected',ndims(Summary.(fn).(fnRA).(fnVA).(varFN)))
            end
            
            %%%
            % if plotting ARSV (and not p(ARSV) with other signals on same
            % axes, convert ARSV to negative values
            if strcmp('ARSV',var2Plot{v}) && varType~=4 && length(var2Plot) > 1
%                 y = -y;
            end
        end

        % if summing across dimensions, do that here
        if bitSumDim
            y = sum(y,2);
        end

        % if plotting more than one signal per sRA/PC, we may need to
        % squeeze y
        if ~bitMatchSignal2SRA || (ismember(var2Plot(v),{'RSV','ISV','ARSV'}))
            y = squeeze(y);
        end
%         end

        % PARTIAL CORRELATIONS -- when plotting RSV/ISV for off-target
        % terms, covering RSV/ISV to RSVp/ISVp, where RSVp = var(P) *
        % r_bP|a, where P is the data projected onto the axes, a is the
        % on-target predictor, b is an offtarget predictor, and r_bP|a is
        % the partial correlation corr(res_b,res_p), where res_b and res_p
        % are the residuals of b and p after controlling for corr(a,b) and
        % corr(a,p), respectively.
%         if bitRSVPartial && ~bitMatchSignal2SRA ...
%                 && ismember(var2Plot(v),{'RSV','ISV','ARSV'}) && ~bitSumDim ...
% 
%             % flag
%             bitRSVPartialUsed = true;
%             
%             % find the RSV that will serve as the on-target predictor
%             if bitUsePC
%                 % in the case of PCs, we will define the on-target
%                 % predictor as that with the highest RSV for the given PC
%                 [~,foo] = max(mean(y));
%                 bitTarget = false(1,size(y,2));
%                 bitTarget(foo) = true;
%                 clear foo
%             else
%                 % in the case of sRAs, the on-target predictor is
%                 % well-defined.
%                 bitTarget = false(1,size(y,2));
%                 bitTarget(ind) = true;
%             end
%             
%             %%% Here we define on-target predictor = a, off-target
%             %%% predictor = b, variance explained by axis = p, r_xy =
%             %%% corr(x,y), r_xy|z = partialcorr(x,y,z). We are intereted
%             %%% specifically in r_bp|a. It is convenient to use an
%             %%% alternative means to compute partialcorr(), as: 
%             %%% rbp|a = (r_bp - r_ap*r_ab) / sqrt((1-r_ap^2)*(1 - r_ab^2))
%             %%% see http://stats.stackexchange.com/questions/76815/multiple-regression-or-partial-correlation
% 
%             % extract variance explained by axis
%             if varType == 1
%                 VE = Summary.(fn).(fnRA).(fnVA).V_RA(:,ind);
%             elseif varType == 2
%                 VE = Summary.(fn).(fnRA).(fnVA).VE_RA(:,ind);
%             end
%             
%             % find r_ap as RSV_a = var(p) * r_ap^2, thus
%             % r_ap = sqrt(RSV_a/var(p))
%             r_ap = sqrt(y(:,bitTarget) ./ VE);
%             
%             % loop through off-target terms
%             for j = find(~bitTarget)
%                 % compute correlation between on and off-target predictor.
%                 % Must transform predictors in the same way as they were
%                 % for computing RSV (see condSeparation.m), where they are
%                 % mean-centered and normalized to unit vectors. (Actually,
%                 % this may not make any difference).
%                 a = Summary.(otdr).codedParams(:,bitTarget);
%                 b = Summary.(otdr).codedParams(:,j);
%                 a = a - mean(a);
%                 b = b - mean(b);
%                 a = a / norm(a);
%                 b = b / norm(b);
%                 r_ab = corr(a,b);
%                 
%                 % compute r_bp as 
%                 % r_bp = sqrt(RSV_b/var(p))
%                 r_bp = sqrt(y(:,j)./VE);
%                 
%                 % compute partial correlation r_bp|a as:
%                 % r_bp|a = (r_bp - r_ap*r_ab) / sqrt((1-r_ap^2)*(1 - r_ab^2))
%                 r_bp_a = (r_bp - r_ap.*r_ab) ./ sqrt((1-r_ap.^2)*(1 - r_ab.^2));
%                 
%                 % recompute RSV_b = var(p) * r_bp|a.^2
%                 y(:,j) = VE .* r_bp_a.^2;
%                 
%                 % compare to newly computed value
%                 if any(abs(y(:,j) - ...
%                         Summary.(fn).(fnRA).(fnVA).([varFN,'_partial'])(:,ind,ind,(Summary.(fn).(fnRA).(fnVA).signalN4Partial(ind,:) == j))) > eps*100)
%                     error('partial correlation computation error')
%                 end
%             
%             end
%             clear VE r_ap r_ab r_bp r_bp_a a b
%         else
%             % flag
%             bitRSVPartialUsed = false;
%             
%         end

        % SPECIAL CASE when plotting probability of variance summed across
        % dimensions, we have to compute that probability de novo, since it
        % is not returned by TDR/oTDR functions
        if varType == 4 && bitSumDim == 1
            
            % at this point, y is the summed absolute variance/RSV/ISV
            % across dimensions
            
            % save y and re-instantiate
            ySumVar = y;
            y = NaN(size(y));

            % determine null distribution based on var2Plot
            if strcmp(var2Plot{v}, 'V')
                y = gamcdf(ySumVar,gamma_V(:,1),gamma_V(:,2),'upper');
            else
                if strcmp(var2Plot{v},'RSV')
                    gammaTemp = gamma_RSV;
                elseif strcmp(var2Plot{v},'ISV')
                    gammaTemp = gamma_ISV;
                else
                    error('Variance "%s" not recognized',var2Plot{v})
                end
                
                % if plotting RSV/ISV, compute p-value relative to null
                % distribution of current signal (given by ind2, which
                % should be the same as i)
%                 for j = 1:size(y,2)
                    y = gamcdf(ySumVar,gammaTemp(:,ind2,1),gammaTemp(:,ind2,2),'upper');
%                 end
                
            end                    
            
            
        end 
        
        % if plotting difference of RSV and RSV bound, find bound and use
        % plotyy
        if bitRSVBoundDiff % && size(y,2) > 1
%             % loop through each offtarget predictor
%             for j = find(~bitTarget(i,:))
%                 % compute correlation between on and off-target
%                 % predictor
%                 % compute correlation between on and off-target predictor.
%                 % Must transform predictors in the same way as they were
%                 % for computing RSV (see condSeparation.m), where they are
%                 % mean-centered and normalized to unit vectors
%                 a = Summary.(otdr).codedParams(:,bitTarget(i,:));
%                 b = Summary.(otdr).codedParams(:,j);
%                 a = a - mean(a);
%                 b = b - mean(b);
%                 a = a / norm(a);
%                 b = b / norm(b);
%                 
%                 rij = corr(a,b);
%                 
%                 % compute upper bound as VE * rij^2 and subtract from
%                 % off-target signal
%                 y(:,j) = y(:,j) - Summary.(fn).(fnRA).(fnVA).VE_RA(:,i)*rij^2;
%             end
            

            % extract variance explained by axis depending on absolute
            % variance or normalized variance
            if varType == 1
                VE = Summary.(fn).(fnRA).(fnVA).V_RA(:,ind);
            elseif varType == 2
                VE = Summary.(fn).(fnRA).(fnVA).VE_RA(:,ind);
            end
            
            % store the on-target RSV for later plotting
            yOnTarg = y(:,bitTarget(i,:));
            
            % replace off-target RSV with difference between off-target RSV
            % for signal B and RSV bound based on signals A and B, which is
            % based on variance VE of projection P of data onto sRA.
            % RSV = r_BP^2 * VE
            % RSV_bound = r_AB^2 * VE
            % RSVDiff = RSV - RSV_bound
            %         = VE * (r_BP^2 - r_AB^2)
            if varType == 4
                % for probabilities, we can only plot p(RSVDiff) for
                % off-target signals
                y(:,~bitTarget(i,:)) = Summary.(fn).(fnRA).(fnVA).pDiffCorr_proj_signal(:,ind,bitTarget(i,:),:);
                
                % replace on-target RSV with NaNs so that there's a place
                % holder for that signal and a line handle
                y(:,bitTarget(i,:)) = NaN;
                
                % clean p==0 values (if any)
                y = cleanPZero(y,varType);
                % plot without ontarget RSV
                foo = plot(x,y);
            else
                % for RSV and RSVE
                y(:,~bitTarget(i,:)) = VE .* Summary.(fn).(fnRA).(fnVA).diffCorr_proj_signal(:,ind,bitTarget(i,:),:);

                % replace on-target RSV with difference between on-target RSV
                % and upperbound (VE).
                y(:,bitTarget(i,:)) = y(:,bitTarget(i,:)) - VE;
                
                % plotyy with difference of RSV-Bound on left and on-target RSV
                % on right for reference
                % clean p==0 values (if any)
                y = cleanPZero(y,varType);
                yOnTarg = cleanPZero(yOnTarg,varType);
                [axTemp,goo1,goo2] = plotyy(x,y,x,yOnTarg);
                % order line handles
                clear foo
                %             foo(bitTarget) = goo1;
                %             foo(~bitTarget) = goo2;
                foo = [goo1;goo2];
                clear goo1 goo2
                % store second axes
                axH2(end+1) = axTemp(2);
                % set second axes to next plot = add
                set(axTemp(2),'NextPlot','add');
                clear axTemp
                
            end
            clear VE

%             
%             foo = plot(x,y);
            
        else
            % clean p==0 values (if any)
            y = cleanPZero(y,varType);
            foo = plot(x,y);
        end
        
        % check number of line handles
        if bitRSVBoundDiff 
            if varType == 4 && length(foo) ~= length(ind2)
                % special case for probability of RSV Bound Diff. 
                error('When plotting probability of RSVBoundDiff, number of line handles (%d) should match number of indecies "ind2" (%d)',length(foo),length(ind2))
            elseif varType ~= 4 && length(foo) ~= length(ind2) + 1
                error('When including RSVBoundDiff, number of line handles (%d) should match number of indecies "ind2" (%d) + 1',length(foo),length(ind2))
            end
        else
            if length(foo) ~= length(ind2)
                error('Number of line handles (%d) should match number of indecies "ind2" (%d)',length(foo),length(ind2))
            end
        end
        
        % collect handles (each row is a different axis)
        h(axN,nLine(axN)+1:nLine(axN)+numel(foo)) = foo;
        
        % this may have created an axes. if so grab the handle
        if isempty(axH)
            axH = gca;
        end
        
        % set line color 
        if (bitUsePC || bitSumDim) && ~ismember(var2Plot(v),{'RSV','ISV','ARSV'}) % indicates plotting variance, not signal-specific
            % single axis, one color for all PCs
            lc = PCColor;
%         elseif bitUsePC || ~bitMatchSignal2SRA
%             % separate signals for each axes
%             lc = sRAColor(ind2,:); % expect set of colors
%         else
%             % separate axes for each signal
%             lc = sRAColor(i,:);
        else
            lc = sRAColor(ind2,:);
        end
        
        % set linestyle based on type of variance
        if strcmp(var2Plot{v},'ISV')
            ls = '--';
        elseif (strcmp(var2Plot{v},'V') | strcmp(var2Plot{v},'ARSV')) && ismember('RSV',var2Plot)
            % plotting V or ARSV when RSV is also present
            ls = '-.';
        else
            % plotting V when RSV is NOT present, or plotting RSV
            ls = '-';
        end
            
        % may have multiple colors, one for each line
        for ii = 1:size(lc,1)
            % if only one color, apply it to all lines. If multiple colors,
            % assume color(pos) --> line(pos)
            if size(lc,1) > 1
                iii = ii;
            else
                iii = 1:length(foo);
            end
            
            % if plotting Diff(RSV - RSV Bound), change line style to
            % dotted
            if bitRSVBoundDiff && size(lc,1) > 1 % && ~bitTarget(i,ii)
                lsTemp = '-.';
            else
                lsTemp = ls;
            end
            set(foo(iii),'Color',lc(ii,:),'LineWidth',2,'LineStyle',lsTemp);
            
        end
        % if plotting RSV difference (and not the probability of the
        % difference), then specify appearance of extra line for on-target
        % RSV
        if bitRSVBoundDiff && varType ~= 4
            set(foo(size(lc,1)+1),'Color',lc(bitTarget(i,:),:),'LineWidth',2,'LineStyle',ls);
        end
        
        % if plotting PCs on same axes, distinguish by line width
        if bitUsePC && ~bitSepAx && ~bitSumDim
            set(foo,'LineWidth',3-(i-1));
        end
        clear ii iii
        
        if bitRSVBoundDiff && varType ~= 4
            % when including RSV Diff (and not p(RSV Diff)) and on-target RSV,
            % force loop to go through a second time for the on-target RSV
            foo = [ind2, find(bitTarget(i,:))];
            
%             if varType == 4
%                 % when plotting probability of RSV Diff, remove the
%                 % on-target signal from the list of signals used to
%                 % generate the legend
%                 foo = setdiff(ind2,find(bitTarget(i,:)));
%             else
%                 % when including RSV Diff (and not p(RSV Diff)) and on-target RSV,
%                 % force loop to go through a second time for the on-target RSV
%                 foo = [ind2, find(bitTarget(i,:))];
%             end

        else
            foo = ind2;
        end
        
        % legend 
        lineNTemp = 0;
        for ii = foo
            % this should be gauranteed to have at most one
            % line at a time
            lineNTemp = lineNTemp + 1;
            
            % When plotting RSV/ISV, we have to specify for which
            % signals
            if ismember(var2Plot(v),{'RSV','ISV','ARSV'})
                % SPECIAL CASE
                % Epoch 4 projects present trial data on previous trial
                % vectors but measures variance WRT present trial
                % conditions. Therefore, we remove "Prev" from the
                % description of RSV, since the RSV is WRT present trial
                % conditions
                if eN == 4
                    soo = strrep(condLabel{ii},'Prev ','');
                else
                    soo = condLabel{ii};
                end
                signal4Var = sprintf(' %s',soo);
                clear soo
            else
                signal4Var = '';
            end
            
            if bitSumDim
                % sum across dimensions
                legText{axN,nLine(axN)+lineNTemp} = sprintf('%s%s sum(%s 1-%d), mean = %0.1f',var2Plot{v},signal4Var,strrep(fnRA,'_',' '),nRA,mean(y(x>=0,lineNTemp)));
                
            else
                if bitRSVPartialUsed && ~bitTarget(i,ii)
                    % rename RSV/ISV field if using partial correlation for
                    % off-target signal
                    var2PlotTemp = sprintf('%s_{partial}',var2Plot{v});
                elseif lineNTemp > length(ind2) && (bitRSVBoundDiff && varType ~= 4) ...
                        && bitTarget(i,ii)
                    % label extra line as on-target signal 
                    var2PlotTemp = var2Plot{v};
                elseif bitRSVBoundDiff && ~bitTarget(i,ii)
                    % label as difference between off-target signal and
                    % bound
                    var2PlotTemp = sprintf('Diff(%s - RSV Bound)',var2Plot{v});
                    if varType == 4
                        % add prefix for probability
                        var2PlotTemp = ['p',var2PlotTemp];
                    end
                elseif bitRSVBoundDiff && bitTarget(i,ii)
                    if varType == 4 
                        % special flag to skip line handle and legend
                        % entry
                        var2PlotTemp = 'SKIP';
                    else
                        % label as difference between on-target signal and
                        % variance bound (only when not plotting probability of
                        % difference)
                        var2PlotTemp = sprintf('Diff(%s - VE Bound)',var2Plot{v});
                    end
                else
                    var2PlotTemp = var2Plot{v};
                end
                
                % if plotting sRA/PCs on separate axes, then we do not need to
                % include the sRA name/PC num in the legend
                if bitSepAx
                    legText{axN,nLine(axN)+lineNTemp} =  sprintf('%s%s',var2PlotTemp,signal4Var);
                else
                    if bitUsePC
                        if bitUsePC_CC
                            goo = 'CC ';
                        else 
                            goo = '';
                        end
                        legText{axN,nLine(axN)+lineNTemp} =  sprintf('%s%s on %sPC_%d',var2PlotTemp,signal4Var,goo,i);
                        clear goo
                    else
                        legText{axN,nLine(axN)+lineNTemp} =  sprintf('%s%s on sRA %s',var2PlotTemp,signal4Var,condLabel{Summary.(otdr).(fnRA).param4RA(i)});
                    end
                end
            end
        end
            
        
        %%% THE FOLLOWING IS ALL DEFUNCT NOW THAT WE'RE PLOTTING PCS ON
        %%% SEPARATE AXES
%         % if plotting PCs...
%         if bitUsePC
%             if numel(foo) == 1
%                 if bitSumDim
%                     % sum across dimensions
%                     
%                     if ismember(var2Plot(v),{'RSV','ISV','ARSV'})
%                         % when plotting separate lines for each signal
%                         legText{axN,nLine(axN)+1} = sprintf('%s %s sum(PC 1,2,3), mean = %0.1f',var2Plot{v},condLabel{i},mean(y(x>=0)));
%                     else
%                         % when plotting a single line (e.g., for variance) 
%                         legText{axN,nLine(axN)+1} = sprintf('%s sum(PC 1,2,3), mean = %0.1f',var2Plot{v},mean(y(x>=0)));
%                     end
%                 else
%                     % and handle is singular, then we are plotting PCs per each
%                     % loop through condLabel (and not actually plotting the
%                     % signals). 
%                     set(foo,'LineWidth',3-(i-1))
%                     legText{axN,nLine(axN)+1} = sprintf('%s PC_%d',var2Plot{v},i);
%                 end
%             else                
%                 % Otherwise, all three PCs are within handle.
%                 for p = 1:numel(foo)
%                     set(foo(p),'LineWidth',3-(p-1))
%                     legText{axN,nLine(axN)+p} = sprintf('%s PC_%d',var2Plot{v},p);
%                 end
%             end
%         else
%             if bitSumDim && ismember(var2Plot(v),{'RSV','ISV','ARSV'})
%                 % summing across dimensions, plotted by separate signals, report mean
%                 legText{axN,nLine(axN)+1:nLine(axN)+numel(foo)} =  sprintf('%s %s sum(sRA 1,2,3), mean = %0.1f',var2Plot{v},condLabel{i},mean(y(x>=0)));
%             elseif bitSumDim 
%                 % summing across dimensions, plotted as one signal, report mean
%                 legText{axN,nLine(axN)+1:nLine(axN)+numel(foo)} =  sprintf('%s sum(sRA 1,2,3), mean = %0.1f',var2Plot{v},mean(y(x>=0)));
%             else
%                 % if plotting multiple signals' RSV/ISV for each
%                 % dimensions, we have to add separate legend entries. 
%                 
%                 % loop through signals' RSV/ISV. When only plotting 1
%                 % signal, ind2=i.
%                 lineNTemp = 0;
%                 for ii = ind2
%                     lineNTemp = lineNTemp + 1;
%                     % if plotting sRA on separate axes, then we do not need to
%                     % include the sRA in the legend
%                     if bitSepAx
%                         % this should be gauranteed to have at most one
%                         % line at a time
%                         legText{axN,nLine(axN)+lineNTemp} =  sprintf('%s %s',var2Plot{v},condLabel{ii});
%                     else
%                         legText{axN,nLine(axN)+lineNTemp} =  sprintf('%s %s on sRA %s',var2Plot{v},condLabel{ii},condLabel{i});
%                     end
%                 end
%             end
%         end
                
        % update lines per axes, but only for the first pass through the
        % condition (signal) loop
        nLine(axN) = nLine(axN) + numel(foo);
        
        hold on;
        
        clear foo ind2

        yLim = get(gca,'YLim');
      
        % log scale for p-value
        if varType == 4
            set(axH(end),'YScale','log','YDir','reverse','YMinorTick','off')
%             yLim = get(gca,'YLim');
        end
        
    end
    
    % number of lines per axes
%     nLine(axN) = size(h,2);
    
    % loop through RAs
    for i = 1:nRA
        
        % if plotting on separate axes
        if bitSepAx
            
            axes(axH(i));
            
            axN = i;
        else
            axN = 1;
        end
        
        % intialize
        y = [];
        
        
        % add max variance, provided not plotting PCs or plotting PCs but
        % using dPC for max var)
        if (bitMaxVar && (~bitUsePC || bitUsePC4MaxVar)) ...
                || (bitRSVBound || bitRSVBoundDiff) %% && varType<3 % && ...
            % ~(any(ismember(var2Plot{v},{'RSV','ISV','ARSV'})) && varType==1)
            
            
            foo = varMaxFN;
            
            
            % for a reference comparison of maximum variance of a
            % single dimension, we generally use time-varying PC1.
            % However, for RSV and ISV this makes less sense since
            % the upper bound is not max variance, but max RELEVANT
            % variance, which should be estimated by the TDR(t)
            % vector. And yet, we are concerned that RSV of non-orthogonal
            % RAs is not accurate, so we will not enforce using RSV(dRA)
            if bitUseTDR4MaxVar % || ismember(var2Plot{v},{'RSV','ISV','ARSV'})
                
                
                % special case -- we don't have a TDR RSV term for absolute
                % variance (only variance explained)
%                 if varType == 1
%                     error('We don''t have a TDR RSV term for absolute variance (only variance explained). Set varType = 2 or 4')
%                 end
                
                foo = [foo,'_RA'];
                
                % take dim 1 (data times t), dim 2 (regression
                % times for dRA at time t), dim 3 (dRA for signal
                % i), and only when varType ==4, dim 4 (RSV/ISV for
                % signal i).
                yTemp = Summary.(tdr).varAnalysis.(foo);
                
                % set which dRA to take
                nRA2Use = i;
                
                % if DOVE analysis available, extract significance of
                % difference between dRA and sRA
                if isfield(Summary.(fn).(fnRA),'DOVE_sRAvdRA')
                    % get field
                    DOVEfn = var2Plot{v};
                    % note that DOVE is run on absolute variance, which
                    % does not change the p-value but requires that if
                    % var2Plot = 'VE', we must convert to 'V'
                    if strcmp(DOVEfn,'VE')
                        DOVEfn = 'V';
                    end
                    bitSigDiff = Summary.(fn).(fnRA).DOVE_sRAvdRA(i).(DOVEfn).p <= alpha;
                else
                    bitSigDiff = [];
                end
                
                % name for legend
                maxVarName = 'dRA';
            elseif bitUseSerialSRA4MaxVar 
                foo = [foo,'_RA'];
                
                % capture signal. NOTE THAT SERIAL SRA DOES NOT REQUIRE
                % ALIGNMENT! 
                goo = [':',num2cell(repmat(i,1,ndims(Summary.(otdrSerial).(fnRA).varAnalysis.(foo))-1))];
                y = Summary.(otdrSerial).(fnRA).varAnalysis.(foo)(goo{:});

                % if DOVE analysis available, extract significance of
                % difference between dRA and sRA
                if isfield(Summary.(fn).(fnRA),'DOVE_sRAvSerialSRA')
                    % get field
                    DOVEfn = var2Plot{v};
                    % note that DOVE is run on absolute variance, which
                    % does not change the p-value but requires that if
                    % var2Plot = 'VE', we must convert to 'V'
                    if strcmp(DOVEfn,'VE')
                        DOVEfn = 'V';
                    end
                    bitSigDiff = Summary.(fn).(fnRA).DOVE_sRAvdRA(i).(DOVEfn).p <= alpha;
                else
                    bitSigDiff = [];
                end
                
                % name for legend
                maxVarName = 'serial sRA';
                

            elseif bitUsePC4MaxVar                                
                % here we extract the V for a dynamic PC1 so as to
                % show the upperbound of how much a single dimension can
                % capture.

                foo = [foo,'_RA'];
                
                % take dim 1 (data times t), dim 2 (regression
                % times for dRA at time t), dim 3 (dRA for signal
                % i), and only when varType ==4, dim 4 (RSV/ISV for
                % signal i).
                yTemp = Summary.(tdr).varAnalysis_dPC.(foo);
                
                % set which dRA to take
                nRA2Use = 1;
                
                % dummy
                bitSigDiff = [];
                
                % name for legend
                maxVarName = 'dPC';                

%                 warning('Here we wish to plot V/RSV/ISV for a dynamic PC1 so as to show the upperbound')
                %                 if any(ismember(var2Plot{v},{'RSV','ISV','ARSV'}))
                %                     foo = [foo,num2str(i)];
                %                 end
                % %                   foo = [foo,'_dPC1'];
                %                     foo = [foo,'_RA1'];
                %
                %                 if eN == 1
                % %                     y = Summary.(fn).(fnVA).(foo);
                %                     y = Summary.(fn).PCs.(fnVA).(foo);
                %                 else
                % %                     y = Summary.(fn).PP.(fnVA).(foo);
                %                     y = Summary.(fn).PCs.(fnVA).(foo);
                %                 end
                
                %                     if ndims(Summary.(fn).PCs.(fnVA).(varFN)) == 2
                %                         % note -- here we use the 1st PC and find the RSV for
                %                         % the ith signal as projected onto it
                %                         y = Summary.(fn).PCs.(fnVA).(varFN)(:,1);
                %                     elseif ndims(Summary.(fn).PCs.(fnVA).(varFN)) == 3
                %                         % note -- here we use the 1st PC and find the RSV for
                %                         % the ith signal as projected onto it
                %                         y = Summary.(fn).PCs.(fnVA).(varFN)(:,1,i);
                %                     else
                %                         error('number of dimensions (%d) nor as expected',ndims(Summary.(fn).PCs.(fnVA).(varFN)))
                %                     end
                
%             else
%                 error('Unsure type of maximum variance to plot')
            end
            
            % extract variance values related where data times = regression
            % times. NOTE THAT SERIAL SRA DOES NOT REQUIRE ALIGNMENT!
            if bitUseTDR4MaxVar || bitUsePC4MaxVar
                y = TDRAlignDataAndRegressTimes(yTemp,nRA2Use,'SN',i,...
                    'T',Summary.(tdr).Times.all_times,'RT',Summary.(tdr).Times.regressTimes,...
                    'nBinRegress',Summary.(tdr).Times.regressBins);
            end
            
%             
%             % Applies to both dRA and dPC vectors:
%             % need to know how many data timebins per
%             % regression time bins. For now, assume the ratio
%             % between dim 1 and 2:
%             nBinRegress = size(yTemp,1) ./ size(yTemp,2);
%             if nBinRegress ~= round(nBinRegress)
%                 error('Not integer number of time bins per regress bins')
%             end
%             
%             y = NaN(size(yTemp,1),1);
%             % loop through regress bins to extract time bins
%             for ii = 1:size(yTemp,2)
%                 
%                 if ndims(yTemp) == 4
%                     % for RSV/ISV
%                     y(ii*nBinRegress-(nBinRegress-1):ii*nBinRegress) = ...
%                         yTemp(ii*nBinRegress-(nBinRegress-1):ii*nBinRegress,ii,nRA2Use,i);
%                 else
%                     % for VE
%                     y(ii*nBinRegress-(nBinRegress-1):ii*nBinRegress) = ...
%                         yTemp(ii*nBinRegress-(nBinRegress-1):ii*nBinRegress,ii,nRA2Use);
%                 end
%             end
%             
            
            % only plot on first pass through signals when plotting
            % absolute variance and variance explained (since they
            % share a common max variance) OR for all passes when
            % plotting RSV/ISV (since these are different for each
            % signal
            if bitMaxVar &&...
                    (i == 1 || bitUseTDR4MaxVar || bitUseSerialSRA4MaxVar || ismember(var2Plot{v},{'RSV','ISV','ARSV'}) || bitSepAx)...
                    && ~isempty(y) %  && ismember(varType,[1 2]))
                clear foo
                
                foo = plot(x,y);
                
                % line style depends on var
                switch var2Plot{v}
                    case 'RSV'
                        ls = ':';
                    case 'ISV'
                        ls = '-.';
                    otherwise
                        ls = ':';
                end
                
                set(foo,'LineWidth',2,'LineStyle',ls);
                
                if bitUseTDR4MaxVar || bitUseSerialSRA4MaxVar || ismember(var2Plot{v},{'RSV','ISV','ARSV'})
                    set(foo,'Color',sRAColor(Summary.(otdr).(fnRA).param4RA(i),:));
                else
                    set(foo,'Color',textColor);
                end
                
                %                         % update nLine
                %                         if i == 1
                %                             nLine(axN) = nLine(axN) +1;
                %                         end
                % store line handle
                h(axN,nLine(axN)+1:nLine(axN)+numel(foo)) = foo;
                
                % When plotting RSV/ISV, we have to specify for which
                % signals
                if ismember(var2Plot(v),{'RSV','ISV','ARSV'})
                    signal4Var = sprintf(' %s',condLabel{Summary.(otdr).(fnRA).param4RA(i)});
                else
                    signal4Var = '';
                end
                
                % update legend
                legText{axN,nLine(axN)+1:nLine(axN)+numel(foo)} =  ...
                    sprintf('Time-varying %s(%s)%s',var2Plot{v},maxVarName,signal4Var);
                
                % update lines per axes
                nLine(axN) = nLine(axN) + numel(foo);
                
                % add markers to those points that were significantly
                % different from sRA
                if ~isempty(bitSigDiff)
                    foo = plot(x(bitSigDiff),y(bitSigDiff),'*','MarkerSize',12);
                    if bitUseTDR4MaxVar || bitUseSerialSRA4MaxVar || ismember(var2Plot{v},{'RSV','ISV','ARSV'})
                        set(foo,'Color',sRAColor(Summary.(otdr).(fnRA).param4RA(i),:));
                    else
                        set(foo,'Color',textColor);
                    end
                end
                
                clear foo y
                
                %             error('for max variance for RSV and ISV, include TDR instead of PC1');
            end
            %             elseif varType<3 && i==1 && ~(any(ismember(var2Plot{v},{'RSV','ISV','ARSV'})) && varType==1)
            %                 warning('We are skipping plotting of PC1 for Deviation from Chanace analysis because these data are not available yet.');
            
            % Adding UPPER BOUND on off-target RSV given correlation
            % between predictors, and UPPER BOUND for on-target RSV given
            % variance explained
            if (bitRSVBound || bitRSVBoundDiff) && ismember(var2Plot(v),{'RSV','ISV','ARSV'}) ...
                    && ~bitMatchSignal2SRA ... 
                    && ~bitSumDim && ismember(varType,[1 2])
                
                
                %%% (I think this section is entirely moot, since for upper
                %%% bounds we force separate axes and prohibit PCs)
                % When on same figure axes we have to specify for which
                % signals
                if ~bitSepAx
                    if bitUsePC
                        if bitUsePC_CC
                            foo = 'CC ';
                        else
                            foo = '';
                        end
                        signal4Var = sprintf(' on %sPC_%d',foo,i);
                        clear foo
                    else
                        signal4Var = sprintf(' on sRA %s',condLabel{Summary.(otdr).(fnRA).param4RA(i)});
                    end
                else
                    signal4Var = '';
                end
                
                % compute and plot bound only for bitRSVBound (not
                % bitRSVBoundDiff)
                if bitRSVBound
                    
                    
                    % find index of off-target predictors
%                     ind = setdiff(1:length(condLabel),i);
                    
                    % loop through each offtarget predictor
                    for j = find(~bitTarget(i,:))
                        % compute correlation between on and off-target
                        % predictor
                        % compute correlation between on and off-target predictor.
                        % Must transform predictors in the same way as they were
                        % for computing RSV (see condSeparation.m), where they are
                        % mean-centered and normalized to unit vectors
                        a = Summary.(otdr).codedParams(:,bitTarget(i,:));
                        b = Summary.(otdr).codedParams(:,j);
                        a = a - mean(a);
                        b = b - mean(b);
                        a = a / norm(a);
                        b = b / norm(b);
                        
                        rij = corr(a,b);
%                         fprintf('rij = %f for i=%d and j=%d \n',rij,find(bitTarget(i,:)),j)
                        
                        % plot upper bound as RSV * rij^2
                        %                     foo = plot(x,Summary.(fn).(fnRA).(fnVA).(varFN)(:,i,i)*rij^2);

                        % extract variance explained by axis depending on absolute
                        % variance or normalized variance
                        if varType == 1
                            VE = Summary.(fn).(fnRA).(fnVA).V_RA(:,i);
                        elseif varType == 2
                            VE = Summary.(fn).(fnRA).(fnVA).VE_RA(:,i);
                        end
                        
                        % plot upper bound as VE * rij^2
                        foo = plot(x,VE*rij^2);
%                         foo = plot(x,Summary.(fn).(fnRA).(fnVA).VE_RA(:,i)*rij^2);
                        clear VE
                        
                        % apperance
                        set(foo,'LineWidth',2,'LineStyle',':','Color',sRAColor(j,:));
                        
                        % store line handle
                        h(axN,nLine(axN)+1:nLine(axN)+numel(foo)) = foo;
                        
                        
                        % update legend
                        legText{axN,nLine(axN)+1:nLine(axN)+numel(foo)} =  ...
                            sprintf('Upperbound %s %s%s due to corr r = %0.2g',var2Plot{v},condLabel{j},signal4Var,rij);
                        
                        % update lines per axes
                        nLine(axN) = nLine(axN) + numel(foo);
                        
                        clear foo a b rij
                    end
                end
                
                % add variance explained by axes as upper bound for
                % ontarget signal. If plotting RSV Diff (and not
                % proabability of RSV Diff), then this line should go on
                % the alternative axes created by plotyy
                if bitRSVBoundDiff && varType~=4
                    axTemp = axH2(i);
                else
                    axTemp = gca;
                end
                if varType == 1
                    foo = plot(axTemp,x,Summary.(fn).(fnRA).(fnVA).V_RA(:,i));
                elseif varType == 2
                    foo = plot(axTemp,x,Summary.(fn).(fnRA).(fnVA).VE_RA(:,i));
                end      
                clear axTemp
                
                % apperance
                set(foo,'LineWidth',2,'LineStyle',':','Color',PCColor);
                
                % store line handle
                h(axN,nLine(axN)+1:nLine(axN)+numel(foo)) = foo;
                
                % update legend
                legText{axN,nLine(axN)+1:nLine(axN)+numel(foo)} =  ...
                    sprintf('Variance explained%s (upper bound)',signal4Var);
                
                % update lines per axes
                nLine(axN) = nLine(axN) + numel(foo);
                
            end
        end
    end
    
    % loop through axes
    for axN = 1:length(axH)
        
        % call axis
        axes(axH(axN))
    
        % whether chance level goes to ylim(1) or (2) depends on whether
        % plotting VE or p-value
        yLim = get(gca,'YLim');        
        if varType == 1 || varType == 2
            goo = yLim(1);
        else
            goo = yLim(2);
        end
        
        % add chance level -- just do this for one type of variance -- they
        % will are share it:
        if v == 1
            hF = [];
            % does not makes sense to plot variable chance level for RSV and ISV since these
            % have a difference chance level for each signal. OK to plot a flat
            % chance level when plotting p(RSV) or p(ISV), or when plotting
            % RSV on different axes. Not OK to plot chance level for both
            % RSV and ISV unless plotting p(RSV) or p(ISV). Also, does not
            % make sense to plot chance level when summing across
            % dimensions, since chance level is for a single dimension.
            if ~isempty(chY) && ~(any(ismember(var2Plot(v),{'RSV','ISV','ARSV'})) && ...
                    ismember(varType,[1 2]) && ~bitSepAx) && ...
                    ~(all(ismember({'RSV','ISV','ARSV'},var2Plot)) && varType ~= 4) && ...
                    (eN == 1 || varType == 3 || varType == 4) && ... % For now, we don't have chance levels for the other epochs.
                    ~bitSumDim
                
%             if ~isempty(chY) && ~(any(ismember(var2Plot(v),{'RSV','ISV','ARSV'})) && ...
%                     ismember(varType,[1 2]) && ~bitUsePC && length(var2Plot)>1) && ...
%                     (eN == 1 || varType == 3 || varType == 4) % For now, we don't have chance levels for the other epochs.
                % choose which column of chance. always col 1 for varType 4
                % or when plotting Variance (not RSV, ISV)
                if varType == 4 || strcmp(var2Plot{v},'V')
                    col = 1;
                else
                    col = axN;
                end
                
                % DETERMINE APPEARANC OF P-VALUE THRESHOLD
                if bitLine4PThresh
                    % P-THRESH AS LINE                    
%                     hF = line([min(x), max(x)], ones(1,2) * unique(chY(:,col)), ...
%                         'Color', lineColor, 'LineStyle', '--','LineWidth',2);
                    % We used to have an actual line (see above), but more
                    % flexible to plot a curve, which can vary with time
                    % when varType~=4, or be time-invariant when
                    % varType==4.
                    hF = line(x, chY(:,col), ...
                        'Color', lineColor, 'LineStyle', '--','LineWidth',2);
                    
                else
                    % P-THRESH AS BOX 
                    
                    hF = fill([min(x);x;max(x)],[goo;chY(:,col);goo],[0.6 0.6 0.6]);

                    clear col
                    set(hF,'EdgeColor',[0.6 0.6 0.6]);

                    % when exporting to EPS, Illustrator, etc, turn the patch to
                    % transparent to prevent the patch from being broken-up
                    % into multiple triangles. Also set the edges to red so it
                    % can be distinguished easily
                    if bit4EPS
                        set(hF,'FaceAlpha',0,'EdgeColor','r','LineWidth',3);
                    end
                end
                
                % move p-value threshold behind other elements
                foo = get(gca,'Children');
                foo(foo==hF) = [];
                foo = [foo;hF];
                set(gca,'Children',foo);
                clear foo
                
                %%% NO LONGER ADD CHANCE LEVEL TO THE LEGEND
%                 % add hF to handles
%                 h(a,nLine(axN)+1) = hF;
%                 % add legend
%                 legText{a,nLine(axN)+1} = sprintf('Chance, p = %g',alpha);
%                 
%                 % update lines per axes, but only for the first pass through the
%                 % condition (signal) loop
%                     nLine(axN) = nLine(axN) + numel(hF);
            end
        end
        
        yLimAll(end+1,:) = yLim;
    
    end
    
    
end

% set yLim depending on varType
if varType == 4
    % probability scale
    foo = [min(yLimAll(:,1)) 1];
else
    foo = [min(yLimAll(:,1)) max(yLimAll(:,2))];
end

% at this point, all lines should be plotted, so we can set limits
set([axH axH2],'box','off','FontSize',16,'XLim',[min(x), max(x)])
% set ylim separately in case axH2 (plotyy)
set(axH,'YTickMode','auto','YLim',foo);
% loop through secondary axes and set ylim, ylim mode, and color
yLimAll2 = NaN(length(axH2),2);
for i = 1:length(axH2)
    yLimAll2(i,:) = get(axH2(i),'YLim');
end
for i = 1:length(axH2)
    set(axH2(i),'YLim',[-max(yLimAll2(:,2)) max(yLimAll2(:,2))],'YColor',lineColor);
    % change ticks to just bottom, 0 and top
    foo = get(axH2(i),'YLim');
    set(axH2(i),'YTick',[foo(1),0,foo(2)]);
end

% loop again through variables and axes to plot stuff that depends on limits
% Loop through different types of variance requested
for v = 1:length(var2Plot)

    % loop through axes
    for axN = 1:length(axH)

        % focus axes
        axes(axH(axN));
        
        % get limits
        yLim = get(gca,'YLim');
        
        % mark times at which vectors were selected. Only do this if
        % requested (i.e., bitTimeOfSRA==1) AND plotting sRAs (i.e., fn
        % starts with oTDR summary struct name, usually "oTDRSummary" AND
        % not plotting PCs)
        if bitTimeOfSRA && startsWith(fn,otdr) && ~bitUsePC
            
            if bitSepAx
                % if plotting on separate axes, just show the times for the
                % current signal
                
                sigSet = axN;
            else
                % loop through all signals
                sigSet = [1:length(Summary.(otdr).(fnRA).param4RA)];
            end
            
            for i = sigSet
                
                if varType == 4 % log scale
                    foo = 0.01;
                else
                    foo = 0.02;
                end
                
                if varType == 1 || varType == 2
                    y = yLim(2)-foo*(i-1)*range(yLim);
                else
                    y = 10^(log10(yLim(1))+(foo*(i-1)*-log10(yLim(1))));
                end
                
                % extract time at which sRA was computed
                if isfield( Summary.(fn).(fnRA),'t4RA')
                    sRA_t = Summary.(fn).(fnRA).t4RA{i};
                else                           
                    % old naming convention for field
                    sRA_t = Summary.(fn).(fnRA).(['t',num2str(i)]);
                end                
                
                % convert bin centers to bins edges
                t = [sRA_t-dt/2, sRA_t+dt/2]';
                plot(t,y*ones(size(t)),'Color',sRAColor(Summary.(otdr).(fnRA).param4RA(i),:),'LineWidth',2);
            end
        end
        
        % mark times of key task events
        if bitTimeOfEvent && isfield(Summary.(otdr).Times,'timeIntervals')
            
            % determine alignment field to use
            if strcmpi(alignName,'fixation')
                foo = 'Fix';
            else
                foo = alignName;
            end
            ff = ['align2',foo];
            
            % number of events
            eventNames = fieldnames(Summary.(otdr).Times.timeIntervals.(ff));
            nEvent = length(eventNames);
            % colors
            eventColor = lines(nEvent);
            
            % loop through events
            for i = 1:nEvent
                f = eventNames{i};
                
                % compute mean across neurons -- current using 2.5 and 97.5
                % percentiles
                eventBounds = mean(Summary.(otdr).Times.timeIntervals.(ff).(f).prctile(:,[1 end]));
                % event median
                eventMed = mean(Summary.(otdr).Times.timeIntervals.(ff).(f).prctile(:,Summary.(otdr).Times.timeIntervals.prctile2Use==50));
                
                if varType == 4 % log scale
                    foo = 0.01;
                else
                    foo = 0.02;
                end
                
                if varType == 1 || varType == 2
                    y = yLim(2)-foo*(i-1)*range(yLim);
                else
                    y = 10^(log10(yLim(1))+(foo*(i-1)*-log10(yLim(1))));
                end
                
                % draw line
                plot(eventBounds,y*ones(size(eventBounds)),'Color',eventColor(i,:),'LineWidth',2);
                % mark median (just in case line is very narrow)
                plot(eventMed,y,'Marker','s','Color',eventColor(i,:),'MarkerFaceColor',eventColor(i,:));
                % label event
                text(mean(eventBounds),y,f,'Color',eventColor(i,:),'FontSize',12,'HorizontalAlignment','center',...
                    'VerticalAlignment','bottom');
            end
        end
        
    end
    
end % type of variance 

% loop through axes, adding legend, lines, etc
for a = 1:length(axH)
    axes(axH(a));
    % find entries that should be skipped
    bitSkip = contains(legText(a,:),'SKIP');
    legend(h(a,~bitSkip),legText(a,~bitSkip),'box','off','Color','none',...
        'EdgeColor','k','FontSize',14);
    clear bitSkip
    
    % line at x = 0
    line([0 0],[min(yLimAll(:,1)) max(yLimAll(:,2))],'Color',lineColor);
    
    % when plotting difference between RSV and bound, include line at y=0
    if bitRSVBoundDiff
        line([min(x), max(x)],[0 0],'Color',lineColor);
    end
end

ylabel(axH(1),yLabel,'FontSize',18);
xlabel(axH(ceil(length(axH)/2)),xLabel,'FontSize',18);

%%% add annotation
dt = round(diff(Summary.(otdr).sRA.t4RA{1}(1:2)),2);
figTitleHand = annotation('textbox',[0 0 1 .05]);
if isfield(Summary.(otdr).sRA,'magThresh') && ~isempty(Summary.(otdr).sRA.magThresh)
    figAnnotStr = {sprintf('sRA vs. Time, present trial, sRA selected as R2 > %g',...
        Summary.(otdr).sRA.magThresh);
        [datestr(now),', ',mfilename,'.m']
        };    
else
    figAnnotStr = {sprintf('sRA vs. Time, present trial, sRA selected as: Benefit [%g %g], Choice [%g %g], Exp Rwd [%g %g]',...
        Summary.(otdr).sRA.t4RA{1}(1)-dt/2,Summary.(otdr).sRA.t4RA{1}(end)+dt/2,...
        Summary.(otdr).sRA.t4RA{2}(1)-dt/2,Summary.(otdr).sRA.t4RA{2}(end)+dt/2,...
        Summary.(otdr).sRA.t4RA{3}(1)-dt/2,Summary.(otdr).sRA.t4RA{3}(end)+dt/2);
        [datestr(now),', ',mfilename,'.m']
        };
end
set(figTitleHand,'String',figAnnotStr,'Interpreter','none','FontSize',12);
set(figTitleHand,'VerticalAlignment','top','HorizontalAlignment','right');
set(figTitleHand,'LineStyle','none');
set(figTitleHand,'Color','r');

end % end main function

%% SUBFUNCTION for cleaning up p = 0 values
function y = cleanPZero(y,varType)
    % only clean if varType == 4 (i.e., probabilities)
    if varType == 4
        % replace p = 0 values with eps value, unless there are non-zero
        % values smaller than zero, in which case throw an error so we can
        % sort this out.
        bitNonZero = y ~= 0;
        if any(y(bitNonZero) < eps) & any(~bitNonZero)
            error('Not expecting non-zero p-value less than EPS. Figure out new solution')
        end
        y(~bitNonZero) = eps;
    end
end
