function [out,fo] = stabilityFitStep(theta,t,varargin)
% [out,fo] = stabilityFitStep(theta,t,[name,value pairs])
% TODO: update docstring
%
% Function to assess stability of dRA representations from TDR based on
% STEP FUNCTION fits to the angle between dRAs (THETA -- MxM matrix, where
% M is the number of time points) defined at time T (vector of time stamps
% for elements in THETA).
%
% Returns :
% OUT -- struct with extracted times T, angles THETA, fit params B, and fit
% param 95% confidence intervals BINT, where each field is a matrix with
% trial times in rows. Optionally includes fit parameters to the surrogate
% data BS, when surrogate data provided (see optional thetaS).
% FO -- struct where each element has a fit object for each timepoint
%
% Compare to CB_stability_decay, which fits expotential decay functions.
%
% Daniel Kimmel, 2016 Jun 07

%% define parameters

% specify fit method. Can use exhaustive parameter search ('paramSearch')
% via stabilityMeasure(), or gradient descent ('gradDesc') via
% sqPulseFit().
fitMethod = 'paramSearch';

% logical to fit 2 sequential (non-overlapping) steps of
% independent height. Default = FALSE
bit2Step = false; 

% logical on fit formula. We can either constrain the horizontal order of
% two terms (TRUE), in which case "c" refers to the width:
%     a*(1/(1+exp(-%f*(x-b))) - 1/(1+exp(-%f*(x-(b+c)))))
% Or we can constrain the endpoint of the pulse (FALSE), in which case "c"
% refers to the endpoint:
%     a*(1/(1+exp(-%f*(x-b))) - 1/(1+exp(-%f*(x-(c)))))
% In the future, if we really cared, we could use FMINCON() to constrain
% both "c" and "b+c".
% Irrelevant in case of fitMethod = 'paramSearch'
bitConstrainOrder = false;

% logical on whether to remove peak angle (where t = time of dRA)
bitRemovePeak = true;

% set angle offset. This is the amount by which theta is offset as theta* =
% angleOffset - theta. As such, identical angles become angleOffset. When
% set to 90deg, this means that theta* = 0deg refer to orthogonal angles.
% If max(theta) > 90, you may wish to increase angleOffset. However, this
% changes the interpretation of theta* = 0deg to refer to angles that are
% more similar but OPPOSITE.
% Set angleOffset to EMPTY to use the max angle for any given row of theta
% as angleOffset. Note in this case, the offset may be <90deg.
angleOffset = 90;

% set angle offset minimum. In the event that angleOffset is empty, the
% offset will be either the max angle for a given row of theta or
% angleOffsetMin, whichever is greater. Set angleOffset min to -Inf to use
% the max angle regardless. (If left empty, will be set to -Inf)
angleOffsetMin = 90;

% logical on whether to reflect angles about 90deg. When true, theta will
% be transformed to theta-2*max(theta-90,0).
bitReflectTheta = false;

% set scalar to apply to angles such that theta* = angleScale*(angleOffset
% - theta). Generally, scalar will be 1 (default) or -1 to invert the
% angles. This is useful when attempting to fit angles > 90deg.
angleScale = 1;

% define range to include, inclusive. data OUTSIDE these bounds will be
% considered outliers and excluded. Refer to theta in original coordinates,
% before any transformation applied below.
dataRange2Incl = [-Inf Inf];

% define time range to include in fits as [min max] inclusive. Reference t
% input.
tRange2Incl = [0 Inf];

% surrogate data set of theta (dRA time x eval time x repetition)
thetaS = [];

% number of repetitions of fit before determining best fit
nIter = 100;

% Goodness of fit parameter to extract from surrogate fits
gof2Use = 'adjrsquare'; % 'adjrsquare', 'rsquare'

%% collect optionally provided parameters
warnopts(assignopts(who, varargin));


%% basic vars
dt = t(2)-t(1);
% timeRange = t(end) - t(1) + dt;

if isempty(angleOffsetMin)
    angleOffsetMin = -Inf;
end

% determine if angles are folded about 90 deg (acos(abs(dotProd))) or span
% a range > 90deg (acos(dotProd)).
% Currently, we do not use this flag. We could implement it below to force
% angleOffset to be 90deg when all angles are <= 90deg, but this is in
% effect implemented by the default value of angleOffset above -- it would
% only be in rare instances where angleOffset was left empty that this
% would be violated.
bitAngleFolded = max(theta(:)) <= 90;
if bitAngleFolded && isempty(angleOffset) && isinf(angleOffsetMin)
    warning('ALL THETA < 90 DEG BUT ANGLE OFFSET NOT SPECIFIED. This will allow angle offset to be <90 and vary for each row of theta.')
end


%% fit curve

% instantiate
fo = struct([]);

% Store vars for output
out.dataRange2Incl = dataRange2Incl;
out.bitRemovePeak = bitRemovePeak;
out.angleScale = angleScale;
out.aBound = NaN([size(theta,1),2]);
out.bitReflectTheta = bitReflectTheta;

out.tMaster = t; % store so that output has originals
out.thetaMaster = theta; % store so that output has originals
out.t = NaN(size(theta));
out.theta = NaN(size(theta));
out.gof = NaN(size(theta,1),1);
out.gof_var = NaN(size(theta,1),1);
out.gof_min = NaN(size(theta,1),1);
if strcmp('paramSearch',fitMethod) && bit2Step
    out.b = NaN(size(theta,1),3,2);
else
    out.b = NaN(size(theta,1),3);
    out.bint = NaN([size(out.b),2]);
end

% transform theta such that values > 90 are reflected about 90deg.
if bitReflectTheta
    theta = theta - 2*max(theta-90,0);
end

% if bitRemovePeak
%     out.dataRange2Incl = [-Inf angleOffset];
% else
%     out.dataRange2Incl = [-Inf Inf];
% end
posFit = [];

bitT2Incl = t >= tRange2Incl(1) & t <= tRange2Incl(2);

for i = 1:size(theta,1)
    
    % skip if time point out of range
    if ~bitT2Incl(i)
        continue
    end

    % excerpt data to fit
    x = theta(i,:);

    % instantiate bitmask of excluded elements
    bitExclude = false(size(x));
    
    % exclude timepoints
    bitExclude = bitExclude | ~bitT2Incl;    
%     % set excluded time points to NaN
%     x(~bitT2Incl) = NaN;
    
    % remove peak
    if bitRemovePeak
        bitExclude(i) = true;
    end
    
    % determine offset if not specified
    if isempty(angleOffset)
        aoTemp = max(max(x),angleOffsetMin);
    else
        aoTemp = angleOffset;
    end
    
    % invert by angle offset and scale by angle scale
    x = out.angleScale * (aoTemp - x);
    out.aBound(i,:) = [0 aoTemp];
    % transform data range, and reorder so that it goes from min to max
    d2iTemp = sort(out.angleScale * (aoTemp - out.dataRange2Incl));

    % apply data range
    bitExclude = bitExclude | roundDec(x)<d2iTemp(1) | ...
        roundDec(x)>d2iTemp(2); 
    
%     % depending on whether angles are folded about 90
%     if bitAngleFolded 
%         % invert so that max is 90
%         x = angleOffset - x;
%         out.aBound(i,:) = [0 angleOffset];
%     else
%         % invert so that max is max angle or 90, whichever is greater
%         foo = max(max(x,angleOffset));
%         x = foo - x;
%         out.aBound(i,:) = [0 foo];
%         clear foo
%     end
    
    % FIT
    switch fitMethod
        case 'paramSearch'
            [start, width, out.b(i,1,:), out.gof(i)] = stabilityMeasure(x, ...
                'aBound',out.aBound(i,:),'bitExclude',bitExclude,...
                'bit2Step',bit2Step);
            % convert start position and width bin count into times:
            if ~any(isnan(start))
                out.b(i,2,:) = t(start);
%             else 
%                 out.b(i,2,:) = NaN;
            end
            if ~any(isnan(width))
                out.b(i,3,:) = width * dt;
%             else
%                 out.b(i,3) = NaN;
            end
        case 'gradDesc'
            [fo(i).f, fo(i).gof] = sqPulseFit(t,x,'aBound',out.aBound(i,:),'nIter',nIter,...
                'bitExclude',bitExclude,'bitConstrainOrder',bitConstrainOrder);
            
            % save fit
            %     fo(i).gofVar = var(gofTemp);
            %     fo(i).gofMin = min(gofTemp);
            fo(i).t = t(~bitExclude);
            fo(i).theta = x(~bitExclude);
            
        otherwise
            error('Fit method %s not recognized',fitMethod)
    end
    
    tTemp = t;
    tTemp(bitExclude) = NaN;
    xTemp = x;
    xTemp(bitExclude) = NaN;

    out.t(i,:) = tTemp;
    out.theta(i,:) = xTemp;
    
    if strcmp(fitMethod,'gradDesc')
        % save lambda params
        if isobject(fo(i).f)
            out.b(i,:) = coeffvalues(fo(i).f);
            if sum(isfinite(xTemp)) - sum(isfinite(out.b(i,:))) > 0
                out.bint(i,:,:) = confint(fo(i).f)';
            end
            % store position of an actual fit object
            posFit = i;
        end
        if isstruct(fo(i).gof)
            out.gof(i) = fo(i).gof.(gof2Use);
            out.gof_var(i) = fo(i).gof.(sprintf('%s_var',gof2Use));
            out.gof_min(i) = fo(i).gof.(sprintf('%s_min',gof2Use));
        end
    end
end

% store for output
out.fitTypeName = 'SquarePulse';
out.fitMethod = fitMethod;
if bitConstrainOrder || strcmp(fitMethod,'paramSearch')
    out.coeffNames = {'Amplitude','Start','Width'};
else
    out.coeffNames = {'Amplitude','Start','Stop'};
end
if strcmp(fitMethod,'paramSearch')
    out.bitConstrainOrder = NaN;
    out.nIter = NaN;
else
    out.bitConstrainOrder = bitConstrainOrder;
    out.nIter = nIter;
end
if ~isempty(posFit)
    out.fitType = type(fo(posFit).f);
    out.fitFormula = formula(fo(posFit).f);
else
    out.fitType = NaN;
    out.fitFormula = NaN;
end
out.bitRemovePeak = bitRemovePeak;
if ~isempty(angleOffset)
    out.angleOffset = angleOffset;
else
    out.angleOffset = 0;
end

%% SURROGATE FITS

if ~isempty(thetaS)
    
    % instantiate
    bS = NaN([size(out.b),size(thetaS,3)]);
    gofS = NaN(size(out.b,1),size(thetaS,3));
    gofS_var = NaN(size(out.b,1),size(thetaS,3)); % variance of GOF across iterations
    gofS_min = NaN(size(out.b,1),size(thetaS,3)); % Minimum of (worse) GOF across iterations
    
    % loop through surrogates
    for s = 1:size(thetaS,3)
        
        % the theta matrix for a given surrogate dataset may have an entire
        % row (or column) where theta=90 deg. While this could, in theory,
        % occur because of a truly orthogonal dRA (that is orthogonal to
        % all the other dRA), it is most likely due to a time point for
        % which the dRA magnitude was 0 (and thus all angles with that dRA
        % are 90 deg). However, we don't have access to the surrogate dRAs
        % at this level (that would be in
        % TDRSummary.angleAnalysis.surrdRA[n]), and so we can't check the
        % dRA magnitude. Instead we rely on an entire row/column of 90-deg
        % angles.
        if any(all(thetaS(:,:,s)==90,1)) || any(all(thetaS(:,:,s)==90,2))
            warning('Skipping angle profile fits to surrogate dataset %d because at least 1 dRA had angles of 90 deg with ALL other dRAs, likely because it was zero magnitude.',s);
            continue
        end
        
        %%% RECURRANT CALL, including same parameters, except no
        %%% plotting and no surrogates.
        goo = reshape(varargin,2,numel(varargin)/2);
%         disp(sprintf('%d fits to theta profile of surrogate data',s));
        [foo] = stabilityFitStep(squeeze(thetaS(:,:,s)),t,goo{:},'thetaS',[]);
        goo = [];
        
        % extract parameters
        bS(:,:,s) = foo.b;

        gofS(:,s) = foo.gof;
        gofS_var(:,s) = foo.gof_var;
        gofS_min(:,s) = foo.gof_min;
        foo = [];
        
    end
    
    out.bS = bS;
    out.gofS = gofS;
    out.gofS_var = gofS_var;
    out.gofS_min = gofS_min;
    clear bS gofS
end

%% clean up -- only needed for gradient descent method

if strcmp(fitMethod,'gradDesc')
    [st,i] = dbstack;

    % only perform following functions after all interative calls have been
    % made
    if length(st) == 1 || ~strcmp(st(i).name,st(i+1).name)

            % constrain start parameter to be >= min time
            bitOOB = out.b(:,2) < min(t)-dt/2;
            if any(bitOOB)
                warning('LIMITING ''START'' parameter to minimum time for %d rows of THETA',sum(bitOOB));
                out.b(bitOOB,2) = roundDec(max(out.b(bitOOB,2),min(t)-dt/2));
                out.bint(bitOOB,2,:) = NaN;
            end
            if isfield(out,'bS')
                out.bS(:,2,:) = roundDec(max(out.bS(:,2,:),min(t)-dt/2)); 
            end

            % if NOT constraining term order, then 3rd param is the STOP
            % position, which we transform to WIDTH parameter
            if ~out.bitConstrainOrder

                % At times, the stop time will be BEFORE the start time. This
                % only happens when order is not constrained and results in a
                % downward deflecting step. For now, we eliminate these
                % entries.
                bitOOO = out.b(:,3) < out.b(:,2);
                if any(bitOOO)
                    warning('STOP time was BEFORE START time for %d rows of THETA, which will be eliminated',sum(bitOOO));
                    out.b(bitOOO,:) = NaN;
                    out.bint(bitOOO,:,:) = NaN;
                    out.gof(bitOOO) = NaN;
                    out.gof_var(bitOOO) = NaN;
                    out.gof_min(bitOOO) = NaN;
                end                
                % likewise, check the same for the surrogates
                if isfield(out,'bS')
                    bitOOO = out.bS(:,3,:) < out.bS(:,2,:);
                    if any(bitOOO(:))
                        warning('In Surrogates, STOP time was BEFORE START time for %d rows of THETA, which will be eliminated',sum(bitOOO(:)));
                        out.bS(repmat(bitOOO,1,size(out.bS,2),1)) = NaN; 
                        out.gofS(squeeze(bitOOO)) = NaN; 
                        out.gofS_var(squeeze(bitOOO)) = NaN; 
                        out.gofS_min(squeeze(bitOOO)) = NaN; 
                    end
                end 

                % replace stop parameter with width parameter:
                % store the stop parameter and it's CI if not already nan from a prior run
                warning('REPLACING ''STOP'' parameter in angle profile fits with WIDTH.')
                out.stopParam = out.b(:,3);
                out.stopParamInt = squeeze(out.bint(:,3,:));
                % compute width, while pegging the start and stop parameters at
                % the min or max times
                out.b(:,3) = roundDec(min(out.b(:,3),max(t)+dt/2) ...
                    - out.b(:,2));
                out.bint(:,3,:) = NaN;
                out.coeffNames{3} = 'Width';
                if isfield(out,'bS')
                    out.stopParamS = squeeze(out.bS(:,3,:));
                    out.bS(:,3,:) = roundDec(min(out.bS(:,3,:),max(t)+dt/2) ...
                        - out.bS(:,2,:));
                end
            else
                % if constraining term order, then 3rd param is the WIDTH
                % parameter, which we limit
                bitOOB = out.b(:,3) + out.b(:,2) > max(t)+dt/2;
                if any(bitOOB)
                    warning('LIMITING ''WIDTH'' parameter such that start + width <= max time for %d rows of THETA',sum(bitOOB));

    %                 % store original width
    %                 out.widthParamOrig = out.b(:,3);
    %                 out.widthParamIntOrig = squeeze(out.bint(:,3,:));

                    out.b(bitOOB,3) = roundDec(min(out.b(bitOOB,3), ...
                        max(t)+dt/2 - out.b(bitOOB,2)));
                    out.bint(bitOOB,3,:) = NaN;
                end
                if isfield(out,'bS')
                    out.bS(:,3,:) = roundDec(min(out.bS(:,3,:),...
                        max(t)+dt/2 - out.bS(:,2,:)));
                end
            end
    end
end