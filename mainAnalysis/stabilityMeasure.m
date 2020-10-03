% [t0Best, wBest, aBest, errBest] = stabilityMeasure(thetaOrig, [name,value])
%
%% Inputs:
% - thetaOrig: a vector of length times x 1 of the angles between a 
%   reference dRA and all other dRAs.
% - Optional name, value pairs
%
%% Outputs:
% - t0Best: starting time of the step function.
% - wBest: duration of the step.
% - aBest: elevation of the step.
% - ErrorBest: the corresponding error
% Note, all the outputs are with respect to all the angles except the angle 
% at index tref
%
% OPTIONAL name, value pairs
% -aBound = 1 x 2 vector defining [lower upper] bounds of amplitude
% parameter aBest. Set to [-Inf Inf] to permit all values.
%
% -bitExclude = logical times x 1 vector specifying the elements of
% theraOrig that will be excluded from the fit.
%
% -bit2Step = logical to fit 2 sequential (non-overlapping) steps of
% independent height. Default = FALSE
%
% -bitPlot = logical on whether to plot figure.
%
% By Gamal Elsayed. 
% Modified by Daniel Kimmel, 8 Jul 2016, to accommodate multiple exclusions
% and implement bounds on step amplitude.
% Modified by Daniel Kimmel, 21 Jul 2016, to expand to searching for 2-step
% solutions.

function [t0Best, wBest, aBest, errBest] = stabilityMeasure(thetaOrig, varargin)

%% default values for optional input params

% define [lower upper] bounds of amplitude parameter A
aBound = [-Inf Inf];

% exclude specific entries
bitExclude = false(size(thetaOrig));

% logical to fit 2 sequential (non-overlapping) steps of independent height
bit2Step = false;

% plot?
bitPlot = false;

%% collect optionally provided parameters
warnopts(assignopts(who, varargin));

%% fit

theta = thetaOrig;
theta(bitExclude) = nan; % nan-out excluded points

% initialize
[t0Best, wBest, aBest, errBest] = deal(NaN);

%%% exit function...
if sum(isfinite(theta)) == 0
    % ...if no datapoints are available
    return
elseif all(theta(isfinite(theta))==0)
    % ...if all finite datapoints are zero (can't have a step!)
    return
end


T = length(theta);
theta = theta(:);              % make it column vector
normTerm = nansum(theta.^2);           % error normalization term

%% 1 vs. 2 step function

if ~bit2Step
    %%% 1-STEP
    minErperW = nan(T,1);               % initialize the minum error vector for each width
    t0BestperW = nan(T,1);              % initialize the starting index of minum error vector for each width
    parfor w = 1:T
        minEr = inf;
        t0Best = NaN; % must be defined given parfor and continue statement
        % set minimum error to highest possible value
        for t0 = 1:T-w+1
            % skip if value of theta is excluded at start point or end
            % point
            if isnan(theta(t0)) || isnan(theta(t0+w-1))
                continue
            end
            stepFn = zeros(T,1);
            stepFn(t0:t0+w-1) = min(aBound(2),max(aBound(1),nanmean(theta(t0:t0+w-1)))); % create a step function, while implementing bounds
            Er = nansum((stepFn - theta).^2)/normTerm;  % measure the error with respect to the step function
            if Er<minEr
                minEr = Er;              % if error is lower than a previous min error set the new min error to be current value
                t0Best = t0;             % set the best step index to be current step value
            end
        end
        minErperW(w) = minEr;           % store min error
        t0BestperW(w) = t0Best;         % store step index of the min error
    end
    
    [errBest, wBest] = min(minErperW);        % calculate the step width that yields the smallest error
    t0Best = t0BestperW(wBest);         % obtain the corresponding step index
    aBest = min(aBound(2),max(aBound(1),nanmean(theta(t0Best:t0Best+wBest-1)))); % obtain the elevation of the step, applying bounds
else
    %%% 2-STEP
    % need at least two finite values of theta
    if sum(isfinite(theta)) > 1
        minErperW1 = nan(T,1);               % initialize the minum error vector for each width
        t01BestperW1 = nan(T,1);              % initialize the starting index of minum error vector for each width
        w2BestperW1 = nan(T,1);              % initialize the starting index of minum error vector for each width
        t02BestperW1 = nan(T,1);              % initialize the starting index of minum error vector for each width
        parfor w1 = 1:T
            minEr = inf;         % set minimum error to highest possible value
            t01Best = NaN; % must be defined given parfor and continue statement
            w2Best = NaN; % must be defined given parfor and continue statement
            t02Best = NaN; % must be defined given parfor and continue statement
            for t01 = 1:T-w1+1
                % skip if value of theta is excluded at start point or end
                % point
                if isnan(theta(t01)) || isnan(theta(t01+w1-1))
                    continue
                end
                
                % loop through 2nd step -- assumes non-overlapping with first
                % step
                for t02 = t01+w1:T
                    
                    % skip if value of theta is excluded at start point
                    if isnan(theta(t02)) % || isnan(theta(t02+w2-1))
                        continue
                    end
                    
                    for w2 = 1:T-t02+1
                        
                        % skip if value of theta is excluded at endpoint
                        if isnan(theta(t02+w2-1))
                            continue
                        end
                        
                        % create step fn
                        [stepFn] = stepFn_2step(theta,t01,w1,t02,w2,aBound,T)
                        Er = nansum((stepFn - theta).^2)/normTerm;  % measure the error with respect to the step function
                        if Er<minEr
                            minEr = Er;              % if error is lower than a previous min error set the new min error to be current value
                            t01Best = t01;
                            w2Best = w2; % must be defined given parfor and continue statement
                            t02Best = t02; % must be defined given parfor and continue statement
                        end
                    end
                end
            end
            minErperW1(w1) = minEr;           % store min error
            t01BestperW1(w1) = t01Best;
            w2BestperW1(w1) = w2Best;
            t02BestperW1(w1) = t02Best;
            
        end
        
        [errBest, wBest] = min(minErperW1);        % calculate the step width that yields the smallest error
        t0Best = t01BestperW1(wBest);         % obtain the corresponding step index
        wBest(2) = w2BestperW1(wBest);
        t0Best(2) = t02BestperW1(wBest(1));
        [~,aBest] = stepFn_2step(theta,t0Best(1),wBest(1),t0Best(2),wBest(2),aBound,T);
        
    else
        [t0Best(1:2), wBest(1:2), aBest(1:2), errBest] = deal(NaN);

    end
end

%% plot summary figure
if bitPlot
    x = 1:T;
    stepFnBest = zeros(T,1);
    stepFnBest(t0Best:t0Best+wBest-1) = aBest;
    figure;
    hold on;
    plot(x,thetaOrig, 'o', 'color', 0.5*[1 1 1]);
    plot(x(~bitExclude),thetaOrig(~bitExclude),'ko')
    plot(x(~bitExclude), stepFnBest(~bitExclude), 'r')
%     ylim([0 90])
%     ylabel('90-theta')
end

end % main function

%% SUB FUNCTION for creasing 2 step function

    function [stepFn,aBest] = stepFn_2step(theta,t01,w1,t02,w2,aBound,T)
        aBest(1) = min(aBound(2),max(aBound(1),nanmean(theta(t01:t01+w1-1))));
        aBest(2) = min(aBound(2),max(aBound(1),nanmean(theta(t02:t02+w2-1))));
        stepFn1 = zeros(T,1);
        stepFn1(t01:t01+w1-1) = aBest(1); % create a step function, while implementing bounds
        stepFn2 = zeros(T,1);
        stepFn2(t02:t02+w2-1) = aBest(2); % create a step function, while implementing bounds
        stepFn = stepFn1 + stepFn2;

    end % nested function

