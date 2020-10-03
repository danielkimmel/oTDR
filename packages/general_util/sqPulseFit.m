function [cfit,gof] = sqPulseFit(x,data,varargin)
%
% [cfit,gof] = sqPulseFit(x,data,varargin)
%
% Fits continuous approximation of square pulse to DATA using difference of
% logistics:
% DATA = a * (logistic(x,shape,b) - logistic(x,shape,b+c)), 
% where logistic(x,shape,p) is defined as:
% Y = 1 / (1 + exp(-shape*(x-p)))
% The shape parameter is hard-coded to be very steep as dx * 1000.
% Note that by modeling the second term horizontal shift parameter as
% (b+c), we define c as a "step width" parameter relative to b, and allow
% the opportunity to contrain c > 0, so that step has to increase from
% baseline (otherwise, if b+c < b, then step will decrease from baseline).
% The downside of this is that b+c cannot be constrained to be less than
% max(x).
%
% ACCEPTS
% X -- vector of indepedent variable, e.g., time
% DATA -- vector of dependent variable
%
% OPTIONAL [name,value] pairs
%
% RETURNS
% CFIT object
% GOF goodness-of-fit structure
%
% Daniel Kimmel, 7 Jun 2016

%% define parameters

% logical on fit formula. We can either constrain the horizontal order of
% two terms (TRUE), in which case "c" refers to the width:
%     a*(1/(1+exp(-%f*(x-b))) - 1/(1+exp(-%f*(x-(b+c)))))
% Or we can constrain the endpoint of the pulse (FALSE), in which case "c"
% refers to the endpoint:
%     a*(1/(1+exp(-%f*(x-b))) - 1/(1+exp(-%f*(x-(c)))))
% In the future, if we really cared, we could use FMINCON() to constrain
% both "c" and "b+c".
bitConstrainOrder = false;

% number of iterations to repeat fit so as to estimate best fit
nIter = 100;

% define [lower upper] bounds of amplitude parameter A
aBound = [-Inf Inf];

% % define range to exclude. data OUTSIDE these bounds will be
% % considered outliers and excluded
% dataRange2Incl = [-Inf Inf];

% exclude specific entries
bitExclude = false(size(data));

%% collect optionally provided parameters
warnopts(assignopts(who, varargin));

%% clean up and checks

% transpose to column vars if needed
if size(x,1) < size(x,2)
    x = x';
end
if size(data,1) < size(data,2)
    data = data';
end
if size(bitExclude,1) < size(bitExclude,2)
    bitExclude = bitExclude';
end
% % data range should be ascending
% if diff(dataRange2Incl) < 0
%     error('dataRange2Incl must be increasing')
% end

%% fit

% approximate shape parameter by fraction of dx
dx = x(2) - x(1);
shape = dx * 1000;

% define model: double logistic as continuous approximation to square pulse
if bitConstrainOrder
    eq = sprintf('a*(1/(1+exp(-%f*(x-b))) - 1/(1+exp(-%f*(x-(b+c)))))',shape,shape);
   % compute upper and lower bounds (amplitude, start, width)
    lb = [aBound(1) min(x)-dx/2 0];
    ub = [aBound(2) max(x)+dx/2 max(x)+dx/2];
else
    eq = sprintf('a*(1/(1+exp(-%f*(x-b))) - 1/(1+exp(-%f*(x-c))))',shape,shape);
   % compute upper and lower bounds (amplitude, start, stop)
    lb = [aBound(1) min(x)-dx/2 min(x)-dx/2];
    ub = [aBound(2) max(x)+dx/2 max(x)+dx/2];
%     lb = [aBound(1) min(x)-2*dx min(x)-2*dx];
%     ub = [aBound(2) max(x)+2*dx max(x)+2*dx];
end

% eliminate imprecision
lb = roundDec(lb);
ub = roundDec(ub);

% outliers
% bitExclude = bitExclude | roundDec(data,5)<dataRange2Incl(1) | ...
%     roundDec(data,5)>dataRange2Incl(2);

ft = fittype(eq,'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a','b','c'});

cfitTemp = cell(nIter);
gofTemp = cell(nIter);
sse = NaN(nIter,1);


% only fit if there are at least 3 points
if sum(isfinite((data(~bitExclude)))) >= 3
    parfor i = 1:nIter
        sp = rand(1,3); % start points
        % scale by data range
    %     sp(1) = sp(1) * range(data(~bitExclude))+min(data(~bitExclude));
    %     sp([2 3]) = sp([2 3]) * range(x) + min(x);

        % fit options
        fo = fitoptions('method','NonlinearLeastSquares',...
            'Lower',lb,...
            'Upper',ub,...
            'Startpoint',sp,...
            'exclude',bitExclude);

        [cfitTemp{i},gofTemp{i}] = fit(x,data,ft,fo);
        sse(i) = gofTemp{i}.sse;
    end

    % best fit
    [~,pos] = min(sse);
    cfit = cfitTemp{pos};
    gof = gofTemp{pos};

else
    cfit = NaN;
    gof = NaN;
    return
end

% store info about gof
sseTemp = NaN(nIter,1);
rsquareTemp = NaN(nIter,1);
adjrsquareTemp = NaN(nIter,1);
for i = 1:nIter
    sseTemp(i) = gofTemp{i}.sse;
    rsquareTemp(i) = gofTemp{i}.rsquare;
    adjrsquareTemp(i) = gofTemp{i}.adjrsquare;
end

% when boundary conditions are hit, r values can be negative. Eliminate
% these:
rsquareTemp(rsquareTemp < 0) = NaN;
adjrsquareTemp(adjrsquareTemp < 0) = NaN;

gof.sse_var = var(sseTemp);
gof.sse_max = max(sseTemp);
gof.rsquare_var = nanvar(rsquareTemp);
gof.rsquare_min = nanmin(rsquareTemp);
gof.adjrsquare_var = nanvar(adjrsquareTemp);
gof.adjrsquare_min = nanmin(adjrsquareTemp);


% figure;
% stem(x,data);
% hold on;
% plot(x,cfit(x),'g');