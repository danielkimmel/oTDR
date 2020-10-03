function [out] = DOVE(v1,v2,rv,varargin)
% [out] = DoVE(v1,v2,rv)
%
% Computes the difference of varianace explained (DoVE) between variance
% explained values in V1 and V2 and computes statistic of this difference
% based on the distribution of chance (or random) differences computed from
% RV, the distribution of chance variance. Units of variance explained
% (absolute difference, percent explained, relevant/irrelevant variance)
% are not important provided the same units are used across V1, V2, and RV.
%
% ACCEPTS
% v1, v2 -- N x 1 vectors of variance explained. 
% rv -- N x S vector of S chance variance values for each observation N.
%       For instance, the values returned by oTDR() -->
%       oTDRSummary.sRA.varAnalysis.varRand_t.
%
% RETURNS
% out -- structure with following fields:
%   .vDiff -- N x 1 vector of differences v1 - v2
%   .p -- N x 1 vector of p-values for vDiff
%   .nullDistModel -- string with name of model fit to random (null)
%       distribution
%   .param -- N x P matrix of parameter values associated with model fit to
%       random (null) distribution
%   .ciSet -- M x 2 matrix of lower (col 1) and upper (col 2) confidence
%       intervals (in terms of percent, 0 to 100) at which 
%       distribution of random differences was evaluted for field .ci
%   .ci -- M x 2 matrix of values of vDiff associated with confidence
%       intervals specified in .ciSet
%
% Daniel Kimmel, 12 Dec 2016

%% Default param values

% string specifying the model to use to fit to random differences of
% variance. Options include 'gamma, 'gauss', 'empirical' (no model fit)
rpdf = 'gamma';

% N x S matrix of random difference which prevents having to compute this
% matrix based on rv. When rd is provided, rv is not necessary.
rd = [];

% M x 2 matrix of lower (col 1) and upper (col 2) confidence intervals (in
% terms of percent, 0 to 100) at which to evaluate the distribution of
% random difference.
ciSet = [0.05 99.5; 2.5 97.5; 5 95; 25 75];

%% collect optionally provided parameters
warnopts(assignopts(who, varargin));

%% check vars and orient

if sum(size(v1) > 1) > 1 || sum(size(v2) > 1) > 1
    error('v1 and v2 must be vectors, not matricies')
end

% orient v1, v2 as column vector
v1 = v1(:);
v2 = v2(:);

% Number of values of variance
N = length(v1);

if length(v2) ~= N
    error('Must provide same number of observations in V1 and V2')
end

% check that random matrix is in correct orientation
if isempty(rd) && size(rv,1) ~= N
    error('Must provide N vectors of random variance for each of N values provided in v1 and v2')
end

% check that random difference matrix, if provided is correct size
if ~isempty(rd) && size(rd,1) ~= N
    error('When provided matrix of random differences, it must have N vectors of random difference for each of N values provided in v1 and v2')
end


%% compute difference of variance
dv = v1 - v2;

%% compute distribution of random difference


% initialize matrix of random differences
% rd = NaN(N,D);

% initialize matrix of parameter fits, pvalue, and confidence intervals
switch rpdf
    case 'gamma'
        P = 2;
    case 'empirical'
        P = 0;
    otherwise
        error('Model name "%s" to fit to random differences of variance not recognized',rpdf);
end
prm = NaN(N,P);
p = NaN(N,1);
ci = NaN(N,size(ciSet,1),2);

% prepare to compute random differences if not provided
if isempty(rd)
    % Number of random values
    S = size(rv,2);
    
    % compute number of pairwise differences
    D = (S^2 - S)/2;
    
    % pairwise differences will be returned as matrix, and we want to extract
    % the upper triangle, so build index to extract it
    ind = triu(true(S),1);
    
    % because holding var for distribution of random differences is so
    % large, try dimensioning it first:
%     foo = NaN(S);
%     rdTemp = NaN(D,1);
else
    D = size(rd,2);
    
    % make sure differences are in absolute value
    rd = abs(rd);
end



% loop through each row in rv
% When having to compute differences for multiple vectors, use a parralel
% for loop
if isempty(rd) && N > 2
    parfor i = 1:N
        
        % call function. Indexing of rd accommodates empty rd
        [p(i),prm(i,:),ci(i,:,:)] = ...
            fitRD(rd(max(i,isfinite(rd)),:),rv(i,:),dv(i),ciSet,ind,rpdf);
        
%         % if not provided, compute all pairwise differences between columns as
%         % ABSOLUTE difference
%         foo = bsxfun(@minus,rv(i,:),rv(i,:)');
%         rdTemp = abs(foo(ind));
%         
%         % fit chance distribution and compute p value
%         [p(i),prmTemp] = sigTest(abs(dv(i)),rdTemp,'upper',rpdf);
%         % for non-empirical PDFs
%         if ~strcmp(rpdf,'empirical')
%             prm(i,:) = prmTemp;
%             prmTemp = num2cell(prmTemp); % to allow passing elements as separate arguments
%             ci(i,:,:) = icdf(rpdf,ciSet/100,prmTemp{:});
%         else
%             % store as temp var since prctile can only accept a vector of
%             % CI's
%             ciTemp = prctile(rdTemp,ciSet(:));
%             ci(i,:,:) = reshape(ciTemp,length(ciTemp)/2,2);
%         end      
    end
else
    for i = 1:N
        
        % call function. Indexing of rd accommodates empty rd
        [p(i),prm(i,:),ci(i,:,:)] = ...
            fitRD(rd(max(i,isfinite(rd)),:),rv(i,:),dv(i),ciSet,ind,rpdf);
        
%         % if not provided compute all pairwise differences between columns as
%         % ABSOLUTE difference
%         if isempty(rd)
%             foo = bsxfun(@minus,rv(i,:),rv(i,:)');
%             rdTemp = abs(foo(ind));
%         else
%             rdTemp = rd(i,:);
%         end
%         
%         % fit chance distribution and compute p value -- use ABSOLUTE
%         % difference
%         [p(i),prmTemp] = sigTest(abs(dv(i)),rdTemp,'upper',rpdf);
%         % for non-empirical PDFs
%         if ~strcmp(rpdf,'empirical')
%             prm(i,:) = prmTemp;            
%             prmTemp = num2cell(prmTemp); % to allow passing elements as separate arguments
%             ci(i,:,:) = icdf(rpdf,ciSet/100,prmTemp{:});
%         else
%             % store as temp var since prctile can only accept a vector of
%             % CI's
%             ciTemp = prctile(rdTemp,ciSet(:));
%             ci(i,:,:) = reshape(ciTemp,length(ciTemp)/2,2);
%         end
    end
end
clear ind

%% collect output variables

out.vDiff = dv;
out.p = p;
out.nullDistModel = rpdf;
out.param = prm;
out.ciSet = ciSet;
out.ci = ci;

end % main function

%% SUBFunction for computing differences, fitting dist, finding p, etc

function [p,prm,ci] = fitRD(rd,rv,dv,ciSet,ind,rpdf)

% if not provided compute all pairwise differences between columns as
% ABSOLUTE difference
if isempty(rd)
    foo = bsxfun(@minus,rv,rv');
    rd = abs(foo(ind));
end

% fit chance distribution and compute p value
[p,prm] = sigTest(abs(dv),rd,'upper',rpdf);
% for non-empirical PDFs
if ~strcmp(rpdf,'empirical')
    prmTemp = num2cell(prm); % to allow passing elements as separate arguments
    ci = icdf(rpdf,ciSet/100,prmTemp{:});
else
    % store as temp var since prctile can only accept a vector of
    % CI's
    ci = reshape(prctile(rdTemp,ciSet(:)),size(ciSet,1),2);
end

end % subfunction