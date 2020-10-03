function [dataTensor, Summary] = preprocess4TDR(dataTensor, params, varargin)
% [dataTensor, Summary] = preprocess4TDR(dataTensor, params, [bitSubtractCommonCond])
% 
% processes the neural response data dataTensor by transforming it to
% z-scores and optionally subtracting the common condition response.
% 
% Accepts
% dataTensor -- neural response as times x neurons x conditions
% params -- optional structure used to z-transform and then remove common
% condition from dataTenor according to fields:
%       .mu -- mean
%       .nrm -- std deviation
%       .commonConds -- matrix of mean response across conditions (neurons *
%       conditions)
%       (Note that all fields are optional. If not provided, they will be
%       computed according to dataTensor. These fields are what is returned
%       by this same function).
% bitSubtractCommonCond -- logical on whether to subtract the common
%       condition response. Default is true.
% 
% Returns
% dataTensor -- neural response transformed (z-transform and common
%       condition optionally removed) and returned as times * neurons *
%       conditions.
% Summary -- structure with same fields as input PARAMS. Values reflect
%       input values when inputs were provided. Otherwise values are
%       computed by function and returned in their respective fields.
% 
% Gamal Elsayed and Daniel Kimmel, 2017 Feb 01

    %% collect optional vars
    
    if length(varargin) > 0
        bitSubtractCommonCond = varargin{1};
        if isempty(bitSubtractCommonCond)
            bitSubtractCommonCond = true;
        end
    else
        bitSubtractCommonCond = true;
    end

    [T, N, C] = size(dataTensor);
    
    %% zcsore
    XN = reshape(permute(dataTensor,[1 3 2]), [], N);
    nrm = 1./(std(XN));
    mu = mean(XN);
%     warning('DK thinks this next commented line should occur after extracting optional input params mu and nrm')
%     XN = bsxfun(@times, bsxfun(@minus, XN, mu), nrm);
    if exist('params', 'var')
        if isfield(params, 'mu')
            mu = params.mu;
        end
        if isfield(params, 'nrm')
            nrm = params.nrm;
        end
    end
    XN = bsxfun(@times, bsxfun(@minus, XN, mu), nrm);
    dataTensor = permute(reshape(XN, T, C, N), [1 3 2]);
    
    %% compute cross-condition and temporal variance on normalized data
    Summary.var_acrossCond = var(dataTensor,[],3);
    Summary.var_acrossTime = squeeze(var(dataTensor,[],1));
    
    %% common condition subtraction
    commonConds = mean(dataTensor, 3);
    if exist('params', 'var')
        if isfield(params, 'commonConds')
            commonConds = params.commonConds;
        end
    end
    if bitSubtractCommonCond
        dataTensor = dataTensor-repmat(commonConds, 1, 1, C);
    end
    %% Summary
    Summary.mu = mu;
    Summary.nrm = nrm;
    Summary.commonConds = commonConds;
end



