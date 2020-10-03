function varIxs = varianceIXs(projBasis, varargin)
% varIxs = varianceIXs(projBasis, varargins)
%
% varargins is a cell array and expects
%       varargin{1} -- dataTensor
%       varargin{2} -- codedParams
%       varargin{3} -- (optional) bitPartial -- logical on whether to ALSO compute partial correlations
%       varargin{4} -- (optional) cc -- common condition response --
%                       matrix of Times x Neurons. Needed to compute
%                       variance across time per condition.
% By Gamal Elsayed and Daniel Kimmel, latest version: 14 Mar 2017

%% collect input vars
    dataTensor = varargin{1};
    codedParams = varargin{2};
    if length(varargin) > 2
        bitPartial = varargin{3}; % logical on whether to ALSO compute partial correlations
    else
        bitPartial = false;
    end
    
    if length(varargin) > 3
        cc = varargin{4}; % common condition response
    else
        cc = [];;
    end
    
%%    
    [T, N, C] = size(dataTensor);
    K = size(codedParams, 2);
    if iscell(projBasis)
        projBasis = projBasis{1}; 
    end
    numBasis = size(projBasis, 2);
    bigA = reshape(permute(dataTensor,[1 3 2]), [], N);
    [~, projBasis] = normVects(projBasis);
    randProj = bsxfun(@minus, bigA, mean(bigA))*projBasis;
    
    
    % compute variance across conditions, as function of time. Since the
    % 2nd dimension will always have size 1, permute such that 1st
    % dimension is time and 2nd dimensions are bases.
    var_t = permute(var(reshape(randProj, T, C, numBasis), 0, 2),[1 3 2]);
    
    % compute variance across time, as function of condition
    if ~isempty(cc)
        % add back common condition response
        dataTensorCC = bsxfun(@plus,dataTensor,cc);
        bigACC = reshape(permute(dataTensorCC,[1 3 2]), [], N);
        randProjCC = bsxfun(@minus, bigACC, mean(bigACC))*projBasis;
        % Compute variance across time. Since the 1st dimension will always
        % have size 1, permute such that 1st dimension is conditions and
        % 2nd dimensions are bases.
        var_cond = permute(var(reshape(randProjCC, T, C, numBasis), 0, 1),[2 3 1]);
    end
    
    % initialize
    RSV_t = NaN(T,numBasis,K);
%     aRSV_t = NaN(T,numBasis,K);
    ISV_t = NaN(T,numBasis,K);
    Cr_t = NaN(T,numBasis,K);
    
    % if doing partial, add dimension such that:
    % 1 -- time points
    % 2 -- regression axes (basis vectors
    % 3 -- target signal (K)
    % 4 -- partial correlation controlling for signal in dimension 3 (K-1)
    if bitPartial
        RSV_t_partial = NaN(T,numBasis,K,K-1);
%         aRSV_t_partial = NaN(T,numBasis,K,K-1);
        ISV_t_partial = NaN(T,numBasis,K,K-1);
        signalN = NaN(K,K-1);
        % difference in sqaured correlation coefficients r2_bp - r2_ab,
        % where r2_bp is correlation between off-target signal b (dimension
        % 4) and r2_ab is correlation between off-target signal b and
        % on-target signal a (dimension 3).
        diffCorr_proj_signal = NaN(T,numBasis,K,K-1); 
    end
    
    % check that number of conditions in data and number of conditions in
    % codedParams match. Also check that there is more than 1 condition, as
    % we cannot compute variance across conditions with only 1 condition.
    if C > 1 && C == size(codedParams,1)
        % this is the main code 
        for n = 1:numBasis
            for k = 1:K
                [RSV_t(:, n, k), ISV_t(:, n, k), Cr_t(:, n, k)] = condSeparation(randProj(:, n), codedParams(:, k));
%                 [RSV_t(:, n, k), aRSV_t(:, n, k), ISV_t(:, n, k), Cr_t(:, n, k)] = condSeparation(randProj(:, n), codedParams(:, k));
                if bitPartial
                    % set of signal index NOT including k
                    pos = setdiff([1:K],k);
                    for k2 = 1:K-1
                        % RSV based on signal given in pos(k2) controlling
                        % for signal k1. Also extracting correlation r_bp
                        % between signal pos(k2) and data, regardless of
                        % signal k1 (i.e., this is not the partial
                        % correlation coefficient). Also extracting the
                        % correlation r_ab between signal pos(k2) and
                        % signal k1.
                        [RSV_t_partial(:, n, k, k2), ...
                            ISV_t_partial(:, n, k, k2),...
                            ~,~,r_bp,r_ab] = ...
                            condSeparation(randProj(:, n), codedParams(:, pos(k2)), ...
                            codedParams(:, k));
%                         [RSV_t_partial(:, n, k, k2), ...
%                             aRSV_t_partial(:, n, k, k2), ...
%                             ISV_t_partial(:, n, k, k2),...
%                             ~,~,r_bp,r_ab] = ...
%                             condSeparation(randProj(:, n), codedParams(:, pos(k2)), ...
%                             codedParams(:, k));
                        
                        % compute difference of correlation coefficients:
                        diffCorr_proj_signal(:, n, k, k2) = r_bp.^2 - r_ab.^2;
                    end
                    % store which positions are referenced in Dim 4
                    signalN(k,:) = pos;
                end
            end
        end
        
    elseif C ~= 1 
        % this is a special check
        error('Number of conditions in data do not match number of conditions in coded parameters or only one condition provided')
    end
    
    % store
    varIxs.var_t = var_t;
    if ~isempty(cc)
        varIxs.var_cond = var_cond;
    end
    varIxs.RSV_t = RSV_t;
%     varIxs.aRSV_t = aRSV_t;
    varIxs.ISV_t = ISV_t;
    varIxs.Cr_t = Cr_t;
    if bitPartial
        varIxs.RSV_t_partial = RSV_t_partial;
%         varIxs.aRSV_t_partial = aRSV_t_partial;
        varIxs.ISV_t_partial = ISV_t_partial;
        varIxs.signalN = signalN; 
        varIxs.diffCorr_proj_signal = diffCorr_proj_signal;
    end
end