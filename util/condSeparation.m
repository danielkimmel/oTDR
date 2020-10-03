function [condSep_t, condSepBar_t, Cr, varm, Cr_nonpartial, Cr_Corig] = condSeparation(bigA_projC, Corig, Corig_control)
% [condSep_t, condSepBar_t, Cr, varm, Cr_nonpartial, Cr_Corig] = condSeparation(bigA_projC, Corig, [Corig_control])
%
% Corig_control -- when included, returns values based on partial
% correlation between projection and Corig controlling for variables in 
% Corig_control, where each column is a different variable and each row is
% an observation corresponding to those in Corig
%
% RETURNS
% condSep_t -- relevant signal variance (RSV), or the portion of variance
%       explained (varm) that is POSITIVELY related to the signal passed in
%       Corig.  
%           RSV = varm * max(0,corr(bigA_projC, Corig))^2
% [DEFUNCT] condSepAnti_t -- anti-relevant signal variance (aRSV), the portion of
%       variance NEGATIVELY related to Corig (Cr < 0).  
%           aRSV = varm * min(0,corr(bigA_projC, Corig))^2
% condSepBar_t -- irrelevant signal variance (ISV), or the portion of
%       variance explained that is not related to the signal passed in
%       Corig.
%           ISV = varm * (1 - RSV - aRSV)
% Cr -- correlation between data bigA_projC and signal Corig, optionally
%       controlling via partial correlation for correlation between data
%       and signal Corig_control (when Corig_control is passed as input)
% varm -- variance in rows of bigA_projC
% Cr_nonpartial -- correlation between data bigA_projC and signal Corig
%       regardless of whether controlling for signal Corig_control. When 
%       Corig_control is not provided, Cr_nonpartial = Cr.
% Cr_Corig -- correlation between signal Corig and control signal
%       Corig_control. Returns NaN when control signal is not provided.
%
% Gamal Elsayed and Daniel Kimmel. Updated 2017 Apr 21 to divide RSV into
% RSV and aRSV.

if nargin <= 2
    Corig_control = [];
end

numConds = length(Corig);
trialL = size(bigA_projC,1)./numConds;
C = bsxfun(@minus, Corig, mean(Corig));
C = C./norm(C);
A_projC = reshape(bigA_projC, trialL, numConds);
A_projC = bsxfun(@minus, A_projC, mean(A_projC,2)); % subtract common condition response


% % condSep_t = ((A_projC*C)/sqrt(length(C))).^2;
% % [~, nullC] = normVects(null(C'));
% % condSepBar_t = 0;
% % for i = 1:size(nullC,2)
% % condSepBar_t = condSepBar_t+((A_projC*nullC(:,i))/sqrt(length(C))).^2;
% % end

% standard correlation
Cr_nonpartial = corr(A_projC',C);

if ~isempty(Corig_control)
    C_control = bsxfun(@minus, Corig_control, mean(Corig_control));
    C_control = bsxfun(@rdivide,C_control,normVects(C_control)');
    
    % compute correlation between the main predictor Corig and the
    % control predictor Corig_control
    Cr_Corig = corr(C,C_control);

    % partial correlation
%     Cr = partialcorr(A_projC',C,C_control);

    % semipartial (or part) correlation:
    % http://faculty.cas.usf.edu/mbrannick/regression/Partial.html 
    % takes form Cr = (rbp - rap*rab) / sqrt((1-rab^2)), were rbp is
    % correlation between projection p and main predictor b (Corig), rap is
    % between p and control predictor a (Corig_control), and rab is between
    % a and b.
    Cr = (Cr_nonpartial - corr(A_projC',C_control)*Cr_Corig) / sqrt((1-Cr_Corig^2));
    
else
    % when not doing partial correlation, Cr = Cr_nonpartial
    Cr = Cr_nonpartial;
    
    % not computed:
    Cr_Corig = NaN;
end

% sigm = var(A_projC, 1,2).^0.5; % GAMAL had var(~,1,~), which finds the second moment of the mean
% sigm = var(A_projC, 0,2).^0.5; % DK things it should be var(~,0,~), which finds the traditional variance about the mean 
% condSep_t = (sigm.*Cr).^2;
% condSepBar_t = (sigm.^2.*(1-Cr.^2));

% regardless of above changes regarding var(~,0,~) vs var(~,1,~), makes
% sense to compute variance instead of std dev:
varm = var(A_projC, 0,2); % DK things it should be var(~,0,~), which finds the traditional variance about the mean 

%%%% COMPUTE RELEVANT AND IRRELEVANT SIGNAL VARIANCE (RSV, ISV)
% Here we changed the definition of RSV (condSep_t) to include only the porition of
% variance explained POSITIVELY by predictor Corig, that is, when Cr >= 0. 
% We defined a new term, anti-relevant signal variance (aRSV, or condSepAnti_t) to capture
% the portion of variance NEGATIVELY related to Corig (Cr < 0). 
% The remaining variance was defined as the irrelevant signal variance
% (ISV, or condSepBar_t), which was unchanged from before: ISV = 1 - RSV - aRSV
% irrelevant variance.
condSep_t = varm .* Cr.^2; 
condSepBar_t = varm - condSep_t; % OLD: varm.*(1-Cr.^2);
% condSep_t = varm .* max(Cr,0).^2; % OLD: varm .* Cr.^2
% condSepAnti_t = varm .* min(Cr,0).^2;


end