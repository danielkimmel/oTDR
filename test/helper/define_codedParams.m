function [codedParams, bitSingleton, paramLabel] = define_codedParams(Data, bitSymmetricChoiceEncoding)
% [codedParams, bitSingleton] = define_codedParams(Data, bitSymmetricChoiceEncoding)
%
% Helper function for oTDR to define the targeted task variables.
% 
% ACCEPTS:
% Data -- 1 x C struct (returned by genDataStruct()) defining data and task
%       variables for each of C conditions.
% bitSymmetricChoiceEncoding -- logical on whether "choice" is encoded symmetrically as [-1 1] (TRUE) or
%       assymetrically as [0, 1] (FALSE). This determines whether expected
%       reward (benefit * choice) varies with benefit for both accept and
%       reject choices, or only for accept choices, respectively.
% 
% RETURNS:
% codedParams -- C x P matrix where each row is a unique condition and each
%       column is a task variable (i.e. regression predictor). NOTE:
%       Predictor scaling to [0, 1] occurs separately within oTDR.
% bitSingleton -- C x 1 logical vector as to whether a given condition
%       was a "singleton" condition, i.e., whether 8 drops of juice was
%       indicated with a single purple icon instead of 8 yellow icons.
% paramLabel -- 1 x P cell array of labels for the task parameters
% (variables)

% extract task variables (i.e., parameters) for each condition
condition_params = vertcat(Data.Predictors);

% define codedParams matrix, where each row is a unique condition and each
% column is a task variable and putative target for oTDR. (NOTE: Parameter
% scaling to [0, 1] occurs within oTDR.)
%
% benefit (number of rewards offered)
codedParams(:,1) = condition_params(:,1);

% choice (accept or reject)
% determine symmetric vs asymmetric reward encoding (see above)
if bitSymmetricChoiceEncoding
    % symmetric encoding: choice is -1 or 1 (will be scaled later to -0.5
    % to 0.5)
    codedParams(:,2) = (condition_params(:,2)==0).*-1+(condition_params(:,2)~=0).*1; 
else
    % asymmetric encoding: choice is 0 or 1 
    codedParams(:,2) = (condition_params(:,2)==0).*0+(condition_params(:,2)~=0).*1;
end

% expected reward (benefit * choice)
codedParams(:,3) = codedParams(:,1).*codedParams(:,2); 

% determine whether condition was a singleton condition
bitSingleton = condition_params(:,3);

paramLabel = {'Benefit','Choice','Expected Reward'};