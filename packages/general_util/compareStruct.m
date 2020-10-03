function [bitEq, fnUneq] = compareStruct(s1,s2,varargin)
%
% [bitEq, fnUneq] = compareStruct(s1,s2,[bitRecursive],[tol],[skipField],[bitVerbose])
% 
% Compares fields in struct s1 to fields with the same name in struct s2.  
% Returns bitEq = 0 when 
% - Matching fields are not equal.
% - Fields in S1 are not present in S2
% Otherwise returns bitEq = 1.
% When bitEq = 0, fnUneq contains field name from S1 that fulfilled one of
% the conditions above.
% Fields in S2 that are not present in S1 are ignored.
% Optionally pass 
% - BITRECURSIVE (Boolean) to force recursion through structs
%       within S1. Default = 1.
% - TOL (scalar) to tolerate numerical differences smaller than TOL.
%       Defaults to 10^-8
% - SKIPFIELD (cell array) -- skips fields whose name exactly matches
%       element of SKIPFIELD. Note that in the case of recursive calls by
%       compareStruct(), skipField will apply to matching fields at all
%       levels of the structure hierarchy.
% - BITVERBOSE (Boolean) to display messages in command window. Default =
%       TRUE
% 
% Daniel Kimmel, 16 Oct 2016

%% collect optional vars

if length(varargin) > 0
    bitRecursive = varargin{1};
else
    bitRecursive = true;
end
if isempty(bitRecursive)
    bitRecursive = true;
end

if length(varargin) > 1
    tol = varargin{2};
else
    tol = 10^-8;
end
if isempty(tol)
    tol = 10^-8;
end

if length(varargin) > 2
    skipField = varargin{3};
else
    skipField = {};
end

if length(varargin) > 3
    bitVerbose = varargin{4};
else
    bitVerbose = true;
end
if isempty(bitVerbose)
    bitVerbose = true;
end

% if length(varargin) > 2
%     ignoreField = varargin{3};
% else
%     ignoreField = {};
% end


%% 

% initialize
bitEq = true;
fnUneq = [];

% get field names
fn1 = fieldnames(s1);
% fn2 = fieldnames(s2);

% in case s1 is not a scalar
if length(s1) ~= length(s2)
    bitEq = false;
    fnUneq = '[top level]';
    return
end

for j = 1:length(s1)
    % loop through s1 fields
    for i = 1:length(fn1)
        % skip field if requested
        if ismember(fn1(i),skipField)
            if bitVerbose
                fprintf('Skipping field %s\n',fn1{i})
            end
            continue
        end
        % return if field is not present in s2
        if ~isfield(s2(j),fn1{i})
            bitEq = false;
            fnUneq = fn1{i};
            return
        end

        % if field is a substruct...
        if isstruct(s1(j).(fn1{i}))
            % ...and not a struct is s2, then return
            if ~isstruct(s2(j).(fn1{i}))
                bitEq = false;
                fnUneq = fn1{i};
                return
            end
            % ...and recursion is ON, then process it
            if bitRecursive
                [bitEq, fnUneqRecur] = compareStruct(s1(j).(fn1{i}),s2(j).(fn1{i}),bitRecursive,tol,skipField,bitVerbose);
                % after a recursive call that is unequal, propagate up the
                % failure
                if ~bitEq
                    fnUneq = [fn1{i},'.',fnUneqRecur];
                    return
                end
            else
                % if recusion is OFF, then we effectively skip the field
                warning('Will not compare field %s because it is a substructure and bitRecursive = FALSE',fn1{i})
            end

        else
            % if field is not a struct, compare directly
            if ~isequaln(s1(j).(fn1{i}),s2(j).(fn1{i}))
                % attempt to compare numerically
                if isnumeric(s1(j).(fn1{i})) 
                    foo = abs(s1(j).(fn1{i}) - s2(j).(fn1{i}));
                    if ~any(foo(:) > tol) % code as "not any" so as to ignore NaNs
                        % ignore the difference
                        continue
                    end
                end
                bitEq = false;
                fnUneq = fn1{i};
                return
            end
        end
    end
end
    