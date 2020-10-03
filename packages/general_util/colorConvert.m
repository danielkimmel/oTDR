function c2 = colorConvert(c1,bitLong)
%
% c2 = colorConvert(c1,[bitLong])
% 
% Converts color specification C1 to altnernative and equivalent color
% specification C2. 
%
% Behahvior depends on input C1. When C1 is a character, C2 is either the
% alternative character spec (i.e., when C1 is short name, C2 is long name)
% or if bitLong==true, then C2 is the RGB triplet. When C2 is an RGB
% triplet, then C2 is either the equivalent short name, or if
% bitLong==true, the equivalent long name. bitLong defaults to false.
%
% E.g.,
%   c2 = colorConvert('r') --> returns 'red'
%   c2 = colorConvert('r',true) -- returns [1 0 0] 
%   c2 = colorConvert('red') -- returns 'r' 
%   c2 = colorConvert([1 0 0],1) -- returns 'red'
%
% Daniel Kimmel, 4 Mar 2017

%% handle input vars
if nargin <= 1
    bitLong = false;
end

%% define lookup table

lut = {...
    [1 1 0], 'y', 'yellow'
    [1 0 1], 'm', 'magenta'
    [0 1 1], 'c', 'cyan'
    [1 0 0], 'r', 'red'
    [0 1 0], 'g', 'green'
    [0 0 1], 'b', 'blue'
    [1 1 1], 'w', 'white'
    [0 0 0], 'k', 'black'
    };

%% process c1 according to type

% define lookup columns
if ischar(c1)
    % process c1 according to length
    if length(c1) == 1
        % set c1 column to short name
        col1 = 2;
        
        if bitLong
            % set c2 column to RGB
            col2 = 1;
        else
            % set c2 column to long name
            col2 = 3;
        end
    else
        % set c1 column to long name
        col1 = 3;

        if bitLong
            % set c2 column to RGB
            col2 = 1;
        else
            % set c2 column to short name
            col2 = 2;
        end
        
    end
elseif isnumeric(c1)
    % c1 is RGB
    col1 = 1;
    
    if bitLong
        % set c2 column to long name
        col2 = 3;
    else
        % set c2 column to short name
        col2 = 2;
    end
        
else
    error('Input color spec c1 must be either character array or RGB vector')
end
    
% look up row based on c1, and return entry for c2.
if col1==1
    % for RGB input, must convert LUT to matrix and use rows
    c2 = lut{ismember(cell2mat(lut(:,col1)),c1,'rows'),col2};
else
    c2 = lut{ismember(lut(:,col1),c1),col2};
end

if isempty(c2)
    error('Input color spec c1 not recognized')
end