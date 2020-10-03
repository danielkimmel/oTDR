function  c = colorCategorical(nColor,varargin)
% c = colorCategorical(nColor)
%
% Returns m x 3 matrix C with RGB color specifications for up to 4
% categories, as specified in nColor.
%
% Optionally accepts PROJECT, string specifying specific project for which
% specific colors have been assigned. Current options: 'CB'
%
% Based on colorbrewer2.org
%
% Daniel Kimmel, 11 April 2013

if length(varargin) > 0 
    project = varargin{1};
else
    project = [];
end

if ~isempty(project)
    switch project
        case 'CB'
            c = [
                24 173 150
                213 124 50
                123 136 192
            ]/255;
        otherwise
            error('Project %s was not recognized',project)
    end
elseif nColor > 4
%     warning('Colors for > 4 categories not supported. You requested %d. Using HSV instead',nColor);
    
    c = hsv(nColor);
else

    c = [
        51 160 44
        178 223 138
        31 120 180
        145 179 193
        ]/255;
    % 166 206 227
    c = c(1:nColor,:);

end
