function  c = colorContinuous(nColor,varargin)
% c = colorContinuous(nColor)
%
% Returns m x 3 matrix C with RGB color specifications for NCOLOR
% continuous values. When NCOLOR > 10, built in MATLAB colormaps are used.
% When NCOLOR <= 10, custom colorblind-safe colormaps are used based on
% colorbrewer2.org
%
% Optionally accepts PROJECT, string specifying specific project for which
% specific colors have been assigned. Currently accepts:
%   'CT' -- context task
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
        case 'CT'
            % no custom set defined
%             c = [
%             ]/255;
        otherwise
            error('Project %s was not recognized',project)
    end
elseif nColor > 10
%     warning('Colors for > 4 categories not supported. You requested %d. Using HSV instead',nColor);
    
    c = hsv(nColor);
else
    % define full color map
    
    % from http://colorbrewer2.org/?type=diverging&scheme=RdBu&n=10
    c = [
        103,0,31
        178,24,43
        214,96,77
        244,165,130
        253,219,199
        209,229,240
        146,197,222
        67,147,195
        33,102,172
        5,48,97
    ]/255;

%     % from http://colorbrewer2.org/?type=diverging&scheme=RdYlBu&n=10
%     c = [
%         165,0,38
%         215,48,39
%         244,109,67
%         253,174,97
%         254,224,144
%         224,243,248
%         171,217,233
%         116,173,209
%         69,117,180
%         49,54,149        
%     ]/255;

    % flip so that lower values are bluer
    c = flipud(c);

    % pare down colormap to nColor, excluding middle values
    nExcl = size(c,1) - nColor;
    c(1+floor(nColor/2):floor(nColor/2)+nExcl,:) = [];

end
