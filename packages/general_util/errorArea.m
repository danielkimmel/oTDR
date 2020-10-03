function [hL, hF] = errorArea(x,y,errL,color,varargin)
% [hL, hF] = errorArea(x,y,errL,color,...)
%
% Plot line with surrounding filled area representing error about the line.
% Provide vectors of equal size: X for horizontal (time) points, Y for
% vertical points, and ERRL for +/- error from line. Plots line and filled
% area in COLOR specified in LineSpec or RGB vector format.  To use default
% axes color, leave COLOR as empty [].
%
% Note that X cannot contain NaN values. However, where ERRL contains NaN
% values, no shading will appear. (This occurs when ERRL represents the
% s.e.m. and the number of trials at a particular timepoint is zero). Where
% Y contains NaN values, no line and no shading will appear.
%
% Returns hL and hF. hL is a scalar with handle for the line. hF is a
% vector with handles to the fill areas -- multiple fill areas are
% generated when finite values of Y and ERRL are interspersed with NaN
% values.
%
% Optionally provide the following as name,value pairs:
% transparency = scalar from 0 (opaque) to 1 (invisible) specifying degree
%    of transparency. Defaults to 0.2.
% lineSpec = string specifying color, marker, and line style. Color is moot
%    since it is overridden by separate COLOR value.
% errU = vector of same size as X, Y, and ERRL that specifies the upper
%    error from line. When not provided, the upper error is assumed to be
%    the same as the lower error specified by ERRL. NOTE that error should
%    always be a positive value. The sign of error is determined by whether
%    it is lower error (less than Y) or upper error (greater than Y).
% bit4Illustrator = logical on whether to eliminate fill color from patch
%    objects and instead plot the outline of the patch area. This allows
%    export to Adobe Illustrator, which cannot handle filled patches.
%
% NOTE: plotting a transparent patch on a Mac causes the axes to disappear.
% The workaround is to increase the line thickness of the axes to >=1.5,
% which this function now automatically implements. The only obvious
% downside is that the tick marks get thicker, too.
%
% Depends on custom functions: warnopts(), assignopts()
%
% Daniel Kimmel, 14 Feb 2009. Updated 30 Apr 2017 to remove shading when
% exporting to Illustrator.

%% default values for optional vars
transparency = 0.2;
lineSpec = [''];
errU = [];
bit4Illustrator = false;

%% collect optionally provided parameters
warnopts(assignopts(who, varargin));

%% check input vars
if isempty(transparency)
    transparency = 0.2;
end

% made x, y, and errL into row vectors
if size(x,1) > size(x,2)
    x = x';
end
if size(y,1) > size(y,2)
    y = y';
end
if size(errL,1) > size(errL,2)
    errL = errL';
end
if ~isempty(errU) & size(errU,1) > size(errU,2)
    errU = errU';
end

% if color is empty, use default color for axes
if isempty(color)
    foo = get(gca,'ColorOrder');
    color = foo(1,:);
end

%% check

if any(errU < 0) || any(errL < 0)
    error('Error provided in errL (and optionally errU) must be positive sign')
end

%% save current plot values
nextPlot = get(gca,'NextPlot');

%% plot fill bounds

% upper and lower boarder of polygon
yL = y-errL;
if ~isempty(errU)
    yU = y+errU;
else
    yU = y+errL;
end
%yU = y(end:-1:1) + errL(end:-1:1);

% find good (non-NaN and finite) sections
bitGood = isfinite(yL);
% find start and end of good sections:
if bitGood(1) == 1
    startGood = 1;
else
    startGood = [];
end
startGood = [startGood, find(diff(bitGood) == 1) + 1];
endGood = find(diff(bitGood) == -1);
% when the vector ends with a good section
if length(endGood) < length(startGood)
    endGood(end+1) = length(bitGood);
end
if length(startGood) ~= length(endGood)
    error('Every "good" fill section must have both a start and end position');
end

% need at least one good section
if isempty(startGood)
    % quit without plotting
    warning('No finite data points passed to errorArea.  Will quit without plotting')
    [hL hF] = deal([]);
    return
end

% loop through each good section, and make fill separately:
for i = 1:length(startGood)
    % forward direction
    xTemp = x(startGood(i):endGood(i));
    yTempFill = yL(startGood(i):endGood(i));
    % reverse direction
    xTemp = [xTemp, [x(endGood(i):-1:startGood(i))]];
    yTempFill = [yTempFill, [yU(endGood(i):-1:startGood(i))]];
    
    % plot fill
    hF(i) = fill(xTemp, yTempFill,color);
    if bit4Illustrator
        set(hF,'EdgeColor',color);
        set(hF,'FaceColor','none');
    else
        set(hF,'EdgeColor','none');
        set(hF,'FaceAlpha',transparency);
    end

    hold on;
end

%% plot center line
hL = plot(x,y,lineSpec);
% set(hL,'LineWidth',2);
set(hL,'Color',color);

%% fix axes appearance bug
 
% plotting a transparent patch on a Mac causes the axes to disappear.
% The workaround is to increase the line thickness of the axes to >= 1.5,
% which this function now automatically implements.

if ismac && transparency < 1 && get(gca,'lineWidth') < 1.5
    set(gca,'lineWidth',1.5);
end
   
%% return axes to original state:
set(gca,'NextPlot',nextPlot);

%% gather handles
% handle = [hL hF];

