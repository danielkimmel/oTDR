function [h,hCB] = plotScatter3(x,y,z,varargin);
%
% [h,hCB] = plotScatter3(x,y,z,[name,value pairs]) 
% 
% Plots 3D scatter of X, Y and Z in current axes. X, Y, and Z are nx1
% vectors.  Optionally colors data points according to some arbitrary value
% passed in 'color' (see below). 
%
% Returns H, a vector of handles to each data point.  Also returns hCB, a
% handle to the optional colorbar.
%
% The following can be passed as name, value pairs:
%
% color -- nx1 vector of raw data values that will be linearly mapped to a
% color scale (see dataRange and colormapName). Instead, can be n x 3 RGB
% matrix specifying exact color of each point.  Lastly, can instead by a
% single character string specifying a single color for all points.
% Defaults to [] -- no color code.
%
% markerSize -- scalar for size of all points. Instead, can be nx1 vector
% specifying marker size for each point. Note, if vector is passed, it will
% supercede highlightSize.  Note that sizes should be 10+ times normal
% Matlab sizes. Default = 40.
%
% bitSig -- nx1 logical vector on whether each point is signficant and
% should be highlighted according to highlight parameters below. Defaults
% to all "false" (not signficant).
%
% marker -- single character string specifying marker type for all points.
% Can be nx1 character array specifying marker of each point, in which case
% it overrides highlightMarker.  Defaults to [], which saves looping
% through each point and plots the default 'o' marker.
% 
% highlightEdge -- color of highlight edge on "significant" points (see
% bitSig). Can be single character string or 1x3 RGB vector.  Set to [] to
% not highlight edge.  Default to [].
% 
% highlightMarker -- single character array specifying marker type of
% "significant" points (see bitSig).  Defaults to [] (no different than
% other points).
%
% highlightSize -- scalar specifying size of significant markers (see
% bitSig).  If different sizes are needed for each point, specify in
% markerSize.  Defaults to [] (no different than specified in markerSize).
% 
% bitFillMarker -- logical on whether to fill markers. Defaults to true.
%
% dataRange -- 1 x 2 vector with the minimum and maximum values to use as
% the color scale in the same units as 'color'.  When 'color' contains RGB
% values, then dataRange only specifyies the labeling of the colorbar.
% Defaults to the actual range passed in color (when color contains data
% values). 
%
% nullH -- scalar value of the null hypothesis which will be labeled on the
% colorbar.  Defaults to [].
% 
% nColor -- number of discrete colors to use in colormap. Defaults to 256.
%
% colormapName -- string with name of MATLAB colormap to use.  Defaults to
% 'jet'.
%
% bitColorbar -- logical on whether to plot a colorbar.  Defaults to true.
%
% datatipLabel -- nx1 cell array of custom text to attach to the datatip of
% each point.  Defaults to {}.
%
% onClickFunction -- contains a single function handle to a callback
% function to be called when the user mouse-clicks on a point in the
% scatter plot.  If the callback function requires data to passed to it, 
% place that data in onClickData. In this case, a cell array will be passed
% to the callback function with the first cell containing the function
% handle and each subsequent cell containing 1 variable to be passed to the
% callback function.  
%
% onClickData -- nx1 cell array of data to be passed to onClickFunction.
% Each cell should contain a cell array with 1 variable per cell to be
% passed to onClickFunction.
%
% VERSION NOTE: the function relies on MATLAB function SCATTER3 and calls
% the "v6" option, which uses an antiquated version of SCATTER3 so as to
% return separate handles to each data point (the current version of
% SCATTER3 does not support this for >100 data points).  When MATLAB stops
% supporting the "v6" option, the present function will have to be modified
% to use PLOT3 to plot each point separately in a loop.
%
% Daniel Kimmel, 18 Oct 2010


%% default values for optional params
marker = []; % will default to 'o'.  can be column of markers.
bitSig = false(size(x));
color = []; % data values, n x 3 RGB matrix, or single char string of color value
markerSize = 40; % can also be vector, supercedes highlight size
highlightEdge = []; % color as string or RGB. leave empty to not highlight edge.
highlightMarker = []; % scalar only. leave empty to match normal marker
highlightSize = []; % can only be scalar or empty.
bitFillMarker = true; % logical on whether to fill markers.

nullH = [];
dataRange = []; % defaults to range passed in color. 2 element row vector. 
%               specifies color range and colorbar label range when color
%               contains data. specifies only colorbar label when color
%               contains RGB values  

datatipLabel = {};
onClickFunction = [];
onClickData = {};

nColor = 256; % number of colors to use
colormapName = 'jet';
bitColorbar = true;

%% collect optional name,value pairs
warnopts(assignopts(who, varargin));

%% construct appearance vars

% if size vars are present and scalar, make them into vectors according to
% bitSig. Only do this if at least one point has bitSig==1.  The size var
% is a special case since we need to pass a vector to scatter3(), whereas
% the other features have to be set after the scatter3() call.
if length(markerSize) == 1 && length(highlightSize) == 1 && any(bitSig)
    markerSize = repmat(markerSize,(size(x)));
    
    markerSize(bitSig) = highlightSize;
elseif length(markerSize) > 1 && length(highlightSize) == 1 
    warning('Custom highlight size IGNORED because markerSize input specifies size for ALL points. Use scalar for markerSize if highlightSize is desired.')
    highlightSize = [];
end

% if marker highlight var is present and scalar, but marker was already
% specified as a vector, warn that highlight will be ignored.
if length(marker) > 1 && length(highlightMarker) == 1 
    warning('Custom highlight marker IGNORED because MARKER input specifies marker for ALL points. Use scalar for MARKER if highlightMarker is desired.')
    highlightMarker = [];
end


%% remove nans

% build index for data points PRIOR to removing NaNs
obsInd = 1:length(x);

bitRemove = isnan(x) | isnan(y) | isnan(z);

% also test if optional vars have nans
if size(color,1) > 1
    bitRemove = bitRemove | any(isnan(color),2);
end
% if length(markerSize) > 1
%     bitRemove = bitRemove | isnan(markerSize);
% end
% if length(marker) > 1
%     bitRemove = bitRemove | isnan(marker);
% end

% now remove entries that have nans in any vary
if size(color,1) > 1
    color(bitRemove,:) = [];
end
if length(markerSize) > 1
    markerSize(bitRemove) = [];
end
if length(marker) > 1
    marker(bitRemove) = [];
end
if length(datatipLabel) > 1
    datatipLabel(bitRemove) = [];
end
if length(onClickData) > 1
    onClickData(bitRemove) = [];
end
x(bitRemove) = [];
y(bitRemove) = [];
z(bitRemove) = [];
bitSig(bitRemove) = [];
obsInd(bitRemove) = [];


%% plotting

% handle for attaching custom datatips:
hdt = datacursormode;
set(hdt,'UpdateFcn',{@scatterDatatipCallback});

% instatiate:
bitObsPlotted = false(size(x));
obsind = cell(1,size(x,2));

% define colormap for figures
eval(['colormap(',colormapName,'(nColor))']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build string for plotting scatter
str = 'h = scatter3(''v6'',x,y,z';
% add term for marker size
if ~isempty(markerSize)
    str = [str,',markerSize'];
end
% add term for data/color
if ~isempty(color)
    str = [str,',color'];
end
% add term for filling-in data points
if bitFillMarker
    str = [str,',''filled'''];
end
% close string
str = [str,');'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT scatter
% note that "v6" functionality will no longer be supported in future
% releases and we'll have to loop through each point plotting it
% separately.
warning('off',['MATLAB:scatter3:DeprecatedV6Argument'])
eval(str);
clear str
warning('on',['MATLAB:scatter3:DeprecatedV6Argument'])

% change marker if requested:
if ~isempty(marker)
    if length(marker) ==  1
        % set marker en masse
        set(h,'Marker',marker);
    else
        % set marker point by point
        for i = 1:length(h)
            set(h(i),'Marker',marker(i));
        end
    end
end

% mark significant clusters en masse
if any(bitSig) && ~isempty(highlightEdge)
    set(h(bitSig),'MarkerEdgeColor',highlightEdge);
end
if any(bitSig) && ~isempty(highlightMarker)
    set(h(bitSig),'Marker',highlightMarker);
end

% bring significant clusters to front
[foo,i] = sort(bitSig,'descend');
kids = get(gca,'Children');
set(gca,'Children',kids(i));

grid off

% Attach custom datatips to set of points
% requires looping through all handles
for i = 1:length(x)
    setappdata(h(i),'obsind',obsInd(i));
    if ~isempty(datatipLabel)
        setappdata(h(i),'datatipLabel',datatipLabel{i});
    end
end

% Assign callback function and attach callback data
if ~isempty(onClickFunction)
    if isempty(onClickData)
        for i = 1:length(x)
            % if data is included to be passed to the callback function...
            set(h(i),'ButtonDownFcn',onClickFunction);
        end
    else
        % if data is included to be passed to the callback function...
        for i = 1:length(x)
            set(h(i),'ButtonDownFcn',{onClickFunction,onClickData{i}});
        end
    end
end

% scale color bar if specific data range is supplied 
if ~isempty(dataRange)
    caxis(dataRange);
end

% draw colorbar only if color data provided
if ~isempty(color) && bitColorbar
    hCB = colorbar;
    % add tick for nullH in colorbar
    if ~isempty(nullH)
        yTick = get(hCB,'YTick');
        if ~ismember(nullH,yTick)
            yTick = sort([yTick nullH]);
            set(hCB,'YTick',yTick);
        end
    end
else
    hCB = [];
end

% turn off data cursor mode
datacursormode;

%% callback function for custom text on datatips
function datatipTxt = scatterDatatipCallback(obj,evt)

target = get(evt,'Target');
i = get(evt,'DataIndex');
pos = get(evt,'Position');

% retrieve attached data
% not all elements in plot necessarily have data attached
if isappdata(target,'obsind')
    obsind = getappdata(target,'obsind');
else
    obsind = [];
end
if isappdata(target,'datatipLabel')
    datatipLabel = getappdata(target,'datatipLabel');
else
    datatipLabel = {};
end

datatipTxt = {...
    ['x: ' num2str(pos(1))]...
    ['y: ' num2str(pos(2))]...
    ['z: ' num2str(pos(3))]
    };

% add data value coded by color
if ~isempty(get(target,'CData'))
    datatipTxt{end+1} = ['cdata: ' num2str(get(target,'CData'))];
end

% add observation number
if ~isempty(obsind)
    datatipTxt{end+1} = ['Obs: ' num2str(obsind(i))];
end

if ~isempty(datatipLabel) 
    if iscell(datatipLabel)
        datatipTxt{end+1} = datatipLabel{i};
    else
        datatipTxt{end+1} = datatipLabel;
    end
end
