function [figTitleHand] = annotateFig(figName,mfile,figH)
%
% [figTitleHand] = annotateFig([figName],[mfile],[figH])
%
% Annotates the current figure with the string passed in FIGNAME, the
% current date/time, and the m-filename passed in MFILE. Optionally accepts
% handle figH to figure. Returns handle FIGTITLEHAND to the annotation.
%
% Daniel Kimmel, 2017 Mar 05

% collect input vars
if nargin < 1
    figName = NaN;
end
if nargin < 2
    mfile = [];
end
if nargin < 3
    figH = gcf;
end

% if figName not provided, get from figure
if isnan(figName)
    figName = get(figH,'Name');
end

% mfile is empty, replace it will calling function. If mfile is char, then
% assume that it is left empty intentionally
if isempty(mfile) && ~ischar(mfile)
    foo = dbstack;
    if length(foo) > 1
        mfile = foo(2).file;
        if ~endsWith(mfile,'.m')
            mfile = [mfile,'.m'];
        end
    end
end

figTitleHand = annotation(figH,'textbox',[0 0 1 .05]);
    figAnnotStr = {figName
        [datestr(now),', ',mfile]
        };    
set(figTitleHand,'String',figAnnotStr,'Interpreter','none','FontSize',12);
set(figTitleHand,'VerticalAlignment','top','HorizontalAlignment','right');
set(figTitleHand,'LineStyle','none');
set(figTitleHand,'Color','r');
