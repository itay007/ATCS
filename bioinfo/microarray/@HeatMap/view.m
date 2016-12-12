function view(obj, varargin)
%VIEW View a HeatMap object in a MATLAB figure window.
%
%   VIEW(HMOBJ) shows a HeatMap object in a MATLAB figure window.
%
%   Examples:
%       hmap = HeatMap(gallery('invhess',20));
%       view(hmap)
%
%   See also HEATMAP, HEATMAP/PLOT.

%   Copyright 2009-2011 The MathWorks, Inc.


if isempty(obj)
    return;
end
figureTag = 'HeatMap';
if nargin > 1
    arg = varargin{1};
    if ishandle(arg) && ~strncmpi(get(arg,'Tag'), figureTag, numel(figureTag)) 
        obj.FigureHandle = varargin{1};
    elseif ischar(arg)
       plot(obj, obj.FigureHandle, arg); 
    end
end

if isempty(obj.FigureHandle) || ~ishandle(obj.FigureHandle)
    figureName = bioinfoprivate.indexedFigureName('HeatMap', 'HeatMap');
    
    obj.FigureHandle = figure('Renderer',         'zbuffer',...
                              'Name',             figureName,...
                              'NumberTitle',      'off',...
                              'Tag',              figureTag,...
                              'IntegerHandle',    'off',...
                              'HandleVisibility', 'callback',...
                              'Visible',          'on');
else
    disableUIControls(obj.FigureHandle);
    rmappdata(obj.FigureHandle,'HeatMapListeners');
end

%== Update toolbar and menu items
updateToolbar(obj);
updateUIMenus(obj);

%== Plot the clustergram;
delete(findall(obj.FigureHandle, 'Type', 'axes'))
plot(obj, obj.FigureHandle);
%== Handle data cursor on heatmap
updateFigureModes(obj, @heatmapDataCursorCB)
end

function updateToolbar(obj)
% helper function to update the toolbar
% Remove save, open edit plot buttons

oldSH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

set(obj.FigureHandle,'toolbar','figure')  % needs to update because uicontrols turn it off

% Fix toolbar options, we keep: ZoomIn,ZoomOut,Pan
hw = findall(obj.FigureHandle, 'type','uitoolbar');
hf = get(hw,'Children');
h1 = findall(hf,'Tag','Exploration.Pan');
h2 = findall(hf,'Tag','Exploration.ZoomOut');
h3 = findall(hf,'Tag','Exploration.ZoomIn');
% % h3 = findall(hf,'Tag','Annotation.InsertColorbar');
h4 = findall(hf,'Tag','Exploration.DataCursor');
h5 = findall(hf,'Tag','Standard.PrintFigure');
delete(setxor(hf,[h1,h2,h3,h4, h5]))

iconfile = fullfile(toolboxdir('matlab'),'icons','tool_colorbar.gif');
hColorbar = bioinfoprivate.createToolbarButton(hw, 2, iconfile,...
    'Insert Colorbar', 'HMInsertColorbar', 'State', obj.Colorbar, 'separator', 'on');
set(hColorbar, 'ClickedCallback', {@insertColorbarCB, obj})

iconfile = fullfile(toolboxdir('bioinfo'),'proteins','icons','reset_view.gif');
hAnnot = bioinfoprivate.createToolbarButton(hw, 2, iconfile,...
    'Annotate', 'HMAnnotateText', 'State', obj.Annotate);
set(hAnnot, 'ClickedCallback', {@annotateTextCB, obj})

set(0,'ShowHiddenHandles',oldSH)
end
%---------------------------------------------------------
function updateUIMenus(obj)
% helper function to set UI menus
% Remove File->Save, Save as etc. menu items

oldSH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

%== Delete figure menus not used
%h1 = findall(obj.FigureHandle,'Type','uimenu', 'Label','&Edit');
h1 = findall(obj.FigureHandle,'Type','uimenu', 'Tag','figMenuEdit');
%h2 = findall(obj.FigureHandle,'Type','uimenu', 'Label','&View');
h2 = findall(obj.FigureHandle,'Type','uimenu', 'Tag','figMenuView');
%h3 = findall(obj.FigureHandle,'Type','uimenu', 'Label','&Insert');
h3 = findall(obj.FigureHandle,'Type','uimenu', 'Tag','figMenuInsert');
delete([h1,h2, h3])

%== Repair "File" menu
%hw = findall(obj.FigureHandle,'Type','uimenu', 'Label','&File');
hw = findall(obj.FigureHandle,'Type','uimenu', 'Tag','figMenuFile');
hf = get(hw,'children');
%h1 = findall(hw,'Label','Expo&rt Setup...');
h1 = findall(hw,'Tag','figMenuFileExportSetup');
%h2 = findall(hw,'Label','Print Pre&view...');
h2 = findall(hw,'Tag','figMenuFilePrintPreview');
%h3 = findall(hw,'Label','&Print...');
h3 = findall(hw,'Tag','printMenu');
delete(setxor(hf,[h1,h2,h3]))

uimenu(hw,'Label','Exit','Separator','on',...
    'Position',4,'Callback', @closeFigure)
set(h1,'Separator','on')

%==Repair "Tools" menu
%hw = findall(obj.FigureHandle,'Type','uimenu','Label','&Tools');
hw = findall(obj.FigureHandle,'Type','uimenu','Tag','figMenuTools');
hf = get(hw,'children');
h1 = findall(hw,'Tag','figMenuZoomIn');
h2 = findall(hw,'Tag','figMenuZoomOut');
h3 = findall(hw,'Tag','figMenuResetView');
h4 = findall(hw,'Tag','figMenuOptions');
h5 = findall(hw,'Tag','figMenuDatatip');
set([h1,h3],'separator','off')
delete(setxor(hf,[h1,h2,h3,h4,h5]))

% Add Colorbar switch
h6 = bioinfoprivate.createMenuItem(obj.FigureHandle, hw, 'Insert Colorbar',...
                                   'HMInsertColorbar',...
                                   'Position', 6,...
                                   'Checked', obj.Colorbar);
set(h6,'Separator','on')
setappdata(obj.FigureHandle, 'HMMenuItems', h6);                              
% Add Annotation switch
h7 = bioinfoprivate.createMenuItem(obj.FigureHandle, hw, 'Annotate',...
                                   'HMAnnotateText',...
                                   'Position', 7,...
                                   'Checked', obj.Annotate);
setappdata(obj.FigureHandle, 'HMMenuItems', h7);

% Repair "Help" menu
bioinfoprivate.bioFigureHelpMenu(obj.FigureHandle, 'HeatMap', 'HeatMap_refpage')
set(0,'ShowHiddenHandles',oldSH)
end

function closeFigure(h, evt) %#ok
%== Remove all axes first
hFig = gcbf;
delete(findall(hFig, 'Type', 'axes'))
close (hFig)
end

%---------------------------------------
%== Data cursor for heatmap
function datacursorLabel = heatmapDataCursorCB(hsrc, event, dcManager, obj) %#ok
% Display data cursor in heatmap image.
datacursorLabel = [];
if strcmpi(get(dcManager, 'SnapToDataVertex'), 'off')
    return;
end

tg = get(event, 'Target');
if strcmpi(get(tg, 'Type'), 'image')
    type = get(tg, 'Tag');
    if strcmpi(type, 'HeatMapImage')
        pos = get(event, 'Position');
        cdata = get(tg, 'CData');
        val = cdata(pos(2), pos(1));
        
        origval = obj.OriginalData(pos(2), pos(1));
        scale = [];
        if ~isempty(obj.Scales)
            switch obj.Standardize
                case 'x'
                    scale = obj.Scales(pos(1));
                case 'y'
                    scale = obj.Scales(pos(2));
            end
        end
        
        if isempty(scale)
            valstr=sprintf('%0.2f(Value:%0.2f)', val, origval);
        else
            valstr=sprintf('%0.2f(Value:%0.2f,Std:%0.2f)', val, origval, scale);
        end
        if isempty(obj.RowLabels)
            rowStr = num2str(pos(2));
        else
            rowStr = obj.RowLabels{pos(2)};
        end
        
        if isempty(obj.ColumnLabels)
            colStr = num2str(pos(1));
        else
            colStr = obj.ColumnLabels{pos(1)};
        end
        datacursorLabel = {valstr, rowStr, colStr};
    else
        datacursorLabel = [];
        return;
    end
else
    datacursorLabel = [];
    return;
end
end
%----------------------------------------------------------
function disableUIControls(hFig)
activateuimode(hFig, '')
end

%--------------
function annotateTextCB(hSrc, hEvt, obj) %#ok
hFig= gcbf;
state = bioinfoprivate.toggleState(hFig, hSrc);

switch state
    case 'on'
        obj.Annotate = true;
    case 'off'
        obj.Annotate = false;
end
end % end of function

function insertColorbarCB(hSrc, ~, obj)
hFig= gcbf;
state = bioinfoprivate.toggleState(hFig, hSrc);

switch state
    case 'on'
        obj.Colorbar = true;
    case 'off'
        obj.Colorbar = false;
end
end % end of function

