function varargout = plot(obj, hFig, varargin)
%PLOT Render data values of a matrix as a heat map.
%
%   PLOT(HMOBJ) renders a heat map of a HeatMap object HMOBJ, in a MATLAB
%   figure window.
%
%   PLOT(HMOBJ,HFIG) renders a heat map of a HeatMap object HMOBJ, in a
%   MATLAB figure window with the handle HFIG.
%
%   H = PLOT(...) returns the handle to the heat map axes. The graphic
%   properties are stored as application data in the axes handle.
%
%   Examples:
%
%       plot(hmobj)
%
%   See also HEATMAP, HEATMAP/VIEW.

%   Copyright 2009-2010 The MathWorks, Inc.


if isempty(obj)
    return;
end

%== Get figure handle
if nargin < 2
    hFig = [];
end

if nargin == 3 && ~isempty(varargin{1})
    if islogical(varargin{1}) % for g683204 plot colorbar
        showColorbar(hFig, bioinfoprivate.opttf(obj.Colorbar))
        return;
    else
        updateHMAxes(obj, varargin{1})
        hHMAxes = obj.HMAxesHandle;
        hFig = obj.FigureHandle;
    end
else
    if isempty(hFig) % Print to figure case
        hFig = figure('Renderer',     'ZBuffer',...
                      'Visible',      'on');
    end
    
    %== Use the specified colormap
    if ~isempty(obj.Colormap)
        set(hFig, 'Colormap',obj.Colormap);
    end
    
    %== Setup HeatMap axes
    hHMAxes = axes('Parent', hFig,...
                   'Units', 'normalized',...
                   'HitTest', 'off');
    
    %== Draw heat map image
    image(obj.Data, 'Parent',       hHMAxes,...
                    'Tag',          'HeatMapImage',...
                    'CDataMapping', 'Scaled');
    
    %== Refine the heatmap axes properties. Visible set to on by IMAGE.
    set(hHMAxes, 'Visible', 'off',...
                 'YDir',    'normal',...
                 'Layer',   'bottom',... 
                 'Tag',     'HeatMapAxes');
             
    %== Setup HeatMap Axes handle for view
    if ishandle(obj.FigureHandle)
        obj.HMAxesHandle = hHMAxes;
    end
        
    %== Scale the heat map
    scaleHeatMap(hHMAxes, obj)
    
    %== HeatMap Axes
    initHMAxes(obj, hHMAxes);
    
    %== Missing data values
    if obj.MissingDataFlag
        plotMissingData(obj, hHMAxes);
    end
    
    %== Setup resize and axes limit listeners
    setupPlotListeners(obj, hFig, hHMAxes);
 
    %== Display annotation text
    if strcmpi(obj.Annotate, 'on')
        createAnnotationText(obj, hHMAxes);
    end
end

if nargout >= 1
    varargout{1} = hHMAxes;
end

if nargout > 1
    varargout{2} = hFig;
end

end % end of plot

%== Helper functions ---------------------------------------
function tickVal = getTickValues(nTicks, tickLimit)
%Return the tick values if the number of ticks are within the limits. If
% over limit return zero.
if nTicks <= tickLimit
    tickVal = 1:nTicks;
else
    tickVal = [];
end
end
%-----------
function showColorbar(hHMAxes, showFlag)
% Insert Colorbar and modify the color bar context menu items.
hmdata = getappdata(hHMAxes, 'HeatMapAxesData');
if ishandle(hHMAxes) && showFlag
    hmPos = get(hHMAxes, 'Position');
    
    %== Create colorbar
    if isempty(hmdata.HMColorbarAxes) ||~ishandle(hmdata.HMColorbarAxes)
        hmdata.HMColorbarAxes = colorbar('peer',     hHMAxes,...
            'location', 'WestOutside',...
            'Tag',      'HeatMapColorbar');
    end
    
    cbarPos = get(hmdata.HMColorbarAxes, 'Position');
    if isempty(hmdata.HMColorbarWidth)
        hmdata.HMColorbarWidth = cbarPos(3);
    end
    switch hmdata.RowLabelsLocation
        case 'left'
            new_cbarPos = [hmPos(1)+ hmPos(3)+ hmdata.HMColorbarWidth*0.5
                cbarPos(2)
                hmdata.HMColorbarWidth*0.7
                cbarPos(4)];
            set(hmdata.HMColorbarAxes, 'YAxisLocation', 'right')
        case 'right'
            if isempty(hmdata.CBarLeftStart)
                startLeft = hmPos(1);
            else
                startLeft = hmdata.CBarLeftStart;
            end
            new_cbarPos = [startLeft-hmdata.HMColorbarWidth*0.95
                cbarPos(2)
                hmdata.HMColorbarWidth*0.7
                cbarPos(4)];
    end
    set(hmdata.HMColorbarAxes, 'Position', new_cbarPos)
    
    updateColorbarContextmenu(hmdata.HMColorbarAxes)
else
    colorbar('peer', hHMAxes,'off');
    hmdata.HMColorbarAxes = [];
    hmdata.HMColorbarWidth = [];
end

setappdata(hHMAxes, 'HeatMapAxesData', hmdata)
positionTickLabels(hHMAxes, 'X')
end

function updateColorbarContextmenu(cbar)
% Set colorbar context menu items "Visible' off.
uic = get(cbar, 'UIContextmenu');
set(findall(uic,'Type','uimenu','Tag','scribe:colorbar:location'), 'Visible', 'off');
set(findall(uic,'Type','uimenu','Tag','scribe:colorbar:colormap'), 'Visible', 'off');
set(findall(uic,'Type','uimenu','Tag','scribe:colorbar:interactivecolormapshift'), 'Visible', 'off');
set(findall(uic,'Type','uimenu','Tag','scribe:colorbar:editcolormap'), 'Visible', 'off');
set(findall(uic,'Type','uimenu','Tag','scribe:colorbar:propedit'), 'Visible', 'off');
set(findall(uic,'Type','uimenu','Tag','scribe:colorbar:mcode'), 'Visible', 'off');
end

function hRectangle = plotMissingData(obj, hHMAxes)
% Plot a gray rectangle over the image to indicate the missing data.

color = [0.8 0.8 0.8];
[ridx, cidx] = find(isnan(obj.OriginalData));

% Get the axes information
xl = get(hHMAxes, 'Xlim');
yl = get(hHMAxes, 'Ylim');

nMarkers = numel(ridx);
hRectangle = zeros(nMarkers, 1);

for i = 1:nMarkers
    recPos = [xl(1)+(cidx(i)-1), yl(1)+ridx(i)-1, 1, 1];
    hRectangle(i) = rectangle('Parent',    hHMAxes,...
                              'Position',  recPos,...
                              'FaceColor', color,...
                              'EdgeColor', color,...
                              'LineStyle', 'none');
    set(get(get(hRectangle(i),'Annotation'), 'LegendInformation'),...
        'IconDisplayStyle','on');
end
end
%-----------------
function initHMAxes(obj, hHMAxes)
% Update the HeatMap axes properties. Called by HeatMap axes related
% property listener and initial plot.

set(hHMAxes,'Fontsize',getFontSize(hHMAxes, obj.TickLimit, 4, true))
%== Determine the size of the image
[nYTicks, nXTicks] = size(obj.Data);

%= Y-axis ticks (rows)
ytickvals = getTickValues(nYTicks, obj.TickLimit); 

%= X-axis ticks (columns)
xtickvals = getTickValues(nXTicks, obj.TickLimit); 

%== Put tick labels on heatmap axes.
set(hHMAxes, 'XTick',         xtickvals,...
             'XTickLabel',    obj.ColumnLabels,...
             'YTick',         ytickvals,...
             'YTickLabel',    obj.RowLabels,...
             'TickLength',    [0 0],...
             'XAxisLocation', obj.ColumnLabelsLocation,...
             'YAxisLocation', obj.RowLabelsLocation);
         
hmdata = getappdata(hHMAxes, 'HeatMapAxesData');

%==
hmdata.TickLimit = obj.TickLimit; 
hmdata.HMTitleAxes = axes('Parent', get(hHMAxes, 'Parent'),...
                          'Position', [0 1 1 1],...
                          'visible', 'off',...
                          'Tag','HeatMapTitleAxes',...
                          'layer', 'bottom',...
                          'HitTest', 'off');
%== Update appdata for heatmap axes         
hmdata.RowLabelsLocation = obj.RowLabelsLocation;
hmdata.HMXTickLabelRot = obj.ColumnLabelsRotate;
hmdata.HMYTickLabelRot = obj.RowLabelsRotate;
%== Update title, axes labels
hmdata.HMAxisTitles = [get(hHMAxes, 'Xlabel'), get(hHMAxes, 'Ylabel')];
hmdata.HMTitleText = createHeatMapTitle(obj, hmdata.HMTitleAxes);
set(hmdata.HMAxisTitles(1), obj.XLabelPVPairs{:})
set(hmdata.HMAxisTitles(2), obj.YLabelPVPairs{:})

hmdata.XMarkerAxes = [];
hmdata.YMarkerAxes = [];
hmdata.XAnnotMarkerHandles = [];
hmdata.YAnnotMarkerHandles = [];

hmdata.HMXAnnotate = obj.ColumnLabelsColor;
hmdata.HMYAnnotate = obj.RowLabelsColor;
if obj.LabelsWithMarkers && ~isempty(obj.ColumnLabelsColor)
    hmdata = doAxisAnnotation(obj, hmdata, 'X', hHMAxes);
end
if obj.LabelsWithMarkers && ~isempty(obj.RowLabelsColor)
    hmdata = doAxisAnnotation(obj, hmdata, 'Y', hHMAxes);
end

hmdata.CBarLeftStart = [];
hmdata.HMColorbarAxes = [];
hmdata.HMColorbarWidth = [];

setappdata(hHMAxes, 'HeatMapAxesData', hmdata);

positionTickLabels(hHMAxes,'X');
positionTickLabels(hHMAxes,'Y');
positionAxes(obj, hHMAxes);
set(hHMAxes, 'Visible','on')
end

function updateHMAxes(obj, prop)
% Update the HeatMap axes properties. Called by HeatMap axes related
% property listener and initial plot.

hmdata = getappdata(obj.HMAxesHandle, 'HeatMapAxesData');
switch prop
    case 'ColumnLabels'
        hmdata = turnOffColorbar(obj.HMAxesHandle, hmdata);
        set(obj.HMAxesHandle, 'XTickLabel',obj.ColumnLabels);
        dim = 'X';
    case 'RowLabels'
        hmdata = turnOffColorbar(obj.HMAxesHandle, hmdata);
        set(obj.HMAxesHandle, 'YTickLabel',obj.RowLabels);
        dim = 'Y';
    case 'ColumnLabelsLocation'
        hmdata = turnOffColorbar(obj.HMAxesHandle, hmdata);
        set(obj.HMAxesHandle,'XAxisLocation', obj.ColumnLabelsLocation);
        dim = 'X';
    case 'RowLabelsLocation'
        hmdata = turnOffColorbar(obj.HMAxesHandle, hmdata);
        set(obj.HMAxesHandle,'YAxisLocation', obj.RowLabelsLocation);
        dim = 'Y';
    case 'ColumnLabelsRotate'
        hmdata.HMXTickLabelRot = obj.ColumnLabelsRotate;
        dim = 'X';
    case 'RowLabelsRotate'
        hmdata.HMYTickLabelRot = obj.RowLabelsRotate;
        dim = 'Y';
    case 'XLabel'
        set(hmdata.HMAxisTitles(1), obj.XLabelPVPairs{:})
        dim = 'X';
    case 'YLabel'
         set(hmdata.HMAxisTitles(2), obj.YLabelPVPairs{:})
        dim = 'Y';
    case 'Title'
        set(hmdata.HMTitleText, obj.TitlePVPairs{:})
        return;
    case 'ColumnLabelsColor'
        dim = 'X';
        hmdata = doAxisAnnotation(obj, hmdata, dim); 
    case 'RowLabelsColor'
        dim = 'Y';
        hmdata = doAxisAnnotation(obj, hmdata, dim);
    case 'LabelsWithMarkers'
            dim = ['X'; 'Y'];
            hmdata = doAxisAnnotation(obj, hmdata, dim(1));
            hmdata = doAxisAnnotation(obj, hmdata, dim(2));
    case 'Annotate'
        annotTextHandles = getappdata(obj.HMAxesHandle, 'HeatMapAnnotText');

        if isempty(annotTextHandles) || ~ishandle(annotTextHandles(1))
            createAnnotationText(obj, obj.HMAxesHandle);
            annotTextHandles = getappdata(obj.HMAxesHandle, 'HeatMapAnnotText');
        end
        
        set(annotTextHandles, 'visible', obj.Annotate);
        if strcmpi(obj.Annotate, 'on')
            resetAnnotTextFontSize(obj.HMAxesHandle)
        end
        return;
    case 'AnnotColor'
        annotTextHandles = getappdata(obj.HMAxesHandle, 'HeatMapAnnotText');
        if isempty(annotTextHandles) || ~ishandle(annotTextHandles(1))
            createAnnotationText(obj, obj.HMAxesHandle);
            annotTextHandles = getappdata(obj.HMAxesHandle, 'HeatMapAnnotText');
        end
        
        set(annotTextHandles, 'Color', obj.AnnotColor);
        return;
    case 'AnnotPrecision'
        createAnnotationText(obj, obj.HMAxesHandle);
        return;
    case 'Colormap'
        set(obj.FigureHandle, 'Colormap', obj.Colormap);
        return;
    case 'DisplayRange'
        scaleHeatMap(obj.HMAxesHandle, obj);
        return;
    case 'Symmetric'
        scaleHeatMap(obj.HMAxesHandle, obj);
        return;
    case 'Colorbar'
        showColorbar(obj.HMAxesHandle, bioinfoprivate.opttf(obj.Colorbar))
        return;
end
setappdata(obj.HMAxesHandle, 'HeatMapAxesData', hmdata);
if isempty(dim)
    return;
elseif size(dim, 1) == 1
    positionTickLabels(obj.HMAxesHandle,dim);
elseif size(dim, 1) == 2
    positionTickLabels(obj.HMAxesHandle,dim(1));
    positionTickLabels(obj.HMAxesHandle,dim(2));
end
%== For reposition axes if needed.
positionAxes(obj, obj.HMAxesHandle);
end

%--------------------------
function hmData = doAxisAnnotation(obj, hmData, dimFlag, hHMAxes)
if nargin < 4
    hHMAxes = obj.HMAxesHandle;
end
switch dimFlag
    case 'X'
        hmData.HMXAnnotate = obj.ColumnLabelsColor;
        annotStruct = obj.ColumnLabelsColor;
        mAxesName = 'XMarkerAxes';
        tLabels = 'ColumnLabels';
        annotAppdataLabel = 'XAnnotMarkerHandles';
        oAxesName = 'YMarkerAxes';
    case 'Y'
        hmData.HMYAnnotate = obj.RowLabelsColor;
        annotStruct = obj.RowLabelsColor;
        mAxesName = 'YMarkerAxes';
        tLabels = 'RowLabels';
        annotAppdataLabel = 'YAnnotMarkerHandles';
        oAxesName = 'XMarkerAxes';
end
        
if (isempty(annotStruct) && ~isempty(hmData.(mAxesName))) || ~obj.LabelsWithMarkers
    delete(hmData.(mAxesName))
    hmData.(mAxesName) = [];
    hmData.(annotAppdataLabel) = [];
elseif ~isempty(annotStruct) && isempty(hmData.(mAxesName)) && obj.LabelsWithMarkers
    hmData.(mAxesName) = createAnnotAxes(hHMAxes, obj.AnnotMarkerRatio, dimFlag); 
    %== Link Xlim and YLim properties for both axes
    hXAnnotLink = linkprop([hHMAxes hmData.(mAxesName) hmData.(oAxesName)], 'XLim');
    hYAnnotLink = linkprop([hHMAxes hmData.(mAxesName) hmData.(oAxesName)], 'YLim');
    setappdata(hHMAxes,'HeatMapXAnnotLimLink',hXAnnotLink )
    setappdata(hHMAxes,'HeatMapYAnnotLimLink', hYAnnotLink )
    
    hmData.(annotAppdataLabel) = plotAxesAnnotMarkers(hmData.(mAxesName),...
                                    obj.(tLabels),  annotStruct, dimFlag);
end
end

function scaleHeatMap(hHMAxes, obj)
%SCALEHEATMAP Update the CLIM in image axes
if obj.Symmetric
    maxval = min(max(abs(obj.Data(:))), obj.DisplayRange);
    minval = -maxval;
else
    maxval = min(max(obj.Data(:)), obj.DisplayRange);
    minval = min(obj.Data(:));
end
set(hHMAxes, 'Clim', [minval,maxval]);
end

function createAnnotationText(obj, hHMAxes)
% Create data value annotation text handles and save them in HMAxes
% appdata.

annotTextHandles = getappdata(hHMAxes, 'HeatMapAnnotText');
if ~isempty(annotTextHandles) && ishandle(annotTextHandles(1,1))
    delete(annotTextHandles)
end
[nY, nX] = size(obj.Data);

fontSize = getFontSize(hHMAxes, obj.TickLimit, 4);
if fontSize == 0
    obj.Annotate = 'off';
    fontSize = 9;
end

%== Display annotation text
annotTextHandles = zeros(nY, nX)';
for i = 1:nY
    for j = 1:nX
        annotTextHandles(j,i) = text(j, i, num2str(obj.Data(i,j), obj.AnnotPrecision),...
            'Parent', hHMAxes,...
            'Color', obj.AnnotColor,...
            'HorizontalAlignment', 'Center',...
            'FontSize', fontSize,...
            'Clipping','on',...
            'HitTest', 'off',...
            'Visible','off');
    end
end

%== Save handles into appdata of HMAxes
setappdata(hHMAxes, 'HeatMapAnnotText', annotTextHandles);
end
%--------------------------
function annotAxes = createAnnotAxes(hHMAxes, markerRatio, dim)
% Create axes for axis annotation markers, and link the axis to the heatmap
% Axes, markerRatio is the ratio between width and height of annotation
% marker axes and heatmap axes.
hmPos = get(hHMAxes, 'Position');
dim = lower(dim);
switch dim
    case 'x'
        aPos = [hmPos(1)
            hmPos(2)-hmPos(4)/markerRatio
            hmPos(3)
            hmPos(4)/markerRatio];
    case 'y'
        aPos = [hmPos(1)- hmPos(3)/markerRatio
            hmPos(2)
            hmPos(3)/markerRatio
            hmPos(4)];
end
annotAxes = axes('Parent', get(hHMAxes, 'Parent'),...
                 'Position', aPos,...
                 'HitTest', 'off');
axis(annotAxes, 'off')
set(hggetbehavior(annotAxes,'Zoom'),'Enable',false) ;
linkaxes([hHMAxes annotAxes], dim);
end

%---------
function fs = getFontSize(imAxes, magicNumber, minLim, noZeroFlag)
%SETFONTSIZE computes the appropriate font size for axes tick labels
% and annotation text.
% Try to keep font size reasonable for x and y tick labels.
if nargin < 4
    noZeroFlag = false;
end
hFig = get(imAxes,'Parent');
if nargin == 1
    magicNumber = 200;
    minLim = 4;
elseif nargin == 2
    minLim = 4;
end

nrows = diff(get(imAxes,'YLim'));
ncols = diff(get(imAxes,'XLim'));
if ncols < magicNumber && nrows < magicNumber
    ratio = max(get(hFig,'Position').*[0 0 0 1])/max(nrows,ncols);
elseif ncols < magicNumber
    ratio = max(get(hFig,'Position').*[0 0 0 1])/ncols;
elseif nrows < magicNumber
    ratio = max(get(hFig,'Position').*[0 0 0 1])/nrows;
else
    ratio = 1;
end

fs = min(9, ceil(ratio/1.5));    % the gold formula
if fs < minLim
    fs = 0;
end

if noZeroFlag && fs == 0
   fs = 9 ;
end
end
%------------
function resetAnnotTextFontSize(hHMAxes)
%== Annotation text font size
annotTextHandles = getappdata(hHMAxes, 'HeatMapAnnotText');
if ~isempty(annotTextHandles) && ishandle(annotTextHandles(1,1))
    fs = getFontSize(hHMAxes);
    if fs > 0
        set(annotTextHandles,'FontSize',fs,...
                             'Visible','on');
    else
        set(annotTextHandles,'Visible','off');
    end
end
end
%------------

function setupPlotListeners(obj, hFigure, hHMAxes)
% Sets up listeners for resizing, axes limit changes.
resizeListener = addlistener(hFigure, 'SizeChange',...
                             @(src,evt) resizeListenerCB(obj, hHMAxes,src, evt));

xLimListener = addlistener(hHMAxes, 'XLim', 'PostSet', ...
                         @(src,evt) hmAxisLimitListenerCB(obj, src, evt));

yLimListener = addlistener(hHMAxes, 'YLim', 'PostSet', ...
                         @(src,evt) hmAxisLimitListenerCB(obj, src, evt));
                     
positionListener = addlistener(hHMAxes, 'Position', 'PostSet', ...
                         @(src,evt) hmPositionListenerCB(src, evt));
                     
cameraViewAngleModeListener = addlistener(hHMAxes, 'CameraViewAngleMode', 'PostSet', ...
                         @(src,evt) hmCameraViewAngleModeListenerCB(src, evt));
                     
listeners = getappdata(hFigure,'HeatMapListeners');
if isempty(listeners)
    listeners={};
end

setappdata(hFigure,'HeatMapListeners',...
[listeners {resizeListener xLimListener yLimListener positionListener cameraViewAngleModeListener}]);
end
%-----------------
function resizeListenerCB(obj, hHMAxes, src, ~)
% Callback for figure window resize
if ~ishandle(hHMAxes)
    children = get(src, 'children');
    hmIdx = strcmpi(get(children, 'Tag'), 'HeatMapAxes');
    hHMAxes = children(hmIdx);
end
if ~isempty(hHMAxes) && ishandle(hHMAxes)
    if strcmpi(obj.Annotate, 'on')
        resetAnnotTextFontSize(hHMAxes);
    end
    positionTickLabels(hHMAxes, 'X')
    positionTickLabels(hHMAxes, 'Y')
    
    positionAxes(obj, hHMAxes)
end
end

%--------------
function hmAxisLimitListenerCB(obj, hSrc,event)
% obj - HeatMap object.

hHMAxes = event.AffectedObject;
set(hHMAxes,'Fontsize',getFontSize(hHMAxes, obj.TickLimit, 4, true))
%== Annotation text font size
if strcmpi(obj.Annotate, 'on')
    resetAnnotTextFontSize(hHMAxes);
end
%== Figure out new limits
switch hSrc.Name
    case 'XLim'
        newLim = get(event.AffectedObject,'XLim');
        if diff(newLim) <= obj.TickLimit
            colLabels = obj.ColumnLabels;
            xtickvals = max(1,ceil(newLim(1))): min(floor(newLim(2)),numel(colLabels));
            xticklabels = colLabels(xtickvals);
        else
            xtickvals = [];
            xticklabels = '';
        end
        set(hHMAxes,'XTick',xtickvals,'XTickLabel',xticklabels);
    case 'YLim'
        newLim = get(event.AffectedObject,'YLim');
        if diff(newLim) <= obj.TickLimit
            rowLabels = obj.RowLabels;
            ytickvals = max(1,ceil(newLim(1))):min(floor(newLim(2)),numel(rowLabels));
            yticklabels =  rowLabels(ytickvals);
        else
            ytickvals = [];
            yticklabels = '';
        end
        set(hHMAxes,'YTick',ytickvals,'YTickLabel',yticklabels);
end
positionAxes(obj, hHMAxes)
positionTickLabels(hHMAxes, 'X', true);
positionTickLabels(hHMAxes, 'Y', true);

annotTextHandles = getappdata(hHMAxes, 'HeatMapAnnotText');
if ~isempty(annotTextHandles) && ishandle(annotTextHandles(1))
    set(annotTextHandles, 'visible', obj.Annotate);
end
end

function hmPositionListenerCB(hSrc,event)
% obj - HeatMap object.
hHMAxes = event.AffectedObject;
hmdata = getappdata(hHMAxes, 'HeatMapAxesData');
if strcmpi(hSrc.Name, 'Position')
    if ~isempty(hmdata.YMarkerAxes)
        hmPos = get(hHMAxes, 'Position');
        maPos = get(hmdata.YMarkerAxes, 'Position');
        set(hmdata.YMarkerAxes, 'Position', [hmPos(1)+hmPos(2)+2
            hmPos(2)
            maPos(3)
            hmPos(4)])
        positionTickLabels(hHMAxes, 'Y')
    end
    
    if ~isempty(hmdata.XMarkerAxes)
        hmPos = get(hHMAxes, 'Position');
        maPos = get(hmdata.XMarkerAxes, 'Position');
        set(hmdata.XMarkerAxes, 'Position', [hmPos(1)
            hmPos(2)-maPos(4)-2
            hmPos(3)
            maPos(4)])
        positionTickLabels(hHMAxes, 'X')
    end
end
end

function hmCameraViewAngleModeListenerCB(~, event)
% HMCAMERAVIEWANGLEMODELISTENER - Add the listener as a workaround for when
% reset view, the CAMERAVIEWANGEMODELISTENER property changed from 'manual'
% to 'auto', and did not trigger reset the axis labels.
hHMAxes = event.AffectedObject;
positionTickLabels(hHMAxes, 'X', true)
positionTickLabels(hHMAxes, 'Y', true)
end
%------------
function positionTickLabels(haxes, axisType, limFlag)
% POSITIONTICKLABEL Position tick labels, title and axes labels.
%
% POSITIONTICKLABEL(HAXES, AXISTYPE, LIMFLAG) position the specified axis
% (AXISTYPE=X or Y) tick labels, axis labels. If LIMFLAG is true, not need
% to recompute the positions.

if nargin < 3
    limFlag = false;
end

% Need to check if is called while HG printdlg, if so nothing is updated
% (G591707)
dbs = dbstack;
if any(strcmp({dbs.file},'printdlg.m'))
    return;
end

hmdata = getappdata(haxes, 'HeatMapAxesData');
switch axisType
    case 'X' % X-axis
        ini = 1; % Along the direction to be deal with
        iniSpan = 3;
        fix  = 2; % Along the direction no need to be changed
        fixSpan = 4;
        if strcmpi(get(haxes, 'XAxisLocation'), 'bottom')
            posDir = 2; % bottom
            textHorzAlign = {'right', 'left'};
        else
            posDir = 4; %top
            textHorzAlign = {'left', 'right'};
        end
        tickType = 'XTick';
        tickLabelType = 'XTickLabel';
        axisLim = 'XLim';
        rot = hmdata.HMXTickLabelRot;
        axisDim = 1;
        tickAppdataLabel = 'XTickLabelTextHandles';
        annotStruct = hmdata.HMXAnnotate;
        markerAxes = hmdata.XMarkerAxes;
        o_markerAxes = hmdata.YMarkerAxes;
    case 'Y' %Yaxis
        ini = 2;
        iniSpan = 4;
        fix  = 1; % Along the direction no need to be changed
        fixSpan = 3;
        if strcmpi(get(haxes, 'YAxisLocation'), 'left')
            posDir = 1; %left
            textHorzAlign = {'right','left'};
        else
            posDir = 3; %right
            textHorzAlign = {'left','right'};
        end
        
        tickType = 'YTick';
        tickLabelType = 'YTickLabel';
        axisLim = 'YLim';
        rot = hmdata.HMYTickLabelRot;
        axisDim = 2;
        tickAppdataLabel = 'YTickLabelTextHandles';
        annotStruct = hmdata.HMYAnnotate;
        markerAxes = hmdata.YMarkerAxes;
        o_markerAxes = hmdata.XMarkerAxes;
end

hAxisTitles = hmdata.HMAxisTitles;

%== Getting axis title handles
textUnit_orig = get(hAxisTitles(1), 'unit');
set(hAxisTitles, 'Unit', 'points')

%== Get tick labels
oldHTickText = getappdata(haxes, tickAppdataLabel);
tickLabels_orig = [];
if ishandle(oldHTickText)
    tickLabels_orig = get(oldHTickText, 'String');
    %== Remove the old text object
    delete(oldHTickText);
    setappdata(haxes, tickAppdataLabel, [])
end

if diff(get(haxes, axisLim)) > hmdata.TickLimit
    return;
end

%== Determine real axis tick value
ticks = get(haxes,tickType);
tickLabels = get(haxes,tickLabelType);

if isempty(tickLabels)
    tickLabels = tickLabels_orig;
end

if isempty(tickLabels)
    return;
end

if ~isempty(tickLabels) && isempty(ticks)
    ticks = getTickValues(numel(tickLabels), hmdata.TickLimit); 
    if isempty(ticks)
        return;
    end
end

%== Add new xTickLabel text object
if axisDim == 1
    x_tick = ticks;
    y_tick = ones(length(ticks), 1);
else
    x_tick = ones(length(ticks), 1);
    y_tick = ticks;
end

hTickText = text(x_tick, y_tick,...
    tickLabels,...
    'Parent', haxes,...
    'Clipping', 'off',...
    'Interpreter', 'none');
%== Rotate the tick text object and set initial font size
tickLabelProps = {'Rotation',rot,...
        'Unit', 'points',...
        'FontUnit', 'points',...
        'HitTest', 'off',...
        'FontSize', 9};
if rot<180
    set(hTickText,'HorizontalAlignment',textHorzAlign{1},...
                  'VerticalAlignment', 'middle',...
                  tickLabelProps{:});
else
    set(hTickText,'HorizontalAlignment',textHorzAlign{2},...
                  'VerticalAlignment', 'top',...
                  tickLabelProps{:});
end

%== Get the size of the text box
if length(hTickText) > 1
    textSizePts = cell2mat(get(hTickText,'extent'));
    textPosPts = cell2mat(get(hTickText,'position'));
else
    textSizePts = get(hTickText,'extent');
    textPosPts = get(hTickText,'position');
end

%== Get the longest value
longestTextPts =  max(textSizePts(:,fixSpan));

%== Get title and axis labels in points
titlesSizePts = cell2mat(get(hAxisTitles, 'extent'));
titlesPosPts = cell2mat(get(hAxisTitles, 'position'));

if ~limFlag
    %== Get axis position values
    looseInset = get(haxes, 'LooseInset');
    defLooseInset = get(get(haxes, 'parent'), 'DefaultAxesLooseInset');
    oPos = get(haxes,'OuterPosition');
    
    if isempty(markerAxes)
        markerPos = zeros(1,4);
    else
        markerPos = get(markerAxes, 'Position');
    end
    
    %== Convert axis position values to points
    looseInsetPts = hgconvertunits(ancestor(haxes, 'Figure'), looseInset,...
        get(haxes, 'Unit'), 'points', get(haxes, 'parent'));
    
    defLooseInsetPts = hgconvertunits(ancestor(haxes, 'Figure'), defLooseInset,...
        get(haxes, 'Unit'), 'points', get(haxes, 'parent'));
    
    oPosPts = hgconvertunits(ancestor(haxes, 'Figure'), oPos,...
        get(haxes, 'Unit'), 'points', get(haxes, 'parent'));
    
    markerPosPts = hgconvertunits(ancestor(haxes, 'Figure'), markerPos,...
        get(haxes, 'Unit'), 'points', get(haxes, 'parent'));
    
    %== Determine how much room is available for tick labels, assuming the
    %axes are squished to zero width.
    maxAllowPts = oPosPts(fixSpan) -...
                  defLooseInsetPts(fix) - ...
                  defLooseInsetPts(fixSpan) -...
                  markerPosPts(fixSpan);
    
    if axisDim == 1
        if (longestTextPts + titlesSizePts(axisDim, fixSpan) + ...
            markerPosPts(fixSpan)) > maxAllowPts    
            longestTextPts = maxAllowPts - titlesSizePts(axisDim, fixSpan)...
                -markerPosPts(fixSpan);
        end
    else
        if (longestTextPts + titlesSizePts(axisDim, fixSpan)...
                + markerPosPts(fixSpan)) > maxAllowPts
            longestTextPts = maxAllowPts -...
                             titlesSizePts(axisDim, fixSpan) -...
                             markerPosPts(fixSpan);
        end
    end
    
    switch posDir
        case 1 %left
            newLIPts = getNewLooseInsetPoints(defLooseInsetPts, looseInsetPts,...
                longestTextPts, titlesSizePts(axisDim, fixSpan),...
                markerPosPts(fixSpan), ini, iniSpan, fix, fixSpan);
        case 2 %bottom
            newLIPts = getNewLooseInsetPoints(defLooseInsetPts, looseInsetPts,...
                longestTextPts, titlesSizePts(axisDim, fixSpan),...
                markerPosPts(fixSpan), ini, iniSpan, fix, fixSpan);
        case 3 %right
            newLIPts = getNewLooseInsetPoints(defLooseInsetPts, looseInsetPts,...
                longestTextPts, titlesSizePts(axisDim, fixSpan),...
                markerPosPts(fixSpan), ini, iniSpan, fixSpan, fix);
            
        case 4 %top with titles and xlabel
            newLIPts = getNewLooseInsetPoints(defLooseInsetPts, looseInsetPts,...
                longestTextPts, titlesSizePts(axisDim, fixSpan),...
                markerPosPts(fixSpan), ini, iniSpan, fixSpan, fix);
    end
    
    %==Convert back to original axes units.
    newLI = hgconvertunits(ancestor(haxes,'figure'), newLIPts, ...
        'points', get(haxes,'Units'),get(haxes,'Parent'));
    
    if newLI(2)+newLI(4)>oPos(4)
        newLI(2) = oPos(4)-newLI(4);
    end
    set(haxes, 'LooseInset', newLI);
    
    %== Get the updated axes position.
    axisPos = get(haxes,'Position');
    axisPosPts = hgconvertunits(ancestor(haxes,'figure'), axisPos, ...
        get(haxes,'Units'),'points',get(haxes,'Parent'));
    
    %== Axis annotation marker
    markerNewPts = zeros(1,4);
    if ~isempty(annotStruct)
        tnLabels = nominal(tickLabels);
        if isempty(markerAxes)
            for iloop = 1:numel(annotStruct.Labels)
                idx = tnLabels == strtrim(annotStruct.Labels{iloop});
                set(hTickText(idx), 'color', annotStruct.Colors{iloop})
            end
        else %by color the text only
            switch posDir
                case 1
                    markerNewPts = [axisPosPts(1)-markerPosPts(3)-2
                        axisPosPts(2)
                        markerPosPts(3)
                        axisPosPts(4)];
                case 2
                    markerNewPts = [axisPosPts(1)
                        axisPosPts(2)-markerPosPts(4)-2
                        axisPosPts(3)
                        markerPosPts(4)];
                case 3
                    markerNewPts = [axisPosPts(1)+axisPosPts(3)+2
                        axisPosPts(2)
                        markerPosPts(3)
                        axisPosPts(4)];
                case 4
                    markerNewPts = [axisPosPts(1)
                        axisPosPts(2)+axisPosPts(4)+2
                        axisPosPts(3)
                        markerPosPts(4)];
            end
            newMarker = hgconvertunits(ancestor(haxes,'figure'), markerNewPts, ...
                'points', get(haxes,'Units'),get(haxes,'Parent'));
            set(markerAxes, 'position', newMarker);
        end
    end
    
    if ~isempty(o_markerAxes) && ishandle(o_markerAxes)
        mpos = get(o_markerAxes, 'Position');
        
        if posDir == 1 || posDir == 3
            mpos(1) = axisPos(1);
            mpos(3) = axisPos(3);
        elseif posDir == 2 || posDir == 4
            mpos(2) = axisPos(2);
            mpos(4) = axisPos(4);
        end
        set(o_markerAxes, 'Position', mpos);
    end
else
    axisPos = get(haxes,'Position');
    axisPosPts = hgconvertunits(ancestor(haxes,'figure'), axisPos, ...
        get(haxes,'Units'),'points',get(haxes,'Parent'));
    if isempty(markerAxes)
        markerNewPts = zeros(1,4);
    else
        markerPos = get(markerAxes, 'Position');
        markerNewPts = hgconvertunits(ancestor(haxes, 'Figure'), markerPos,...
            get(haxes, 'Unit'), 'points', get(haxes, 'parent'));
    end
end

%== Position tick labels
if posDir == 1 || posDir == 2
    labelShiftPts = -3-markerNewPts(fixSpan);
elseif posDir == 3 || posDir == 4
    labelShiftPts = axisPosPts(fixSpan)+markerNewPts(fixSpan) + 3;
end

textPosPts(:,3) = -0.1; % Slight negative Z, so datatip appears in front.
if axisDim == 1
    textPosPts(:,2) = labelShiftPts;
else
    textPosPts(:,1) = labelShiftPts;
end

labelPos = num2cell(textPosPts, 2);

% Set the label positions.  Be sure to specify points units, as
% the units may have changed (for example, during printing).
set(hTickText,{'Position'},labelPos);

%== Haxes clipping is off, turn off the tick labels outside axis range
set(hTickText,'Unit', 'Data');
if length(hTickText) > 1
    tt = cell2mat(get(hTickText, 'position'));
else
    tt = get(hTickText, 'position');
end
alim = get(haxes, axisLim);
idx = tt(:,ini)<=alim(2) & tt(:,ini)>= alim(1);
set(hTickText(~idx), 'visible', 'off')

%== Titles
switch posDir
    case 1
        titlesPosPts(2,1) = -titlesSizePts(2, 3)-longestTextPts-markerNewPts(3);
        titlesPosPts(2,2) = axisPosPts(iniSpan)/2;
        
        titlesPosPts(1,1) = axisPosPts(fixSpan)/2;
    case 2
        titlesPosPts(1,1) = axisPosPts(iniSpan)/2;
        
        titlesPosPts(1,2) = -longestTextPts - titlesSizePts(1, 4) - markerNewPts(4);
        titlesPosPts(2,2) = axisPosPts(fixSpan)/2;
    case 3
        titlesPosPts(1,1) = axisPosPts(fixSpan)/2;
        % Added for G575559 - 08/27/2009
        titlesPosPts(1,2) = -longestTextPts - titlesSizePts(1, 4) - markerNewPts(4);
        
        titlesPosPts(2,1) =  axisPosPts(fixSpan)+ longestTextPts + ...
                             titlesSizePts(2, 3)+markerNewPts(3);
        titlesPosPts(2,2) = axisPosPts(4)/2;
    case 4
        titlesPosPts(1,1) = axisPosPts(3)/2;
        titlesPosPts(1,2) = axisPosPts(fixSpan) + longestTextPts + ...
                            titlesSizePts(1, 4)+markerNewPts(4);
        
        titlesPosPts(2,2) = axisPosPts(iniSpan)/2;
end

titlesPosPts(:,3) = 0.1;
titlePos = num2cell(titlesPosPts, 2);
set(hAxisTitles, {'Position'}, titlePos);

%== Reset to original units
set(hAxisTitles, 'Unit', textUnit_orig);
set(haxes, tickLabelType, '')
setappdata(haxes, tickAppdataLabel, hTickText);
end

function newLIPts = getNewLooseInsetPoints(dfLIPts, lIPts, tickLabelPts,...
                        titleLabelPts, markePts, ini, iniSpan, fix, fixSpan)
% Return adjusted loose inset in points unit. All input much be in points
% unit as well.
% dfLIPts - DefaultLooseInset
% lIPts - LooseInset
% tickLabelPts - Longest tick labels in points
% titleLabelPts - Axis title labels (XLabel, YLabel, titles)

newLIPts([ini iniSpan]) = lIPts([ini iniSpan]);

if tickLabelPts > dfLIPts(fix)-titleLabelPts-markePts
    newLIPts(fix) = dfLIPts(fix) + tickLabelPts + titleLabelPts+markePts;
else
    newLIPts(fix) = dfLIPts(fix);
end

% % newLIPts(fixSpan) = lIPts(fixSpan);
newLIPts(fixSpan) = dfLIPts(fixSpan);
end

function hRecs = plotAxesAnnotMarkers(haxes, tLabels,  annotStruct, dim)
% Plot the axes annotation markers with annotation.
% The rectangle is plotted in data units of axes haxes.

colors = tLabels;
for iloop = 1:numel(annotStruct.Labels)
    idx = nominal(tLabels) == strtrim(annotStruct.Labels{iloop});
    colors(idx) =  annotStruct.Colors(iloop);
    if iloop == 1
        colors(~idx) = {'b'};
    end
end

%==Get the axes information
aXLim = get(haxes, 'Xlim');
aYLim = get(haxes, 'Ylim');
width = diff(aXLim);
height = diff(aYLim);
switch dim
    case 'X'
        aNum = fix(diff(aXLim));
        xStart = aXLim(1):aXLim(2)-1;
        yStart = ones(1, aNum)*aYLim(1);
    case 'Y'
        aNum = fix(diff(aYLim));
        xStart = ones(1, aNum)*aXLim(1);
        yStart = aYLim(1):aYLim(2)-1;
end

hRecs = zeros(aNum, 1);
for i = 1:aNum
    rPos = [xStart(i) yStart(i) width height];
    hRecs(i) = rectangle('Parent', haxes, 'Position', rPos,...
        'FaceColor', colors{i},...
        'LineStyle', 'none',...
        'HitTest', 'off');
    set(get(get(hRecs(i),'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
end
end

%-----------
function hmdata = turnOffColorbar(hHMAxes, hmdata)
%== Turn off colorbar
colorbar('peer', hHMAxes,'off');
hmdata.HMColorbarAxes = [];
hmdata.HMColorbarWidth = [];
end

%--------
function varargout = createHeatMapTitle(obj, HMTitleAxes)             
fs = get(get(HMTitleAxes, 'Parent'), 'defaultaxesfontsize') + 4;
titleXPos = 0.5;
titleYPos = 0.95;             
ht = text(titleXPos, titleYPos-1, obj.TitlePVPairs{2}, 'Parent', HMTitleAxes);
set(ht, 'HorizontalAlignment', 'center',...
        'fontsize',fs);

if length(obj.TitlePVPairs) > 2
   set(ht, obj.TitlePVPairs{3:end}) 
end

if nargout > 0
    varargout{1} = ht;
end
end
