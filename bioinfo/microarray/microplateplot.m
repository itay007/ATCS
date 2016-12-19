function outAxes = microplateplot(data,varargin)
%MICROPLATEPLOT displays a visualization of a microtiter plate
%
%   MICROPLATEPLOT(DATA) displays an image of a microtiter plate with each
%   well colored according to the values in DATA. DATA should be a matrix
%   or DataMatrix object.
%
%   MICROPLATEPLOT(...,'ROWLABELS',ROWLABELS) allows you to specify the
%   labels for the rows. The default values are the first N letters of the
%   alphabet, where N is the number of rows in DATA. If there are more than
%   26 rows, the labels will be AA, AB,...,ZZ. If DATA is a DataMatrix
%   object then the default labels are the row labels of DATA.
%
%   MICROPLATEPLOT(...,'COLUMNLABELS',COLUMNLABELS) allows you to specify
%   the labels for the columns. The default values are 1:M where M is the
%   number of columns in DATA. If DATA is a DataMatrix object then the
%   default labels are the columns labels of DATA.
%
%   MICROPLATEPLOT(...,'TEXTLABELS',TEXT) allows you to specify text to
%   overlay on the wells. TEXT should be a cell array of the same size as
%   DATA.
%
%   MICROPLATEPLOT(...,'TEXTFONTSIZE',FONTSIZE) allows you to specify the
%   font size to be used for TEXTLABELS. If you do not specify
%   'TEXTFONTSIZE' then an appropriate size will be chosen automatically
%   based on the size of the figure.
%
%   MICROPLATEPLOT(...,'MISSINGVALUECOLOR',COLOR) allows you to specify the
%   color to use for missing (NaN) values. Default is black.
%
%   MICROPLATEPLOT(...,'TOOLTIPFORMAT',FORMAT) allows you to specify the
%   format of the text used in the well tooltips, which display the actual
%   value from the input matrix when you click a well. FORMAT should be
%   a format string as used in SPRINTF. The default is 'Value: %.3f', which
%   specifies including three digits to the right of the decimal in fixed
%   point notation.
%
%   H = MICROPLATEPLOT(...) returns the handle to the axes of the plot.
%
%   Examples:
%
%       load microPlateAssay
%       microplateplot(assaydata);
%       colormap(whiteToRed);
%
%       % mark one cell in the assay with an 'X'
%       mask = cell(8,12);
%       mask{5,8} = 'X';
%       microplateplot(assaydata,'TEXTLABELS',mask);
%
% See also: IMAGESC, SPRINTF.

% Copyright 2008-2011 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);
% Note that we pass data to parse_inputs in order to figure out appropriate
% defaults for the labels.
[rowlabels, columnlabels,textlabels,missingvaluecolor,tooltipformat,textFontSize] = parse_inputs(data,varargin{:});

numRows = size(data,1);
numCols = size(data,2);

% set up the plate dimensions
% The sizes are taken from a spec for 96-well plate readers

A1Pos = [14.3 11.2];
wellDiameter = 6.35;
wellStep = 9;
xMax = A1Pos(1)*2+ (numCols-1)*wellStep;
yMax = A1Pos(2)*2+ (numRows-1)*wellStep;
outerDims = [0 0;xMax 0; xMax yMax; 0 yMax];

% We are doing odd things with the axes so clear the current figure.
clf;

% Create the background and label the rows and columns

hAxes = axes('xLim',[0 xMax],'ylim',[0 yMax],...
    'xtick',(A1Pos(1):wellStep:A1Pos(1)+(numCols-1)*wellStep),...
    'ytick',(A1Pos(2):wellStep:A1Pos(2)+(numRows-1)*wellStep),...
    'xticklabel',columnlabels,...
    'yticklabel',rowlabels,...
    'ydir','reverse');

% Draw a white background
patch(outerDims(:,1),outerDims(:,2),'w')
axis tight

% set the figure callbacks
hFig = get(hAxes,'parent');
% clean out any old listeners
if isappdata(hFig,'MicroPlateListeners')
    oldListener = getappdata(hFig,'MicroPlateListeners');
    for i=1:numel(oldListener)
        if isa(oldListener{i},'handle.listener')
            delete(oldListener{i})
        end
    end
    rmappdata(hFig,'MicroPlateListeners');
end
setappdata(hFig,'MicroPlateListeners',{})
set(hFig,'WindowButtonUpFcn',@localInfoOff)
hToolbar = findall(hFig,'type','uitoolbar');
h1 = findall(hToolbar,'Tag','Exploration.Rotate');
h2 = findall(hToolbar,'Tag','Exploration.Brushing');
h3 = findall(hToolbar,'Tag','DataManager.Linking');
h4 = findall(hToolbar,'Tag','Annotation.InsertLegend');
h5 = findall(hToolbar,'Tag','Plottools.PlottoolsOn');
h6 = findall(hToolbar,'Tag','Plottools.PlottoolsOff');
h7 = findall(hToolbar,'Tag','Exploration.DataCursor');

delete([h1,h2,h3,h4,h5,h6,h7]);

%Remove unnecessary items from menus
updateUIMenus(hFig)

% Convert the data values to colors
if isa(data,'bioma.data.DataMatrix')
    coldata = double(data);
    coldata = coldata(:);
else
    coldata = data(:);
end

% create the patches
for x = 0:numCols-1
    for y = 0:numRows -1
        hPatch = drawCircle(A1Pos + [x*wellStep,y*wellStep],...
            wellDiameter , coldata((y+1)+(numRows*x)),missingvaluecolor);
        setappdata(hPatch,'theValue',coldata((y+1)+(numRows*x)));
        if iscell(tooltipformat)
            setappdata(hPatch,'tooltipFormatString',tooltipformat{y+1,x+1});
        else
            setappdata(hPatch,'tooltipFormatString',tooltipformat);
        end
        % set patch callbacks
        set(hPatch,'ButtonDownFcn',@localShowInfo)
    end
end
textColor = [0 0 0];
if isempty(textFontSize)
    fs = 9;
else
    fs = textFontSize;
end
if ~isempty(textlabels)
    textHandles = zeros(numRows,numCols);
    for i = 1:numRows
        for j = 1:numCols
            textHandles(i,j) = text(A1Pos(1)+((j-1)*wellStep),A1Pos(2)+((i-1)*wellStep),.1,textlabels{i,j},...
                'color',textColor,'horizontalAlignment','center',...
                'fontsize',fs,'clipping','on','visible','off');
        end
    end
    
    set(textHandles,'visible','on');
    
    set(hAxes,'UserData',textHandles);
    if isempty(textFontSize)
        setupFigureResizeFontsizeListener(hFig,hAxes,'MicroPlateListeners');
    end
    setupXLimYlimFontsizeListener(hAxes,'MicroPlateListeners');
end

% Return that handle to the axes
if nargout > 0
    outAxes = hAxes;
end

% -----------------------------------------------------------
% Draw a circle -- this could probably be a nested function
function h = drawCircle(center, diameter,colour,missingvaluecolor)

x = center(1) + (diameter/2)*cos(0:.03:2*pi);
y = center(2) + (diameter/2)*sin(0:.03:2*pi);
if isnan(colour)
    colour = missingvaluecolor;
end
h = patch(x,y,colour);

% -----------------------------------------------------------
% Callbacks for adding datatips. There is probably a better way of doing
% this but this does the trick.
% -----------------------------------------------------------
function localInfoOff(fig,eventdata) %#ok<*INUSD>
hText = findobj(fig,'tag','WellDataTip');
if ~isempty(hText) && ishandle(hText)
    delete(hText)
end

% -----------------------------------------------------------
function localShowInfo(hSrc,eventdata)
% check whether to show info
value = getappdata(hSrc,'theValue');
formatString = getappdata(hSrc,'tooltipFormatString');
% setup axes
ax = get(hSrc,'Parent');
axUnits = get(ax,'Units');
set(ax,'Units','pixels');

% axPos = get(ax,'Position');
axXLim = get(ax,'XLim');
axYLim = get(ax,'YLim');


% point relative to data
ax_cp = get(ax,'CurrentPoint');
x = round(ax_cp(1,1));
y = round(ax_cp(1,2));
str = {sprintf(formatString,value)};

htext = text(x,y,str,'Visible','off','Clipping','off','FontName','FixedWidth');

% give it an off white background, black text and grey border
set(htext, 'BackgroundColor', [1 1 0.933333],...
    'Color', [0 0 0],...
    'EdgeColor', [0.8 0.8 0.8],...
    'Tag','WellDataTip',...
    'interpreter','none');

% determine the offset in pixels
set(htext,'position',[x y]);
set(htext,'Units','pixels')
pixpos = get(htext,'Position');

offsets = [0 0 0];

% determine what quadrant the pointer is in
quadrant=[x y]<[mean(axXLim) mean(axYLim)];

if ~quadrant(1)
    set(htext,'HorizontalAlignment','Right')
    offsets(1) = -2;
else
    set(htext,'HorizontalAlignment','Left')
    offsets(1) = 16;
end
if ~quadrant(2)
    set(htext,'VerticalAlignment','Bottom')
    offsets(2) = 2;
else
    set(htext,'VerticalAlignment','Top')
    offsets(2) = -2;
end

set(htext,'Position',pixpos + offsets);

% show the text
set(htext, 'Visible','on')

set(ax,'Units',axUnits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rowlabels, columnlabels,textlabels,missingvaluecolor,tooltipformat,textFontSize] = parse_inputs(data,varargin)
% Parse the varargin parameter/value inputs

% Check that we have the right number of inputs
if rem(nargin,2)== 0
    error(message('bioinfo:microplateplot:IncorrectNumberOfArguments', mfilename));
end

% The allowed inputs
okargs = {'rowlabels', 'columnlabels','textlabels','missingvaluecolor','tooltipformat','textfontsize'};


% Set default values
% rowLabels
numRows = size(data,1);
if isa(data,'bioma.data.DataMatrix')
    rowlabels = data.RowNames;
else
    % These are alphabetic, A-Z unless we need AA-ZZ or AAA-ZZZ,...
    
    alphaLabels = cellstr(('A':'Z')');
    if numRows <= 26
        rowlabels = alphaLabels(1:numRows);
    elseif numRows >26 && numRows <= 26^2
        rowlabels = strcat(repmat(alphaLabels,1,26)',repmat(alphaLabels,1,26));
        rowlabels = rowlabels(1:numRows);
    else
        % if we have more than 26^2 rows then we don't try to label them
        rowlabels = repmat({''},size(data));
    end
end
% columnLabels
numCols = size(data,2);
if isa(data,'bioma.data.DataMatrix')
    columnlabels = data.ColNames;
else
    
    columnlabels = cellstr(num2str((1:numCols)'));
end
% textLabels
textlabels = [];
textFontSize = [];

% missingValueColor
missingvaluecolor = [ 0 0 0];

% Tooltip format text
tooltipformat = 'Value: %.3f';
% Loop over the values
for j=1:2:nargin-1
    % Lookup the pair
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % rowlabels
            if iscellstr(pval)
                rowlabels = pval;
            else
                error(message('bioinfo:microplateplot:RowLabelsNotCellStr'));
            end
            if numel(rowlabels) ~= numRows
                error(message('bioinfo:microplateplot:RowLabelsWrongSize'));
            end
        case 2  % columnlabels
            if iscellstr(pval)
                columnlabels = pval;
            else
                error(message('bioinfo:microplateplot:ColumnLabelsNotCellStr'));
            end
            if numel(columnlabels) ~= numCols
                error(message('bioinfo:microplateplot:ColumnLabelsWrongSize'));
            end
        case 3  % textlabels
            if iscell(pval)
                mask = cellfun(@isempty,pval);
                pval(mask) = {''};
            end
            if iscellstr(pval)
                textlabels = pval;
            else
                error(message('bioinfo:microplateplot:TextLabelsNotCellStr'));
            end
            if ~isequal(size(textlabels),size(data))
                error(message('bioinfo:microplateplot:TextLabelsWrongSize'));
            end
        case 4  % missingvaluecolor
            missingvaluecolor = pval;
        case 5 % tooltipformat
            if ~ischar(pval) && ~iscellstr(pval)
                error(message('bioinfo:microplateplot:TooltipFormatNotChar'));
            end
            tooltipformat = pval;
            if ~ischar(tooltipformat) && ~isequal(size(tooltipformat),size(data))
                error(message('bioinfo:microplateplot:TooltipFormatWrongSize'));
            end
        case 6 % fontsize
            if ~isempty(pval) && (~isnumeric(pval) || ~isscalar(pval) || pval < 0 || pval ~= floor(pval))
                error(message('bioinfo:microplateplot:IncorrectTextFontSize'));
            end
            textFontSize = pval;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setupXLimYlimFontsizeListener(hAxes,listenerName)
% Sets up listeners for zooming so we can detect if we need to change the fontsize
xLimListener = addlistener(hAxes,'XLim',...
    'PostSet',@localXLimYlimFontsizeListener);
yLimListener = addlistener(hAxes,'YLim',...
    'PostSet',@localXLimYlimFontsizeListener);
hFig = get(hAxes,'Parent');
listeners = getappdata(hFig,listenerName);

setappdata(hFig,listenerName,{listeners{:} xLimListener yLimListener});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localXLimYlimFontsizeListener(hSrc,event)%#ok

hAxes = event.AffectedObject;

textHandles = get(hAxes,'UserData');
fs = getBestFontSize(hAxes,textHandles);
if fs > 0
    set(textHandles,'fontsize',fs,'visible','on');
else
    set(textHandles,'visible','off');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fs = getBestFontSize(hAxes,textHandles)
% Guess the best font size for the text in the axes

xLim = get(hAxes,'Xlim');
yLim = get(hAxes,'Ylim');
nrows = size(textHandles,1);
ncols = size(textHandles,2);
hFig = get(hAxes,'Parent');
figHeight = get(hFig,'Position')*[0 0 0 1]';
ratio = max(nrows,ncols)/max(range(yLim),range(xLim));
fs = ceil(ratio*figHeight/2);    % the gold formula
if fs < 4
    fs = 4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setupFigureResizeFontsizeListener(hFig,hAxes,listenerName)
% Set up listeners for resizing if we need to change the fontsize

PostPositionListener = addlistener(hFig,'SizeChange',...
    @(x,y)localFigureResizeFontsizeListener(x,y,handle(hAxes)));

listeners = getappdata(hFig,listenerName);

setappdata(hFig,listenerName,{listeners{:} PostPositionListener});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localFigureResizeFontsizeListener(hSrc,event,hAxes) %#ok
textHandles = get(hAxes,'UserData');
fs = getBestFontSize(hAxes,textHandles);
if fs > 0
    set(textHandles,'fontsize',fs,'visible','on');
else
    set(textHandles,'visible','off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateUIMenus(hFig)

oldSH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

%Remove unnecessary menu items from Tool menu.
%hw = findall(hFig,'Type','uimenu','Label','&Tools');
hw = findall(hFig,'Type','uimenu','Tag','figMenuTools');
hc = get(hw,'children');
h1 = findall(hw,'tag','figMenuZoomIn');
h2 = findall(hw,'tag','figMenuZoomOut');
h3 = findall(hw,'tag','figMenuPan');
h4 = findall(hw,'tag','figMenuResetView');
h5 = findall(hw,'tag','figMenuOptions');
delete(setxor(hc,[h1,h2,h3,h4,h5]));

h6 = findall(hw,'tag','figMenuOptionsDataBar');
h7 = findall(hw,'tag','figMenuOptionsDatatip');
delete([h6 h7]);

%Remove unnecessary menu items from Insert menu.
hLegend = findall(hFig,'tag','figMenuInsertLegend');
delete(hLegend);
hZlabel = findall(hFig,'tag','figMenuInsertZLabel');
delete(hZlabel);

%add separator to colorbar
hColorbar = findall(hFig,'tag','figMenuInsertColorbar');
set(hColorbar,'separator','on');

% Repair "Help" menu
bioinfoprivate.bioFigureHelpMenu(hFig, 'Microplateplot', 'microplateplot_refpage')
 
set(0,'ShowHiddenHandles',oldSH)

