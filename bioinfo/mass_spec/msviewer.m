function  varargout=msviewer(MZ,Y,varargin) 
%MSVIEWER explores a spectrum or a set of spectra with a GUI
%
%   MSVIEWER(MZ,Y) creates a GUI to display and explore the spectra in Y.
%   MZ is the mass-charge column vector with the same size as the vector of
%   ion intensity samples Y (or spectrum). Y can be a matrix with several
%   spectrograms, all sharing the same MZ scale.
%
%   MSVIEWER(...,'MARKERS',M) sets a list of marker positions to allow its
%   exploration and easy navigation. M is a column vector with MZ values.
%   Default is M=[].
%
%   MSVIEWER(...,'GROUP',G) sets the class label for every spectrogram. A
%   different color is used for every class. G can be a numeric vector or a
%   cell array of strings with the same number of elements as spectrograms
%   in Y.
%
%   MSVIEWER(...,'YLABEL',LABEL) labels the y-axis label of the spectra.
%   The default label is 'Relative Intensity'.
%
%   Example:
%
%       load sample_lo_res
%       msviewer(MZ_lo_res, Y_lo_res)
%
%   See also MSALIGN, MSBACKADJ, MSHEATMAP, MSLOWESS, MSNORM, MSPREPRODEMO,
%   MSRESAMPLE, MSSGOLAY.

%   Copyright 2003-2011 The MathWorks, Inc.



% validate required inputs

% check inputs
bioinfochecknargin(nargin,2,mfilename);

if ~isnumeric(Y) || ~isreal(Y)
    error(message('bioinfo:msviewer:IntensityNotNumericAndReal'))
end
if ~isnumeric(MZ) || ~isreal(MZ) || ~isvector(MZ)
    error(message('bioinfo:msviewer:MZNotNumericAndReal'))
end
if size(MZ,1) ~= size(Y,1)
    error(message('bioinfo:msviewer:NotEqualNumberOfSamples'))
end

numSpectrograms = size(Y,2);

% Set defaults
B=[]; % markers
classIdsProvided = false;
y_label = 'Relative Intensity';

% validate optional inputs
nvarargin = numel(varargin);
if nvarargin
    if rem(nvarargin,2)
        error(message('bioinfo:msviewer:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'markers','group','ylabel'};
    for j=1:2:nvarargin
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:msviewer:UnknownParameterName', pname))
        elseif length(k)>1
            error(message('bioinfo:msviewer:AmbiguousParameterName', pname))
        else
            switch(k)
                case 1 % markers
                    B = pval;
                    if ~isnumeric(B) || ~isreal(B) || ~isvector(B)
                        error(message('bioinfo:msviewer:MarkersNotNumericAndReal'))
                    end
                    B=B(:);
                    if max(B)>max(MZ) || min(B)<min(MZ)
                        error(message('bioinfo:msviewer:MarkersOutOfRange'))
                    end
                case 2 % class id
                    ID = pval(:);
                    classIdsProvided = true;
                    if iscell(ID)
                        if ~iscellstr(ID)
                            error(message('bioinfo:msviewer:InvalidCellForID'))
                        end
                        if ismember('',ID)
                            error(message('bioinfo:msviewer:EmptyCellForID'))
                        end
                        [ID,validGroups] = grp2idx(ID); 
                    elseif isnumeric(ID) || islogical(ID)
                        if islogical(ID)
                            ID = double(ID);
                        end
                        if ~isreal(ID) || any(isnan(ID)) || any(rem(ID,1))
                            error(message('bioinfo:msviewer:InvalidNumericForID'))
                        end
                        [ID,validGroups] = grp2idx(ID); 
                    else
                        error(message('bioinfo:msviewer:InvalidTypeForID'))
                    end
                    %check ID has the correct size
                    if numel(ID)~=numSpectrograms
                        error(message('bioinfo:msviewer:InvalidSizeForID'))
                    end
                case 3 % ylabel
                    y_label = pval;
            end % end of switch
        end
    end
end

%Spectra data
figName = 'Mass Spectra Viewer';
uData.xaxis.label = 'M/Z';
uData.yaxis.label = y_label; %'Ion Intensity (counts/sec.)';
uData.xdata = MZ;
uData.ydata = Y;
uData.numspectra = numSpectrograms;
uData.limits.xlim = [min(MZ) max(MZ)];
intensityMinMax = [min(min(Y)), max(max(Y))];
% pad by 5% so that the top and bottom of the spectra are clear.
paddingFactor = .05 * [-1 1] * diff(intensityMinMax);
uData.limits.ylim = intensityMinMax + paddingFactor;
[uData.classid, uData.validGroups] = grp2idx(ones(numSpectrograms,1));
uData.classIdsProvided = classIdsProvided;
if classIdsProvided;
    uData.classid = ID;
    uData.validGroups = validGroups;
end
[uData.iconroot, uData.iconrootMATLAB] = msvicondir;

% Markers
uData.abstol = 0.005;
uData.markers = [];
for i=1:length(B)
    markers=sort(B);
    uData.markers(i).mzvalue = markers(i);
    uData.markers(i).max_markerline=[];
    uData.markers(i).oax_markerline=[];
    uData.markers(i).markerlabel=[];
end
uData.markerlist = [];
uData.importmarkers = [];
uData.importcanceled = true;

%Patch for highlight the x or y zoom background
uData.cursorposition=[];

% Control ui components
uData.upbutton = [];
uData.downbutton = [];
uData.delbutton = [];
uData.addbutton = [];
uData.upmenu = [];
uData.downmenu = [];
uData.delmenu = [];
uData.addmenu = [];
uData.markertoolpanel = [];
uData.importmarkermenu = [];
uData.exportmarkermenu = [];
uData.togglemenugroup = [];
uData.togglebuttongroup=[];
uData.markermodebutton=[];
uData.makermodemenu=[];
uData.tempfig = [];
% uData.bottomPanel=[];

% UI component size
uData.mpw = 100; % marker panel width
uData.paxh = 35; % init panner axes height
uData.ahedge = 65; % estimated horizonatal edges around axes boxes
uData.avedge = 30; % estimated vertical edges around axes boxes
uData.minfigwidth  = 500; % keeps marker tool panel nice
uData.minfigheight = 300;
uData.lastfigposition = [];

uData.pointer = 0;  % == -1 watch, 0 arrow/drag indicators, 1 zoom
uData.zoompointer = 'glassplus';
uData.spectralines.sh = []; %no spectra line objects yet
uData.spectralines.oh = []; % no overview spectrum yet

% Create invisible figure
screensize = get(0,'screensize');
fp = get(0,'defaultfigureposition');
fw = min(fp(3)+ 130, screensize(3)-30);% figure width
fh = fp(4); % figure height
fp = [fp(1)-(fw-fp(3))/2, fp(2)+fp(4)-fh, fw, fh]; %bring figure to center of the screen

hFig = figure('Toolbar','none',...
    'Menubar','none',...
    'HandleVisibility','callback',...
    'IntegerHandle','off',...
    'NumberTitle','off',...
    'Tag','msviewer',...
    'Name', figName,...
    'Visible','off',...
    'Position', fp,...
    'Userdata', uData,...
    'DeleteFcn', @deleteFcn);

createToolbar(hFig);
createMarkerPanel(hFig);
createMenubar(hFig);
if classIdsProvided;
    createBottomPanel(hFig);
end

uData = get(hFig, 'userdata');
uData.contextmenu.cmhandle = uicontextmenu('parent',hFig);
figure(hFig)
% ====================================================================
% Create Main axes
mainaxes = axes('parent', hFig, 'units','pixels',...
    'box','on',...
    'handlevisibility','callback', ...
    'tag','mainaxes',...
    'units','pixels',....
    'xlim', uData.limits.xlim,...
    'ylim', uData.limits.ylim,...
    'buttondownfcn', '',...
    'uicontextmenu', uData.contextmenu.cmhandle);
set(get(mainaxes,'xlabel'),...
    'string',uData.xaxis.label,...
    'tag','mainaxesxlabel')
set(get(mainaxes,'ylabel'),...
    'string',uData.yaxis.label,...
    'tag','mainaxesylabel')
% create a copy that will be underneath the main axes, and
% will be used as a border during panning operations to prevent
% background erasemode from clobbering the main axes plot box.
temp = copyobj(mainaxes,hFig);
mainaxes_border = mainaxes;
mainaxes = temp;
set(mainaxes_border,'xtick',[],'ytick',[], 'visible','off',...
    'tag','mainaxes_border')

uData.mainaxes = mainaxes;
uData.mainaxes_border = mainaxes_border;

set(hFig,'userdata',uData)
%=====================================================================
% create overview axes and patch
msoverview('init', hFig)
uData = get(hFig, 'userdata');
oviewaxes = uData.overview.oviewaxes;
oviewpatch = uData.overview.oviewpatch;
%----------------------------------------------------------------

uData.overview.oviewaxes = oviewaxes;
uData.overview.oviewpatch = oviewpatch;

% plot lines
[msh, osh]=getspectraline(hFig);
for i=1:numSpectrograms
    set(msh(i),'xdata',MZ,'ydata',Y(:,i),'visible','on');
end
avg_Y = sum(Y,2)/numSpectrograms;
set(osh,'xdata',MZ,'ydata',avg_Y,...
    'visible','on');

% create Add Marker contextmenu
uData.contextmenu.addmarker = uimenu(uData.contextmenu.cmhandle,...
    'Label','Add Marker', ...
    'Tag', 'add marker cmenu',...
    'Callback',{@manageMarkerTools, 'msvmarkers', 'add', hFig});
set(hFig, 'userdata', uData);

% add markers if B is not empty
if isempty(uData.markers)
    set([uData.downbutton,uData.downmenu,...
        uData.upbutton, uData.upmenu,...
        uData.delbutton, uData.delmenu,...
        uData.exportmarkermenu], 'Enable', 'off');
else
    uData = msvmarkers([], [], 'init', hFig);
end

set(hFig, 'userdata', uData);
set(hFig,'ResizeFcn', @resizeMSViewer);
set(hFig,'windowbuttonmotionfcn', {@motion_cb, hFig})

if nargout == 1
    varargout{1} = hFig;
end

%=============== end of MSViewer ===============================
function createMenubar(hfig)
ud = get(hfig, 'userdata');
%====File menu Items
filemenu = uimenu(hfig, 'Label','&File','Tag','file menu');

ud.importmarkermenu=uimenu(filemenu,'Label','&Import Markers from Workspace...',...
    'Callback',{@importMarkersFromWS, hfig},...
    'Tag','importMarkersItem');

ud.exportmarkermenu=uimenu(filemenu,'Label','Export Markers to Workspace...',...
    'Callback',{@exportMarkersToWS,hfig},...
    'Tag','exportMarkersItem');

uimenu(filemenu,'Label','&Print to Figure',...
    'Callback',@prinToFigure,...
    'Separator','on',...
    'Tag','printToFigItem');
% print options.
uimenu(filemenu, 'Label','Pa&ge Setup...', 'Tag','pageSetupItem',...
    'Separator', 'on',...
    'Callback', {@msprintfcn, hfig, 'pagesetup'});
uimenu(filemenu, 'Label','Print Pre&view...', 'Tag','printPreviewItem',...
    'Callback', {@msprintfcn, hfig, 'preview'});
uimenu(filemenu, 'Label','&Print...', 'Tag','printItem',...
    'Callback', {@msprintfcn, hfig, 'print'});
% close
uimenu(filemenu,'Label','&Close','Tag','closeviewer',...
    'Accelerator','W',...
    'Separator','on',...
    'Callback',@closeMSViewer);

%====Tools menu Items
toolmenu = uimenu(hfig, 'Label','&Tools','Tag','tools menu');

ud.addmenu=createToolsMenuItem(hfig, toolmenu,'&Add Marker...','add marker');

ud.delmenu = createToolsMenuItem(hfig, toolmenu,'&Delete Marker','delete marker');

ud.downmenu = createToolsMenuItem(hfig, toolmenu,'Next Marker','next marker',...
    'Separator', 'on');
ud.upmenu = createToolsMenuItem(hfig, toolmenu,'Previous Marker','previous marker');

ud.makermodemenu = uimenu(toolmenu, 'Label', 'Move Marker',...
    'Tag', 'move marker',...
    'Separator','on',...
    'Callback', {@markermodeButtonChange, hfig});
set(ud.makermodemenu, 'Checked', 'on');

zoomXYItem = createToolsMenuItem(hfig, toolmenu,'Zoom XY','zoom xy');
zoomXItem  = createToolsMenuItem(hfig, toolmenu,'Zoom X','zoom x');
zoomYItem = createToolsMenuItem(hfig, toolmenu,'Zoom Y','zoom y');

createToolsMenuItem(hfig, toolmenu,'Zoom Out','zoom out');

createToolsMenuItem(hfig, toolmenu,'Reset View','reset view', ...
    'Separator','on');

%====Window menu
matlab.ui.internal.createWinMenu(hfig);

%====Help menu Items
helpmenu = uimenu(hfig, 'Label','&Help','Tag','help menu');
bioinfostandardhelp(helpmenu);

ud.togglemenugroup = [zoomXYItem, zoomXItem, zoomYItem];
set(hfig, 'userdata', ud);
%-----------------------------------------------
function menuH = createToolsMenuItem(hFig, parent,labelname, tag, varargin)
% CREATETOOLMENUITEM Creates a Tool menu Item and connects its behavior
%    to the behavior of other ui components with the same tag.
menu_tag = [tag, ' menu'];
comp = findobj(hFig,'tag',tag);
if isempty(comp)
    menuH = [];
    return;
end

menuH = uimenu(parent,'Label',labelname,'Tag', menu_tag, varargin{:});

if strcmpi(get(comp,'type'),'uipushtool') % for toolbar pushbuttons
    set(menuH,'Callback',get(comp,'ClickedCallback'));
elseif strcmpi(get(comp,'type'),'uitoggletool') % for toolbar toggles
    set(menuH,'Callback',get(comp,'onCallback'));
elseif strcmpi(get(comp,'type'),'uicontrol') %  for marker push buttons
    set(menuH,'Callback',get(comp,'Callback'));
end


%-------------------------------
function toolbar = createToolbar(hfig)
ud =get(hfig, 'userdata');
toolbar =  uitoolbar(hfig);

markerModeIcon = makeToolbarIcon(fullfile(ud.iconroot, 'tool_marker_hand.gif'));
markerModeTool = uitoggletool(toolbar,...
    'Cdata',markerModeIcon,...
    'TooltipString','Move marker',...
    'Tag', lower('move maker'),...
    'ClickedCallback',@markermodeButtonChange,...
    'OnCallback',@markermodeButtonChange,...
    'OffCallback',[]);

set(markerModeTool, 'state', 'on');
ud.markermodebutton = markerModeTool;
%===============================================================
zoomXYIcon = makeToolbarIcon(fullfile(ud.iconroot, 'tool_zoom_z.gif'));
zoomXYTool = createToolbarToggleItem(hfig, toolbar, zoomXYIcon,...
    'Zoom XY',...
    {@manageToggleTools, 'msvzoom', 'mousezoom', hfig},...
    []);

zoomXIcon = makeToolbarIcon(fullfile(ud.iconroot, 'tool_zoom_x.gif'));
zoomXTool = createToolbarToggleItem(hfig, toolbar, zoomXIcon,...
    'Zoom X',...
    {@manageToggleTools, 'msvzoom', 'zoominx', hfig},...
    []);

zoomYIcon = makeToolbarIcon(fullfile(ud.iconroot, 'tool_zoom_y.gif'));
zoomYTool = createToolbarToggleItem(hfig, toolbar, zoomYIcon,...
    'Zoom Y',...
    {@manageToggleTools, 'msvzoom', 'zoominy', hfig},...
    []);

zoomOutIcon = makeToolbarIcon(fullfile(ud.iconroot, 'tool_zoom_out.gif'));

zoomOutTool = createToolbarPushItem(toolbar, zoomOutIcon,...
    {@managePushTools, 'msvzoom', 'zoomout', hfig},...
    'Zoom out'); %#ok<NASGU>
%======================================================================
resetViewIcon = makeToolbarIcon(fullfile(ud.iconroot, 'tool_fittoview.gif'));
resetViewTool = createToolbarPushItem(toolbar,resetViewIcon,...
    {@managePushTools, 'msvzoom', 'resetview', hfig},...
    'Reset view');
set(resetViewTool,'Separator','on')
%===================================================================
helpIcon = makeToolbarIcon(fullfile(ud.iconrootMATLAB, 'helpicon.gif'));
helpTool = createToolbarPushItem(toolbar,helpIcon,...
    'helpview(fullfile(docroot,''toolbox'',''bioinfo'', ''bioinfo.map''),''msviewer_refpage'')',...
    'Help');
set(helpTool,'Separator','on')

ud.togglebuttongroup = [zoomXYTool, zoomXTool, zoomYTool];

set(hfig, 'userdata', ud);

%-------------------------------------------------------------------
function item = createToolbarToggleItem(hfig, toolbar,icon,tooltip,oncallback,offcallback)
item = uitoggletool(toolbar,...
    'Cdata',icon,...
    'TooltipString',tooltip,...
    'Tag', lower(tooltip),...
    'ClickedCallback',{@toggleButtonChange, hfig},...
    'OnCallback',oncallback,...
    'OffCallback',offcallback);

%-------------------------------------------------------------------
function item = createToolbarPushItem(toolbar,icon,callback,tooltip)

item = uipushtool(toolbar,...
    'Cdata',icon,...
    'TooltipString',tooltip,...
    'Tag',lower(tooltip),...
    'ClickedCallback',callback);

%-------------------------------------------------------------------
function createBottomPanel(hfig)
ud = get(hfig, 'userdata');
figbgcolor = get(0,'defaultFigureColor');
ud.bottompanel = uipanel('Parent',hfig,...
    'Tag','BottomPanel',...
    'Bordertype', 'etchedin',...
    'Units','pixels',...
    'backgroundcolor', figbgcolor);
ud.bottomaxes = axes('Parent', ud.bottompanel,...
    'Tag', 'bottomaxes',...
    'Position', [0,0,1,1],...
    'visible', 'off');
set(hfig, 'userdata', ud);

function createMarkerPanel(hfig)
% Create the marker tools panel.
%
figbgcolor = get(0,'defaultFigureColor');
hFigPos = getpixelposition(hfig);
ud = get(hfig, 'userdata');

ud.markertoolpanel = uipanel('Parent', hfig, ...
    'Units', 'pixels', ...
    'Position', [5, 20, ud.mpw, hFigPos(4)-20],...
    'BackgroundColor', figbgcolor,...
    'BorderType', 'none');

uicontrol('Parent', ud.markertoolpanel,...
    'Style', 'Text',...
    'Units', 'Normalized',...
    'Position', [.01 .9 1 .05],...
    'String', 'Markers (M/Z) :',...
    'BackgroundColor',figbgcolor);

ud.markerlist = uicontrol('Parent', ud.markertoolpanel,...
    'Style', 'Listbox',...
    'Units', 'Normalized',...
    'Position', [.01 .25 1 .65],...
    'BackgroundColor', 'w',...
    'Callback', {@manageMarkerTools, 'msvmarkers', 'lists', hfig});
% control buttons
upIcon = makeToolbarIcon(fullfile(ud.iconroot, 'uparrow.gif'));
ud.upbutton = createPushButtons(ud.markertoolpanel,upIcon, 'previous marker',...
    [.16 .17 .25 .05], {@manageMarkerTools, 'msvmarkers', 'previous', hfig});

downIcon = makeToolbarIcon(fullfile(ud.iconroot, 'downarrow.gif'));
ud.downbutton = createPushButtons(ud.markertoolpanel,downIcon, 'next marker',...
    [.56 .17 .25 .05], {@manageMarkerTools, 'msvmarkers', 'next', hfig});

ud.addbutton = createPushButtons(ud.markertoolpanel,[], 'Add Marker...',...
    [.07 .1 .90 .05], {@manageMarkerTools, 'msvmarkers', 'add', hfig});

ud.delbutton = createPushButtons(ud.markertoolpanel,[], 'Delete Marker',...
    [.07 .03 .90 .05], {@manageMarkerTools, 'msvmarkers', 'delete', hfig});
set(hfig, 'userdata', ud);
% ------------------------------------------
function item = createPushButtons(panel, icon, label, position, callback)
idx=strfind(label, '...');
if ~isempty(idx)
    tag = lower(label(1:idx-1));
else
    tag = lower(label);
end
item = uicontrol('Parent', panel,...
    'Style', 'pushbutton',...
    'Callback', callback,...
    'Units', 'Normalized',...
    'Position', position,...
    'Tag', tag);
if isempty(icon)
    set(item, 'String', label);
else
    set(item, 'CData', icon);
end

% ----------------------------------------------------
function resizeMSViewer(hSrc,event)  %#ok<INUSD>
hfig = gcbf;
ud = get(hfig, 'Userdata');
figPos = getpixelposition(hfig);

is_window_docked = strcmpi('docked', get(hfig,'WindowStyle'));

if ~is_window_docked
    % check and set to minimum figure size if too small
    if ~isempty(ud.lastfigposition)
        figPos = setminfigsize(figPos,ud.lastfigposition,[ud.minfigwidth ud.minfigheight]);
    end

    ud.lastfigposition = figPos;
    set(hfig, 'userdata', ud);
    setpixelposition(hfig,figPos);
end

% Make sure maiaxes and panaxes resize well as window
% shrinks/grows.
if ud.classIdsProvided
    bottomPanelPos = [1 1 figPos(3) 20];
    setpixelposition(ud.bottompanel,bottomPanelPos);
    vfactor = 3.5;
    voffset = 0;
else
    vfactor = 3;
    voffset = 20;
end

mainAxesPos = [ud.mpw + ud.ahedge, ud.paxh + (vfactor-1)*ud.avedge, figPos(3)- ud.mpw - (ud.ahedge+20), figPos(4)- ud.paxh - vfactor*ud.avedge];
oviewAxesPos = [ud.mpw + ud.ahedge,  ud.avedge-voffset, figPos(3)- ud.mpw - (ud.ahedge+20), ud.paxh];
markerToolPos = getpixelposition(ud.markertoolpanel);
markerToolPos(4) = figPos(4) - 20;

setpixelposition(ud.mainaxes,mainAxesPos);
setpixelposition(ud.overview.oviewaxes,oviewAxesPos);
setpixelposition(ud.mainaxes_border,mainAxesPos);
setpixelposition(ud.markertoolpanel,markerToolPos);

set(hSrc, 'userdata', ud);
%----------------------------------------------------------------------
function [msh, osh]=getspectraline(hfig)
% idx - the index of Y data set to be plotted
% msh - vector of line handles for mainaxes
% osh - line handle fot oviewaxes

ud = get(hfig, 'userdata');
N = ud.numspectra;
if isempty(ud.spectralines.sh) % get new lines
    msh = line(zeros(2,N),zeros(2,N),'parent',ud.mainaxes, 'visible', 'off');
    osh = line(zeros(2,1),zeros(2,1),'parent',ud.overview.oviewaxes, 'visible', 'off');
    set(osh, 'Color', [0 0 1])
    ud.spectralines.sh = msh;
    ud.spectralines.oh = osh;
    set(ud.spectralines.sh, 'uicontextmenu', ud.contextmenu.cmhandle)
    set(ud.spectralines.oh, 'buttondownfcn', @overviewdown)
    set(hfig, 'userdata', ud);
else
    msh = ud.spectralines.sh;
    osh = ud.spectralines.oh;
    set(osh, 'tag', 'overviewline');
end

[M, I] = unique(ud.classid);
M=length(M);
% color the spectra
cmap=hsv(M);
cmap(1,:)=[ 0  0 1];
cmap(2,:)=[ 1  0 0];
cmap(3,:)=[ 0  0.8 0];
% set line colors
for i = 1:N
    set(msh(i), 'color', cmap(ud.classid(i),:))
end

% set group legend
if ud.classIdsProvided
    %                 lh =[];
    %                 for j = 1:M
    %                     lh(end+1) = msh(I(j));
    %                 end
    lh =zeros(M, 1);
    for j = 1:M
        lh(j) = msh(I(j));
    end
    hleg = legend(ud.bottomaxes,lh, ud.validGroups, 'Orientation', 'horizontal');
    tmpPos = get(hleg, 'position');
    hlegPos = [0.5 - tmpPos(3)/2, 0.5 - tmpPos(4)/2, tmpPos(3), tmpPos(4)];
    set(hleg, 'Position', hlegPos, 'box', 'off')
end
%-------------------------------------------------------------------------
function managePushTools(hSrc,event, fun, action, hfig)  %#ok<*INUSL>
msvswitch(fun, [], [], action, hfig)
%----------------------------------
function manageToggleTools(hSrc,event, fun, action, varargin)
% Manage mutually exclusive itmes on toolbar and menus
if nargin > 4
    hfig = varargin{1};
else
    hfig=gcbf;
end

ud = get(hfig, 'userdata');
type = get(hSrc,'type');

if strcmpi(type,'uimenu')
    % caller is the menu item
    tg = strrep(get(hSrc,'tag'),' menu','');
    hSrcToolbarPeer = findall(ud.togglebuttongroup,'tag', tg);

    if strcmp(get(hSrcToolbarPeer,'State'),'on')
        set(hSrcToolbarPeer,'State','off');
    else
        set(hSrcToolbarPeer,'State','on');
        msvswitch(fun, [], [], action, hfig)
    end
    toggleButtonChange(hSrcToolbarPeer,event, hfig);

else
    % caller is the toolbar button
    msvswitch(fun,[], [], action, hfig)
end

function markermodeButtonChange(hSrc, event, varargin)
if nargin > 2
    hfig = varargin{1};
else
    hfig=gcbf;
end
ud=get(hfig, 'userdata');
hZoomTool = findobj(ud.togglebuttongroup, 'State', 'on');
set(hfig,'windowbuttonmotionfcn', {@motion_cb, hfig})
if ~isempty(hZoomTool)
    set(hZoomTool, 'state', 'off');
    toggleButtonChange(hZoomTool,[], hfig)
    ud=get(hfig, 'userdata');
end
set(ud.markermodebutton, 'state', 'on');
set(ud.makermodemenu, 'Checked', 'on');
set(hfig, 'userdata', ud)

%------------------------------------------
function toggleButtonChange(hSrc,event, varargin) 
if nargin > 2
    hfig = varargin{1};
else
    hfig=gcbf;
end
ud=get(hfig, 'userdata');
menu_tag = [get(hSrc,'tag'),' menu'];
hSrcMenuPeer = findobj(ud.togglemenugroup,'tag',menu_tag);
if ~isempty(hSrcMenuPeer)
    if strcmp(get(hSrc,'State'),'on')
        otherMenu = ud.togglemenugroup(ud.togglemenugroup~=hSrcMenuPeer);
        set(ud.togglebuttongroup(ud.togglebuttongroup~=hSrc),'State','off')
        set(otherMenu,'Checked','off');
        set(hSrcMenuPeer,'Checked','on');
        set(ud.markermodebutton, 'state', 'off');
        set(ud.makermodemenu, 'Checked', 'off');
    else
        set(ud.togglemenugroup,'Checked','off');
        set(ud.markermodebutton, 'state', 'on');
        set(ud.makermodemenu, 'Checked', 'on');
    end
end
if strcmp(get(hSrc,'State'),'off')
    msvswitch('msvzoom', [], [], 'leavezoom', hfig)
    ud=get(hfig, 'userdata');
    set(ud.mainaxes, 'ButtonDownFcn', [])
end
set(hfig, 'userdata', ud)
%----------------------------------
function manageMarkerTools(hSrc, event, fun, action, varargin) 
% Manage marker tools menus, buttons, and contextmenu
if nargin > 4
    hfig = varargin{1};
else
    hfig = gcbf;
end
tag = lower(get(hSrc, 'tag'));
ud = get(hfig, 'userdata');
midx = 0;

if strcmp(tag,'add marker') || strcmp(tag,'add marker menu')
    % caller is the menu item or push button
    [assignmarker_success, midx, ud] = getMarkerInput(ud);
    set(hfig, 'userdata', ud)
elseif strcmpi(tag, 'add marker cmenu') % context menu
    p = get(ud.mainaxes, 'CurrentPoint');
    p = p(1,1:2);
    xlim = get(ud.mainaxes, 'Xlim');
    if xlim(1) > p(1)
        return;
    else
        ud.markers(end+1).mzvalue = p(1);
    end
    markers = zeros(size(ud.markers)); 
    for i =1:length(ud.markers)
        markers(i)= ud.markers(i).mzvalue;
    end
    set(hfig, 'userdata', ud)
    assignmarker_success= 1;
else
    assignmarker_success= 2;
end

if assignmarker_success ==1
    if midx ~= 0
        action='focus';
        msvswitch(fun, [], [], action, hfig, midx)
    else
        msvswitch(fun, [], [], action, hfig)
    end
elseif assignmarker_success == 2
    msvswitch(fun, [], [], action, hfig)
end

%----------------------------------------
function [assignmarker_success, idx, ud] = getMarkerInput(ud)
% Display the input dialog using INPUTDLG from MATLAB

assignmarker_success = 0;
idx=0;
dlg_title = 'Add Marker';
prompt = sprintf('(Hint: You can also add markers by right-clicking on the axis.)\n\nEnter M/Z x-axis value:\n');

num_of_lines_per_prompt = 1;

while assignmarker_success ~= 1
    markervalue = inputdlg(prompt, dlg_title, num_of_lines_per_prompt);
    if isempty(markervalue)
        return;
    end
    markervalue = str2double(markervalue{:});
    if isempty(markervalue) || isnan(markervalue) || ~isreal(markervalue)
        err_str = sprintf('Marker M/Z value must be numeric and real');
        assignmarker_success = 0;
        uiwait(errordlg(err_str, dlg_title));
    elseif markervalue > max(ud.xdata) || markervalue < min(ud.xdata)
        err_str = sprintf('Marker M/Z value must be within the MZ scale: %0.2f - %0.2f',...
            ud.limits.xlim(1), ud.limits.xlim(2));
        assignmarker_success = 0;
        uiwait(errordlg(err_str, dlg_title));
    else
        x = zeros(size(ud.markers));
        for i =1:length(ud.markers)
            x(i)= ud.markers(i).mzvalue;
        end
        existMZ = find(abs(x-markervalue) < ud.abstol);
        if isempty(existMZ)
            ud.markers(end+1).mzvalue = markervalue;
        else
            idx=existMZ(1);
        end

        assignmarker_success = 1;
    end
end

%=================== Callback functions ==============================
function motion_cb(hSrc, event, hfig) 
%         hfig=gcbf;
ud = get(hfig,'userdata');
flag = ud.pointer;

switch flag
    case -1,  % wait mode
        setptr(hfig,'watch'),
    case 0,  % pointer mode
        marker_curs = msvmarkers([], [], 'motion', hfig);
        if marker_curs == 1
            setptr(hfig,'hand')
            return
        end

        oview_curs = msoverview('motion',hfig);
        if oview_curs == 1
            setptr(hfig,'hand')
            return
        end

        set(hfig,'pointer','arrow')

    case 1,  % zoom mode
        p = get(hfig,'currentpoint');  p = p(1,1:2);
        mp = getpixelposition(ud.mainaxes);

        %mouse is in main panel:
        if pointsinrect(p(1,1:2),[mp(1) mp(1)+mp(3) mp(2) mp(2)+mp(4)])
            setptr(hfig,ud.zoompointer)
            return
        end

        oview_curs = msoverview('motion',hfig);
        if oview_curs == 1
            setptr(hfig,'hand')
            return
        end

        setptr(hfig,'arrow')
end
%------------------------------------------------
function prinToFigure(hSrc, event)  %#ok<INUSD>
hfig = gcbf;
msprinttofigure(hfig)

%-----------------------------------------------
function closeMSViewer(obj,evt) %#ok<INUSD>
hfig = gcbf;
close(hfig);

%---------------------------------------------
function deleteFcn(obj,event) %#ok<INUSD>
hfig = gcbf;
if ~isempty(hfig) && ishandle(hfig)
    ud = get(hfig,'userdata');
    if ~isempty(ud.tempfig) && ishandle(ud.tempfig)
        delete(ud.tempfig)
    end
    if ~isempty(ud.spectralines)
        delete(ud.spectralines.sh)
        delete(ud.spectralines.oh)
    end

    delete(hfig);
end

%============== Marker functions =============================
function varargout = msvmarkers(varargin)
%Marker management function.

if nargin < 3
    action = 'add';
else
    action = varargin{3};
end

fig=varargin{4};
if strcmpi(action, 'init') == 1
    action = 'add';
    %     fig=varargin{4};
    ud = get(fig,'userdata');
else
    %     fig = gcbf;
    ud = get(fig,'userdata');
    [idx, mzx] = findselectedmarker(ud);  %#ok<NASGU>
end

switch lower(action)
    %-----------------------------------------------------------------
    % idx = msvmarkers('motion',fig)
    % returns marker index if currentpoint is over markers(index), 0 else
    case 'motion'
        if isempty(idx)
            idx = 0;
            %         elseif length(idx) > 0
        else
            idx=1;
        end

        varargout{1} = idx;
        % add marker to the mainaxes and oviewaxes
    case 'add'
        markerline_props = {
            'visible','on',...
            'linewidth', 1, ...
            'color',[1 0 1],... % m
            'LineStyle', '--',...
            'handlevisibility','callback'};
        markerlabel_props = {
            'parent', ud.mainaxes,...
            'visible','off',...
            'linewidth', 1, ...
            'Rotation', 0,...
            'Editing', 'off',...
            'FontSize', 8,...
            'handlevisibility','callback'};
        marker_count=zeros(length(ud.markers));
        for i=1:length(ud.markers)
            if ~isfield(ud.markers(i), 'max_markerline')
                ud.markers(i).max_markerline = [];
                ud.markers(i).oax_markerline = [];
                ud.markers(i).markerlabel = [];
            end

            if isempty(ud.markers(i).max_markerline)
                ud.markers(i).cmenu = uicontextmenu('parent', fig);
                xd = [ud.markers(i).mzvalue, ud.markers(i).mzvalue];
                yd = ud.limits.ylim;
                hml = line(xd, yd,...
                    'parent', ud.mainaxes,...
                    markerline_props{:},...
                    'UIContextMenu', ud.markers(i).cmenu,...
                    'buttondownfcn',{@markerdown, fig});
                uimenu(ud.markers(i).cmenu, 'Label', 'Delete Marker',...
                    'Callback', {@msvmarkers, 'delete', fig});

                yd_label = get(ud.mainaxes, 'ylim');

                hmlabel = text( double(xd(1)), double(yd_label(2)),num2str(ud.markers(i).mzvalue,'%0.2f'),...
                    markerlabel_props{:},...
                    'buttondownfcn',{@markerdown, fig});

                yd = get(ud.overview.oviewaxes, 'ylim');
                hol = line(xd, yd,...
                    'parent', ud.overview.oviewaxes,...
                    'buttondownfcn', @overviewdown,...
                    markerline_props{:});

                ud.markers(i).max_markerline = hml;
                ud.markers(i).oax_markerline = hol;
                ud.markers(i).markerlabel = hmlabel;
                marker_count(i)=i;
            end
        end
        idx = marker_count == 0;
        marker_count(idx) = [];

        if ~isempty(marker_count)
            ud = focusonmarker(ud, marker_count(1)); % set the first one to focus
            ud = managemarkerlist(fig, ud, marker_count(1));
        end

        set(fig,'userdata',ud)
        if nargout ==1
            varargout{1} = ud;
        end
    case 'delete'
        if isempty(idx)
            val = get(ud.markerlist, 'Value');
            marker_idx=findmarkeridex(ud, val);
        else
            marker_idx=idx(1);
        end

        delete([ud.markers(marker_idx).max_markerline,...
            ud.markers(marker_idx).oax_markerline,...
            ud.markers(marker_idx).markerlabel])
        ud.markers(marker_idx)=[];
        ud = managemarkerlist(fig, ud);
        set(fig, 'userdata', ud)
    case 'next'
        ud = managemarkerlist(fig, ud,'down');
        set(fig, 'userdata', ud)
    case 'previous'
        ud = managemarkerlist(fig, ud,'up');
        set(fig, 'userdata', ud)
    case 'lists'
        val = get(ud.markerlist, 'Value');
        if ~isempty(val)
            marker_idx=findmarkeridex(ud,val);
            ud = centermarker(fig, ud, marker_idx);
            ud=focusonmarker(ud, marker_idx);
            ud = managemarkerlist(fig, ud, marker_idx);

            set(fig, 'userdata', ud)
        end
    case 'focus'
        marker_idx=varargin{5};
        ud=focusonmarker(ud, marker_idx);
        ud = managemarkerlist(fig, ud, marker_idx);
        set(fig, 'userdata', ud)
end

%---------------------------------------------------------
function [sx, orderIds] = getindexedmarkerlist(ud)
x = zeros(size(ud.markers));
for i =1:length(ud.markers)
    x(i)= ud.markers(i).mzvalue;
end
[sx, orderIds]=sort(x);

function marker_idx = findmarkeridex(ud, listidx)
marker_idx = 0;
[sx, orderIds] = getindexedmarkerlist(ud);
if ~isempty(sx)
    marker_idx=orderIds(listidx);
end

function markerdown(hSrc, event, fig) 
%MSVMarkerdown Button down function for ruler lines.
ud = get(fig,'userdata');
active_zoomobj = findobj(ud.togglebuttongroup, 'State', 'on');
if isempty(active_zoomobj)
    marker_idx=0;
    num_found = findselectedmarker(ud);
    if ~isempty(num_found)
        marker_idx=num_found(1);
    end
    %xxx         if strcmp(get(gcf,'pointer'),'arrow')
    if strcmp(get(fig,'pointer'),'arrow')
        % not a hand - just return
        return
    end

    save_wbmf = get(fig,'windowbuttonmotionfcn');
    setptr(fig,'closedhand')
    set(ud.markers(1).max_markerline,'userdata',0)

    ud=focusonmarker(ud, marker_idx);
    ud = managemarkerlist(fig, ud, marker_idx, 'markerdown');
    set(fig,'windowbuttonmotionfcn',{@movemarker, ud, marker_idx})
    set(fig,'windowbuttonupfcn',['msvud = get(gcf,''userdata''); '...
        'set(msvud.markers(1).max_markerline,''userdata'',1), clear msvud'])
    waitfor(ud.markers(1).max_markerline,'userdata',1)

    % OK we're back from the drag
    setptr(fig, 'hand')
    set(fig,'windowbuttonmotionfcn',save_wbmf);
    set(fig,'windowbuttonupfcn','');
    %         set(fig, 'userdata', ud);
end

%--------------------------------------------------------------------------
function ud=focusonmarker(ud, marker_idx)
for i=1:length(ud.markers)
    if i==marker_idx
        set(ud.markers(i).max_markerline, 'LineStyle', '-', 'LineWidth', 1.5);
        set(ud.markers(i).oax_markerline, 'LineStyle', '-', 'LineWidth', 1.5);
        xdata=get(ud.markers(i).max_markerline, 'Xdata');
        xlim = get(ud.mainaxes,'xlim');
        ylim = get(ud.mainaxes,'ylim');
        if ud.markers(i).mzvalue > xlim(1) && ud.markers(i).mzvalue < xlim(2)
            set(ud.markers(i).markerlabel, 'visible', 'on',...
                'Edgecolor', 'none',...
                'BackgroundColor', [0.8 0.8 0.6])
            positionmarkerlabel(ud.markers(i).markerlabel, xdata(1), ylim(2))
        end
    else
        set(ud.markers(i).max_markerline, 'LineStyle', '--', 'LineWidth', 1);
        set(ud.markers(i).oax_markerline, 'LineStyle', '--', 'LineWidth', 1);
        set(ud.markers(i).markerlabel, 'visible', 'off')
    end
end

%--------------------------------------------------------------
function movemarker(hSrc, event, ud, marker_idx) 
% Windowbuttonmotionfcn for dragging a marker line

p = get(ud.mainaxes,'currentpoint');
p = p(1,1:2);
if isnan(p)  % case in which p is outside of axes
    return
end

x=p(1);
xlim = get(ud.mainaxes,'xlim');
ylim =get(ud.mainaxes, 'ylim');
if x<xlim(1)
    x = xlim(1);
elseif x>xlim(2)
    x = xlim(2);
end

% now set the marker values
set(ud.markers(marker_idx).max_markerline, 'XData', [x, x]);
set(ud.markers(marker_idx).oax_markerline, 'XData', [x, x]);
set(ud.markers(marker_idx).markerlabel,'String', num2str(x, '%0.2f'));
positionmarkerlabel (ud.markers(marker_idx).markerlabel, x, ylim(2))
ud.markers(marker_idx).mzvalue=x;

% handle the list
ud = managemarkerlist(hSrc, ud, marker_idx, 'moving');
set(hSrc, 'userdata', ud)
%-------------------------------------------------------------------
function ud = centermarker(hfig, ud, marker_idx)  
xlim = get(ud.mainaxes, 'xlim');
ylim = get(ud.mainaxes, 'ylim');
if (xlim == ud.limits.xlim)
    return
end
delta = diff(xlim)/2;
x1=ud.markers(marker_idx).mzvalue - delta;
x2=(ud.markers(marker_idx).mzvalue + delta);
xlim = [x1,x2];
xlim = fitwithinXlimit(ud, xlim);
set(ud.mainaxes, 'YLim', ylim, 'XLim', xlim);
msoverview('zoom', xlim, ylim, hfig)
%---------------------------------------------
function [num_found, varargout] = findselectedmarker(ud)
p = get(ud.mainaxes,'currentpoint');

x = zeros(size(ud.markers));
for i =1:length(ud.markers)
    x(i)= ud.markers(i).mzvalue;
end
% 3.5 is the "magic number" (used to be 6.5)
mpos = get(ud.mainaxes,'position');
five_pix = 3.5*diff(get(ud.mainaxes,'xlim'))/mpos(3);

num_found = find(abs(p(1)-x)<=five_pix);

% Output a vector of marker mz values
if nargout>1
    varargout{1}=x;
end

function xlim = fitwithinXlimit(ud, xlim)
if xlim(1) < ud.limits.xlim(1)
    xlim(1) = ud.limits.xlim(1);
end
if xlim(2) > ud.limits.xlim(2)
    xlim(2) = ud.limits.xlim(2);
end


%------------------------------------------------------
function positionmarkerlabel (hLabel, x, y)
l_extent=get(hLabel, 'Extent');
xlabel = x-l_extent(3)/2;
if xlabel <=0
    xlabel = x;
end
ylabel=y+l_extent(4)/1.5;
set(hLabel, 'position', [xlabel,ylabel])


%---------------------------------------------
function [ud, varargout] = managemarkerlist(hfig, ud, varargin)
mode_type = 'leavezoom';
is_leaving_zoom = false;
is_moving = false;

if nargin == 3 && ischar(varargin{1}) && strcmpi(varargin{1}, 'inzoom')
    mode_type = varargin{1};
end

if nargin == 4
    if islogical(varargin{2})
        is_leaving_zoom = varargin{2};
    else
        is_moving = true;
    end
end

switch mode_type
    case 'inzoom'
        set([ud.downbutton,ud.downmenu,...
            ud.upbutton, ud.upmenu,...
            ud.delbutton, ud.delmenu,...
            ud.addbutton, ud.addmenu,...
            ud.importmarkermenu,...
            ud.exportmarkermenu,...
            ud.markerlist], 'Enable', 'off');
    case 'leavezoom'
        if ~isempty(ud.markers)
            [sx, orderIds]=getindexedmarkerlist(ud);
            val = get(ud.markerlist, 'Value');

            if nargin == 2 % delete
                if val == length(sx)
                    val = 1;
                elseif val > length(sx)
                    val = length(sx);
                end

                marker_idx = orderIds(val);
                if length(sx) ==1
                    marker_idx=1;
                end
                if ~is_leaving_zoom
                    ud = centermarker(hfig, ud, marker_idx);
                end
                ud=focusonmarker(ud, marker_idx);
            elseif nargin == 3
                if isnumeric(varargin{1}) % add marker
                    marker_idx = varargin{1};
                    val = find(ud.markers(marker_idx).mzvalue == sx);
                    if length(val) > 1
                        val=val(1);
                    end
                    ud = centermarker(hfig, ud, marker_idx);
                    ud=focusonmarker(ud, marker_idx);
                else % up & down
                    type = varargin{1};
                    if strcmpi(type, 'down') == 1 % down
                        val = val+1;
                    elseif strcmpi(type, 'up') == 1 % up
                        val = val-1;
                    end
                    marker_idx = orderIds(val);
                    ud = centermarker(hfig, ud, marker_idx);
                    ud=focusonmarker(ud, marker_idx);
                end
            elseif nargin ==4
                if is_moving % add marker
                    marker_idx = varargin{1};
                    val = find(ud.markers(marker_idx).mzvalue == sx);
                    if length(val) > 1
                        val=val(1);
                    end
                end
            end

            set(ud.markerlist, 'String', sx);
            set(ud.markerlist, 'Value', val);
            % Disable buttons and menus
            mlist_str=cellstr(get(ud.markerlist, 'String'));
            set([ud.delbutton, ud.delmenu, ud.exportmarkermenu], 'Enable', 'on');
            if length(mlist_str) == 1
                set([ud.downbutton,ud.downmenu,...
                    ud.upbutton, ud.upmenu], 'Enable', 'off');
            else
                if val==length(mlist_str)
                    set([ud.downbutton,ud.downmenu], 'Enable', 'off');
                    set([ud.upbutton,ud.upmenu], 'Enable', 'on');

                elseif val==1
                    set([ud.upbutton, ud.upmenu], 'Enable', 'off');
                    set([ud.downbutton,ud.downmenu], 'Enable', 'on');
                else
                    set([ud.downbutton,ud.downmenu,...
                        ud.upbutton, ud.upmenu,], 'Enable', 'on');
                end

            end
        else
            sx = [];
            set(ud.markerlist, 'String', sx);
            set([ud.downbutton,ud.downmenu,...
                ud.upbutton, ud.upmenu,...
                ud.delbutton, ud.delmenu,...
                ud.exportmarkermenu], 'Enable', 'off');
        end
        set([ud.addbutton,ud.addmenu,...
            ud.markerlist, ud.importmarkermenu], 'Enable', 'on');
end

% Output a vector of marker mz values
if nargout>1
    varargout{1}=sx;
end

%=================== Overviews Patchs ====================================
function varargout = msoverview(varargin)
%msoverview  overview tool for ms spcetra in display.

if nargin < 1
    action = 'init';
else
    action = varargin{1};
end

switch lower(action)
    %-----------------------------------------------------------------
    % cursoronpatch = msoverview('motion',fig)
    %  returns 1 if cursor is over patch, 0 else
    %
    case 'motion'
        fig = varargin{2};
        ud = get(fig,'userdata');
        oviewaxes = ud.overview.oviewaxes;
        p = get(oviewaxes,'currentpoint');
        oviewpatch = ud.overview.oviewpatch;
        xd = get(oviewpatch,'xdata');
        yd = get(oviewpatch,'ydata');
        varargout{1} = pointsinrect(p(1,1:2),[xd([1 2]) yd([1 3])]);

        %-----------------------------------------------------------------
        % overview  -or-  overview('init',fig)
        %  startup code - adds overview to fig
        %
    case 'init'
        fig = varargin{2};
        ud = get(fig,'userdata');

        % create overview axes and patch
        oviewaxes = axes( ...
            'parent',fig,...
            'tag', 'overviewaxes', ...
            'box', 'off', ...
            'units','pixels',...
            'xlim', ud.limits.xlim,...
            'ylim', ud.limits.ylim,...
            'xtick', [], ...
            'ytick', [] );

        edgecolor = get(oviewaxes,'xcolor');

        pc = get(oviewaxes,'color');
        if ~ischar(pc)  % might be 'none', in which case don't set x,ycolor of axes
            set(oviewaxes,'xcolor',pc)
            set(oviewaxes,'ycolor',pc)
        else
            fc = get(fig,'color');
            set(oviewaxes,'xcolor',fc)
            set(oviewaxes,'ycolor',fc)
        end

        ud.overview.oviewaxes = oviewaxes;

        xlim = get(ud.mainaxes,'xlim');
        ylim = get(ud.mainaxes,'ylim');

        xd = xlim([1 2 2 1 1]);
        yd = ylim([1 1 2 2 1]);
        oviewpatch = patch( ...
            xd, ...
            yd, ...
            [1 1 1],...
            'facecolor',pc, ... 
            'parent',oviewaxes,...
            'tag','overviewpatch',...
            'edgecolor',edgecolor,...
            'buttondownfcn', @overviewdown);

        set(oviewpatch,'facecolor','none')
        ud.overview.oviewpatch = oviewpatch;

        set(fig,'userdata',ud)

        % overview('close',fig)
        %  shutdown code - removes overview from msviewer
        %   Inputs:
        %     fig - figure handle of the viewer
        %
    case 'close'
        fig = varargin{2};
        delete(findobj(fig,'tag','oviewaxes')) % delete oviewaxes (and all children)
        set(fig,'userdata',ud) %#ok<NODEF>

        % msoverview('update',fig)
        %   assume xlimits of panaxes are correctly set to full view
        %   reset patch to limits of mainaxes
        %
    case 'update'
        fig = varargin{2};

        ud = get(fig,'userdata');
        xlim = get(ud.mainaxes,'xlim');
        ylim = get(ud.mainaxes,'ylim');
        setpdata(ud.overview.oviewpatch,xlim,ylim)

        % msoverview('zoom', xlim,ylim)
        %   set patch limits based on limits input
    case 'zoom'
        fig = varargin{4};
        xlim = varargin{2};
        ylim = varargin{3};

        ud = get(fig,'userdata');
        setpdata(ud.overview.oviewpatch,xlim,ylim)

end

%--------------------------------------------------------------------
function overviewdown(varargin) 
%OVERVIEWDOWN Button down function for overview lines and patch.

fig = gcbf;
ud = get(fig,'userdata');

if ~strcmp(get(fig,'pointer'),'custom') || ...
        ~((ud.pointer==0) || (ud.pointer==1))
    % we're not actually on the patch, or we're in the wrong
    % mouse mode (only pan if in pointer or zoom modes).
    return
end

if ~isempty(ud.spectralines)
    %     marker_idx=[];
    for i=1:length(ud.markers)
        if strcmpi(get(ud.markers(i).markerlabel, 'visible'), 'on')
            set(ud.markers(i).markerlabel, 'visible', 'off');
            %             marker_idx=i;
            break;
        end
    end

    if msvpan('Ax',ud.mainaxes,...
            'Bounds',ud.limits,...
            'BorderAxes',ud.mainaxes_border,...
            'DirectFlag',0,...
            'OviewPatch',ud.overview.oviewpatch,...
            'DynamicDrag',1)
        ud=get(fig, 'userdata');
        xlim = get(ud.mainaxes,'xlim');
        ylim = get(ud.mainaxes,'ylim');
        idx = get(ud.markerlist,'Value');
        focus_marker=findmarkeridex(ud, idx);
        for i=1:length(ud.markers)
            if ud.markers(i).mzvalue > xlim(1) && ud.markers(i).mzvalue < xlim(2)
                if i==focus_marker
                    set(ud.markers(i).markerlabel, 'visible', 'on');
                    xdata=get(ud.markers(i).max_markerline, 'XData');
                    positionmarkerlabel(ud.markers(i).markerlabel, xdata(2), ylim(2))
                end
            else
                set(ud.markers(i).markerlabel, 'visible', 'off');
            end
        end

        set(fig, 'userdata', ud)
    end
end

%=================== Pan functions =======================================
function varargout = msvpan(varargin)
%MSVPAN  Pan Function for MSViewer Axes
%  Usage: call this function either as the buttondownfcn of the lines in
%   an axes that you want to pan, or as the windowbuttondownfcn when you
%   click inside an axes.

% OUTPUTS:
%   returns 1 if limits changed, 0 else
%   Copy from Sptoolui\panfcn.m

ax = [];
bounds = [];
directflag = 1;
borderaxes = [];
oviewpatch = [];
oviewaxes = [];  % parent of overviewpatch
dynamicdrag = 1;
data = [];
transform = '';
immediate = 1;
invisible = [];
userhand = [];
pointer = 'closedhand';
interimpointer = '';

% parse parameters - could use some error checking
for i=1:2:length(varargin)
    param = lower(varargin{i});
    value = varargin{i+1};  %#ok<NASGU>
    eval([param ' = value;'])
end

% finish default value setting
if isempty(ax)
    ax = gca;
end
fig = get(ax,'parent');
ptr = getptr(fig);
if ~isempty(interimpointer)
    setptr(fig,interimpointer)
end
if isempty(userhand)
    userhand = ax;
end
if ~isempty(oviewpatch)
    oviewaxes = get(oviewpatch,'parent');
end
% initialize "data" structure array from axes line objects here.
% This costly computation will be skipped if you pass in the "data"
% structure array yourself.
% This array has fields:
%    .data - matrix
%    .columns - index vector; indicates which columns of .data are
%        actually displayed as lines
%    .h - line handles - one per column of .data(:, data.columns)
%    .xdata - either the same size as .data or empty if .Fs >0
%    .Fs - scalar; -1 indicates xdata matrix contains x information,
%         if >0, xdata is created evenly spaced
%    .t0 - scalar; used with .Fs to compute xdata
if isempty(data)
    h = findobj(ax,'type','line');
    h = setdiff(h,invisible);
    if isempty(h)
        setptr(fig,ptr)
        return
    end
    dataStruct = struct('h',[],'data',[],'columns',0,'xdata',[],'t0',[],'Fs',[]);
    data = repmat(dataStruct,length(h),1);
    for i=1:length(h)
        data(i).h = h(i);
        data(i).data = get(h(i),'ydata');
        data(i).data = data(i).data(:);  % make into column
        data(i).columns = 1;
        data(i).xdata = get(h(i),'xdata');
        data(i).xdata = data(i).xdata(:);  % make into column
        xdf=diff(data(i).xdata);   % see if this line is evenly spaced
        if length(xdf)>1  &&  all(xdf==xdf(1)) &&  xdf(1)~=0
            data(i).t0 = data(i).xdata(1);  data(i).Fs = 1./xdf(1);
        else
            data(i).t0 = 0; data(i).Fs = -1;
        end
    end
end

save_userhand_tag = get(userhand,'tag');
set(userhand,'tag','userhand')

% get current point of mouse click
if directflag
    p = get(ax,'currentpoint');
else
    p = get(oviewaxes,'currentpoint'); %#ok<UNRCH>
end
p = p(1,1:2);
np = p;

if immediate
    [oldptr,save_visible] = start_motion(fig,ax,borderaxes,invisible,pointer);
end

xlim = get(ax,'xlim');
ylim = get(ax,'ylim');

% To get fast panning without axes redrawing, we never change
% the xlim and ylim of ax during the dragging.  Instead, we
% set the xdata and ydata of the lines during dragging and
% set the xlim and ylim when we are done.

actual_xlim = xlim;
actual_ylim = ylim;
logscale = [strcmp(get(ax,'xscale'),'log') strcmp(get(ax,'yscale'),'log')];

save_callBacks = ...
    installCallbacks(userhand,fig,...
    {'windowbuttonmotionfcn', 'windowbuttonupfcn'},{'motion', 'up'});

done = 0;
moved = 0;

while ~done
    event = waitForNextEvent(userhand);
    switch event
        case 'motion'
            if ~moved,
                if ~immediate
                    [oldptr,save_visible] = start_motion(fig,ax,borderaxes,...
                        invisible,pointer);
                end
                moved = 1;
            end

            p = np;
            [xlim,ylim,oxlim,oylim,np] = ...
                draglims(directflag,ax,oviewaxes,xlim,ylim,bounds,p,logscale);

            if ~isequal([xlim ylim],[oxlim oylim])
                if ~isempty(oviewpatch)
                    setpdata(oviewpatch,xlim,ylim)
                end
                if dynamicdrag
                    doDynamicDrag(xlim,ylim,actual_xlim,actual_ylim,....
                        data,transform)
                end
            end

        case 'up'
            done = 1;
    end
end

if immediate || moved
    p = np;
    [xlim,ylim,oxlim,oylim,np] = ...
        draglims(directflag,ax,oviewaxes,xlim,ylim,bounds,p,logscale);  %#ok<NASGU>

    stop_motion(ax,borderaxes,invisible,save_visible)
    if dynamicdrag && moved
        % need to restore x and y data
        for i=1:length(data)
            if data(i).Fs > 0
                xdata = data(i).t0 + (0:size(data(i).data,1)-1)/data(i).Fs;
            else
                xdata = data(i).xdata;
            end
            for j=1:length(data(i).columns)
                y = data(i).data(:,data(i).columns(j));
                if ~isempty(transform)
                    y = feval(transform,y);
                end
                set(data(i).h(j),'xdata',xdata,'ydata',y)
            end
        end
    end

    set(ax,'xlim',xlim,'ylim',ylim)
    set(fig,oldptr{:})
end

set(fig,{'windowbuttonmotionfcn' 'windowbuttonupfcn'},save_callBacks)
set(userhand,'tag',save_userhand_tag)

if nargout>0
    varargout{1} = moved;
end
set(fig,ptr{:})
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [oldptr,save_visible] = start_motion(fig,ax,borderaxes,...
                                            invisible,pointer)
oldptr = getptr(fig);
setptr(fig,pointer)
switchToBorderAxes(ax,borderaxes)

save_visible = get(invisible,{'visible'});
set(invisible,'visible','off')

function stop_motion(ax,borderaxes,invisible,save_visible)
set(invisible,{'visible'},save_visible)
switchToMainAxes(ax,borderaxes)

function [xlim,ylim,oxlim,oylim,np] = draglims(directflag,ax,...
    oviewaxes,xlim,ylim,bounds,p,logscale)
% This function takes the current point information and the old limits,
% and calculates the new xlims and ylims.
if directflag
    np = get(ax,'currentpoint');  % new point
else
    np = get(oviewaxes,'currentpoint');  % new point
end
np = np(1,1:2);

if any(isnan(np)), np = p; end

oxlim = xlim;
oylim = ylim;

if directflag  % drag line in main axes
    if logscale(1)
        xlim = 10.^ (log10(oxlim) - log10(np(1)/p(1)));
    else
        xlim = oxlim - (np(1)-p(1));
    end
    if logscale(2)
        ylim = 10.^ (log10(oylim) - log10(np(2)/p(2)));
    else
        ylim = oylim - (np(2)-p(2));
    end
else  % drag overview patch
    if logscale(1)
        xlim = 10.^ (log10(oxlim) + log10(np(1)/p(1)));
    else
        xlim = oxlim + (np(1)-p(1));
    end
    if logscale(2)
        ylim = 10.^ (log10(oylim) + log10(np(2)/p(2)));
    else
        ylim = oylim + (np(2)-p(2));
    end
end

if ~isempty(bounds)
    xlim1 = inbounds(xlim,bounds.xlim,logscale(1));
    ylim1 = inbounds(ylim,bounds.ylim,logscale(2));

    if ~isequal(xlim,xlim1)
        if directflag
            if logscale(1)
                np(1) = 10.^( log10(np(1)) - log10(xlim1(1)/xlim(1)) );
            else
                np(1) = np(1)-(xlim1(1)-xlim(1));
            end
        else
            if logscale(1)
                np(1) = 10.^( log10(np(1)) + log10(xlim1(1)/xlim(1)) );
            else
                np(1) = np(1)+(xlim1(1)-xlim(1));
            end
        end
        xlim = xlim1;
    end
    if ~isequal(ylim,ylim1)
        if directflag
            if logscale(2)
                np(2) = 10.^( log10(np(2)) - log10(ylim1(1)/ylim(1)) );
            else
                np(2) = np(2)-(ylim1(1)-ylim(1));
            end
        else
            if logscale(1)
                np(2) = 10.^( log10(np(2)) + log10(ylim1(1)/ylim(1)) );
            else
                np(2) = np(2)+(ylim1(1)-ylim(1));
            end
        end
        ylim = ylim1;
    end
end

function saveCallbacks = installCallbacks(h,fig,callbackList,valueList)
% installCallbacks
%   inputs:
%      h - handle of object which will be changed by callbacks
%      fig - handle of figure
%      callbackList - list of figure callbacks in cell array
%           elements are e.g., 'windowbuttonmotionfcn'
%      valueList - same length as callbackList - cell array containing
%           values  (string or numeric) for h's userdata
%   outputs:
%      saveCallbacks - cellarray of what the callbacks were before

saveCallbacks = cell(1,length(callbackList));
for i=1:length(callbackList)
    if ischar(valueList{i})
        vstr = ['''' valueList{i} ''''];
    else
        vstr = num2str(valueList{i});
    end
    if 0,  % if problems with fig not being gcf, set this to 1
        figstr = ['hex2num(''' sprintf('%bx',h) ''')']; %#ok<UNRCH>
    else
        figstr = 'gcf';
    end
    str = ['set(findall(' figstr ',' ...
        '''tag'',''' get(h,'tag') '''),''userdata'',' ...
        vstr ')'];
    saveCallbacks{i} = get(fig,callbackList{i});
    set(fig,callbackList{i},str)
end

%
function event = waitForNextEvent(h)
% waitForNextEvent

set(h,'userdata',0)
waitfor(h,'userdata')
event = get(h,'userdata');


function switchToBorderAxes(ax,borderaxes)
% This function hides the main axes and shows the border axes.
% IF borderaxes is NOT EMPTY:
% 		It hides the main axes by setting its x and ycolor to the
% 		color of the axes, instead of turning the mainaxes invisible -
% 		 this is because if we just set visible to off, any lines
% 		 in background erasemode would erase in the FIGURE's background
% 		 color, not the axes' (and that would be bad).
% ELSE
%       Just turns off the axes ticks.
if ~isempty(borderaxes)
    ax_color = get(ax,'color');
    if ischar(ax_color)
        % use figure's color in case axes's color is 'none'
        ax_color = get(get(ax,'parent'),'color');
    end
    set(ax,'xcolor',ax_color,'ycolor',ax_color)
    set(get(ax,'xlabel'),'color',get(borderaxes,'xcolor'))
    set(get(ax,'ylabel'),'color',get(borderaxes,'ycolor'))
    set(borderaxes,'visible','on')
end
set(ax,'xtick',[],'ytick',[])

function switchToMainAxes(ax,borderaxes)
% This function hides the border axes and shows the main axes.
% It undoes what switchToBorderAxes does.
% Assumes xcolor and ycolor of ax should be the same as those of borderaxes
if ~isempty(borderaxes)
    set(ax,'xcolor',get(borderaxes,'xcolor'),...
        'ycolor',get(borderaxes,'ycolor'))
    set(borderaxes,'visible','off')
end
set(ax,'xtickmode','auto','ytickmode','auto')


function doDynamicDrag(xlim,ylim,actual_xlim,actual_ylim,data,...
    transform)


% needed by setDataWithTransform:
xl = [xlim(:)' actual_xlim];
yl = [ylim(:)' actual_ylim];

hh=[]; xx={}; yy={};
for i=1:length(data)
    if data(i).Fs>0
        % translate interval to integer indices
        xlim1 = xlim*data(i).Fs + (1-data(i).Fs*data(i).t0);
        xlim1 = [floor(xlim1(1)) ceil(xlim1(2))];
        ind = max(1,xlim1(1)):min(xlim1(2),size(data(i).data,1));
        if ~isempty(ind)
            if ind(1) == 0, ind(1) = []; end
            x = (ind-1)/data(i).Fs + data(i).t0;
            for j = 1:length(data(i).h)
                y = data(i).data(ind,j);
                if ~isempty(transform)
                    y = feval(transform,y);
                end
                hh = [hh; data(i).h(j)];  %#ok<AGROW>
                xx = {xx{:} x};
                yy = {yy{:} y};
            end
        else
            hh = [hh; data(i).h(:)];  %#ok<AGROW>
            for j = 1:length(data(i).h)
                xx = {xx{:} []};
                yy = {yy{:} []};
            end
        end
    else
        for j = 1:length(data(i).h)
            x = data(i).xdata(:,j);
            y = data(i).data(:,j);
            if ~isempty(transform)
                y = feval(transform,y);
            end
            hh = [hh; data(i).h(j)]; %#ok<AGROW>
            xx = {xx{:} x};
            yy = {yy{:} y};
        end
    end
end
% wait until end to set data, so all lines are erased, and THEN
% redrawn.  This prevents background erasures from erasing other
% previously drawn lines.
setDataWithTransform(hh,xl,yl,xx',yy')

function setDataWithTransform(h,xlim,ylim,xd,yd)
%setDataWithTransform Set xdata and ydata of lines transforming from one
%   set of limits to another.
%   Inputs:
%      h - vector of line handles
%      xlim, ylim - limits of mainaxes
%         if xlim and ylim are 4 elements long, uses the 3rd and 4th
%         elements as the interval in which to map the data.
%      xd - vector or cell array.  Each element of the cell
%            will get mapped from xlim to [0 1] or xlim(3:4)
%      yd - vector or cell array.  Each element of the cell
%            will get mapped from ylim to [0 1] or ylim(3:4)

if length(xlim)<=2
    xslope = 1/diff(xlim);
    xintercept = -xlim(1)*xslope;
else
    xslope = (xlim(4)-xlim(3))/(xlim(2)-xlim(1));
    xintercept = xlim(3)-xlim(1)*xslope;
end
if length(ylim)<=2
    yslope = 1/diff(ylim);
    yintercept = -ylim(1)*yslope;
else
    yslope = (ylim(4)-ylim(3))/(ylim(2)-ylim(1));
    yintercept = ylim(3)-ylim(1)*yslope;
end

if ~iscell(xd)
    xd = {xd};
end
if ~iscell(yd)
    yd = {yd};
end
for j=1:length(h)
    xd{j} = xslope*xd{j}+xintercept;
    yd{j} = yslope*yd{j}+yintercept;
end
set(h,{'xdata'},xd,{'ydata'},yd)

%=================== Zoom functions ======================================
function msvzoom(varargin)
%MSVZOOM  MSViewer zoom function.
%  Contains callbacks for Zoom button group of MSViewer.
%    zoomxy
%    zoomoutxy
%    zoominx
%    zoominy
%    zoomin
%    zoomout
if nargin > 3
    fig = varargin{4};
else
    fig = gcbf;
end
ud = get(fig,'userdata');
downtype = 'mousedown';
action = varargin{3};
if nargin > 4
    downtype = varargin{5};
end
%xxx=============================
switch action
    case 'mousezoom'
        %         ud = managemarkerlist(ud, 'inzoom');
        %         set(fig,'userdata',ud)
        set(fig,'windowbuttondownfcn',{@msvzoom, 'mousedown', fig, 'zoomxydown'})
        ud.pointer = 1;
        ud.zoompointer = 'cross';
        %         set(ud.contextmenu.addmarker, 'visible', 'off')
        set(fig,'userdata',ud)
        setptr(fig,'arrow')

    case 'resetview'
        msoverview('zoom',ud.limits.xlim,ud.limits.ylim, fig)
        set(ud.mainaxes,'xlim',ud.limits.xlim,'ylim',ud.limits.ylim)
        %         markermodeButtonChange([], [], fig)

    case 'zoominx'
        %         ud = managemarkerlist(ud, 'inzoom');
        %         set(fig,'userdata',ud)
        set(fig,'windowbuttondownfcn', {@msvzoom, 'mousedown', fig, 'zoomxdown'});
        ud.pointer = 1;
        ud.zoompointer = 'lrdrag';
        set(fig,'userdata',ud)
        setptr(fig,'arrow')

    case 'zoominy'
        %         ud = managemarkerlist(ud, 'inzoom');
        %         set(fig,'userdata',ud)
        set(fig,'windowbuttondownfcn',{@msvzoom, 'mousedown', fig, 'zoomydown'})
        ud.pointer = 1;
        ud.zoompointer = 'uddrag';
        set(fig,'userdata',ud)
        setptr(fig,'arrow')

    case 'zoomout'
        xlim = get(ud.mainaxes,'xlim');
        ylim = get(ud.mainaxes,'ylim');

        xlim = .5*[3 -1]*xlim' + [0 diff(xlim)*2];
        xlim = [max(xlim(1),ud.limits.xlim(1)), min(xlim(2),ud.limits.xlim(2))];

        ylim = .5*[3 -1]*ylim' + [0 diff(ylim)*2];
        ylim = [max(ylim(1),ud.limits.ylim(1)), min(ylim(2),ud.limits.ylim(2))];
        msoverview('zoom', xlim, ylim, fig)
        set(ud.mainaxes,'xlim',xlim)
        set(ud.mainaxes,'ylim',ylim)
        %         markermodeButtonChange([], [], fig)

        %-------------- these are self callbacks:
    case 'mousedown'
        save_wbmf = get(fig,'windowbuttonmotionfcn');
        %         hZoomTool = [];
        switch get(fig,'selectiontype')
            case 'normal'     % do zoom
                ud.justzoom = get(fig,'currentpoint');
                set(fig,'userdata',ud)

                pstart = get(fig,'currentpoint');

                % don't do anything if click is outside mainaxes_port
                mp=getpixelposition(ud.mainaxes);
                %Click is outside of main panel:
                if ~pointsinrect(pstart,[mp(1) mp(1)+mp(3) mp(2) mp(2)+mp(4)])
                    return
                end

                set(fig,'windowbuttonmotionfcn', '');
                start = get(ud.mainaxes, 'CurrentPoint');
                if isnan(start)  % case in which start is outside of axes
                    return
                end
                ud.cursorposition=start;

                switch downtype
                    case 'zoomxydown'
                        r=rbbox([pstart 0 0],pstart);

                        oldxlim = get(ud.mainaxes,'xlim');
                        oldylim = get(ud.mainaxes,'ylim');

                        if any(r([3 4])==0)
                            % zero width or height, or just a click - stay in zoom
                            % mode and try again
                            return

                        else
                            % zoom to the rectangle dragged
                            set(fig,'currentpoint',[r(1) r(2)])
                            p1 = get(ud.mainaxes,'currentpoint');
                            set(fig,'currentpoint',[r(1)+r(3) r(2)+r(4)])
                            p2 = get(ud.mainaxes,'currentpoint');

                            xlim = [p1(1,1) p2(1,1)];
                            ylim = [p1(1,2) p2(1,2)];
                        end
                        newxlim = xlim;
                        newylim = ylim;
                        newxlim = fitwithinXlimit(ud, newxlim);
                        if diff(newxlim) > 0
                            msoverview('zoom', newxlim,oldylim, fig)
                            set(ud.mainaxes,'xlim',newxlim)
                        else
                            newxlim = oldxlim;
                        end
                        if diff(newylim) > 0
                            msoverview('zoom', newxlim,newylim, fig)
                            set(ud.mainaxes,'ylim',newylim)
                        end

                        %                         hZoomTool = findobj(ud.togglebuttongroup, 'Tag','zoom xy');
                    case 'zoomxdown'
                        set(fig, 'userdata', ud);
                        set(fig,'windowbuttonupfcn', {@zoom_buttonup})

                        ylim = get(ud.mainaxes,'ylim');
                        xd = start(1,[1 1 1 1 1]);
                        yd = ylim([1 1 2 2 1]);
                        bg_patch =  patch(xd,yd,[1 1 1],...
                            'FaceAlpha', 0.25,...
                            'FaceColor',  [0.85 0.95 1],...
                            'Parent', ud.mainaxes,...
                            'edgecolor',[0 0 1],...
                            'tag', 'msvbgpatch','visible', 'off');

                        set(fig, 'WindowButtonMotionFcn', {@zoomx_move, bg_patch, start(1,1), fig});
                        set(ud.overview.oviewpatch, 'userdata', 0);
                        waitfor(ud.overview.oviewpatch,'userdata',1) %button up
                        ud=get(fig, 'userdata');
                        if ud.cursorposition(1,1) ~= start(1,1)
                            end_xlim=sort([start(1,1), ud.cursorposition(1,1)]);
                            end_xlim = fitwithinXlimit(ud, end_xlim);
                            if diff(end_xlim) > 0
                                set(ud.mainaxes,'xlim',end_xlim)
                                msoverview('zoom', end_xlim, get(ud.mainaxes, 'ylim'), fig);
                            end
                        end
                        ud.cursorposition =[];
                        delete(bg_patch)
                        %                         hZoomTool = findobj(ud.togglebuttongroup, 'Tag', 'zoom x');
                    case 'zoomydown'
                        set(fig, 'userdata', ud);
                        set(fig,'windowbuttonupfcn', {@zoom_buttonup})
                        xlim = get(ud.mainaxes,'xlim');
                        xd = xlim([1 2 2 1 1]);
                        yd = start(1,[2 2 2 2 2]);

                        bg_patch =  patch(xd,yd,[1 1 1],...
                            'FaceAlpha', 0.25,...
                            'FaceColor', [0.85 0.95 1],...
                            'Parent', ud.mainaxes,...
                            'edgecolor',[0 0 1],'visible', 'off');

                        set(fig, 'WindowButtonMotionFcn', {@zoomy_move, bg_patch, start(1,2), fig});
                        set(ud.overview.oviewpatch, 'userdata', 0);
                        waitfor(ud.overview.oviewpatch,'userdata',1)
                        ud=get(fig, 'userdata');

                        if ud.cursorposition(1,2) ~= start(1,2)
                            end_ylim=sort([start(1,2), ud.cursorposition(1,2)]);
                            set(ud.mainaxes,'ylim',end_ylim)
                            msoverview('zoom', xlim, end_ylim, fig);
                        end
                        % OK we're back from the drag
                        ud.cursorposition =[];
                        delete(bg_patch)
                        %                         hZoomTool = findobj(ud.togglebuttongroup, 'Tag', 'zoom y');
                end
                %                 set(hZoomTool, 'state', 'off');
                %                 toggleButtonChange(hZoomTool,[])
                set(fig,'windowbuttonmotionfcn',save_wbmf);
                set(fig,'windowbuttonupfcn','');
                set(fig,'userdata',ud)
                set(fig,'currentpoint',ud.justzoom)

            otherwise
                % do nothing!
        end
        %         set(hZoomTool, 'state', 'off');
        %         toggleButtonChange(hZoomTool,[])
    case 'leavezoom'
        set(fig,'windowbuttondownfcn','')
        ud.pointer = 0;
        set(ud.contextmenu.addmarker, 'visible', 'on')
        ud = managemarkerlist(fig, ud, 'leavezoom', true); %enable marker action
        set(fig,'userdata',ud)
end
ud=get(fig, 'userdata');
xlim = get(ud.mainaxes,'xlim');
ylim = get(ud.mainaxes,'ylim');
idx = get(ud.markerlist,'Value');

focus_marker=findmarkeridex(ud, idx);

for i=1:length(ud.markers)
    if ud.markers(i).mzvalue > xlim(1) && ud.markers(i).mzvalue < xlim(2)
        if i==focus_marker
            set(ud.markers(i).markerlabel, 'visible', 'on');
            xdata=get(ud.markers(i).max_markerline, 'XData');
            positionmarkerlabel(ud.markers(i).markerlabel, xdata(2), ylim(2))
        end
    else
        set(ud.markers(i).markerlabel, 'visible', 'off');
    end
end
set(fig, 'userdata', ud);
%-------------------------------------------------------------
function zoomx_move(hSrc, event, bg_patch, x_start, hfig) 
% Windowbuttonmotionfcn for dragging zoom along x axis
ud=get(hSrc, 'userdata');
p = get(ud.mainaxes,'currentpoint');
if isnan(p)  % case in which start is outside of axes
    return
end
x=p(1,1);
xlim = sort([x_start x]);

if diff(xlim)>0
    setpdata(bg_patch, xlim,get(ud.mainaxes, 'ylim'));
    set(bg_patch, 'visible', 'on')
    msoverview('zoom', xlim, get(ud.mainaxes, 'ylim'), hfig);
end
ud.cursorposition=p;
set(hSrc, 'userdata', ud)

function zoomy_move(hSrc, event, bg_patch, y_start, hfig)
% Windowbuttonmotionfcn for dragging zoom along x axis
ud=get(hSrc, 'userdata');
p = get(ud.mainaxes,'currentpoint');
if isnan(p)  % case in which start is outside of axes
    return
end
y=p(1,2);
ylim = sort([y_start y]);
ud.cursorposition=p;
if diff(ylim) > 0
    setpdata(bg_patch, get(ud.mainaxes, 'xlim'), ylim);
    set(bg_patch, 'visible', 'on')
    msoverview('zoom', get(ud.mainaxes, 'xlim'), ylim, hfig);
end
set(hSrc, 'userdata', ud)

function zoom_buttonup(hSrc, event, zoomtag)  %#ok<INUSD>
% Windowbuttonmotionfcn for dragging zoom along x axis
set(hSrc, 'WindowButtonMotionFcn', '')
ud=get(hSrc, 'userdata');
currentp = get(ud.mainaxes,'currentpoint');
if isnan(currentp)  % case in which start is outside of axes
    return
end

ud.cursorposition = currentp;
set(ud.overview.oviewpatch, 'userdata',1)
% get out of zoom mode
set(hSrc, 'userdata', ud)
%=================== Print functions =====================================
function msprinttofigure(hFigure)
%MSPRINTTOFIGURE Print Mass spectra to a figure.
%  MSPRINTTOFIGURE(hFigure) prints the Mass spectra on main
%  axes to  a new figure.  The new figure window is centered on the screen.

old_fig = hFigure;
old_axes = findobj(old_fig, 'type', 'axes', 'tag', 'mainaxes');
old_fig_pos = get(old_fig,'Position');

fig_width = old_fig_pos(3);
fig_height = old_fig_pos(4);

screensize = get(0,'screensize');

fig_left   = round((screensize(3) - fig_width) / 2);
fig_bottom = round((screensize(4) - fig_height) / 2);

fig_position = [fig_left fig_bottom fig_width fig_height];

new_fig = figure('Visible', 'off', ...
    'Units', 'pixels', ...
    'Position', fig_position);

new_axes = axes('Parent',new_fig,...
    'Units',get(old_axes,'Units'));
newpos = get(new_axes,'Position');
delete(new_axes)

h_axes = copyobj(old_axes, new_fig);
xlim = get(old_axes, 'Xlim');
ylim = get(old_axes,'Ylim');

set(h_axes, 'Position', newpos, ...
    'HandleVisibility','on',...
    'XLim', xlim,...
    'YLim', ylim);

lineobjs = findobj(h_axes, 'type', 'line');
set(lineobjs,'ButtonDownFcn', '',...
    'UicontextMenu','');
set(new_fig, 'Visible', 'on');

%-------------------------------------------------------------------
function msprintfcn(hSrc, event, h_fig, action) 
%MSPRINTFCN Create and display a print preview of MSViewer spectra.

ud = get(h_fig, 'userdata');

if ~isempty(ud.tempfig) && ishandle(ud.tempfig)
    delete(ud.tempfig)
    ud.tempfig = [];
end

ud.tempfig = figure('Visible','off');
copyobj(ud.mainaxes, ud.tempfig);

switch action
    case 'pagesetup'
        pagesetupdlg(ud.tempfig);
    case 'preview'
        printpreview(ud.tempfig);
    case 'print'
        printdlg(ud.tempfig); %#ok<MCPRD>
        close(ud.tempfig);
end

set(h_fig, 'userdata', ud)

%=================== Help functions ======================================
function bioinfostandardhelp(helpmenu)
% BIOINFOSTANDARDHELP Add Toolbox, Demos, and About to help menu.
uimenu(helpmenu, 'Label', 'Bioinformatics Toolbox Help', ...
    'Callback', 'helpview(fullfile(docroot,''toolbox'',''bioinfo'',''bioinfo.map''),''bioinfo_product_page'')');
uimenu(helpmenu, 'Label', 'MS Viewer Help', ...
    'Callback', 'helpview(fullfile(docroot,''toolbox'',''bioinfo'', ''bioinfo.map''),''msviewer_refpage'')');
uimenu(helpmenu, 'Label', 'Examples', ...
    'Callback', @(varargin) demo('toolbox','bioinformatics'), ...
    'Separator', 'on');
tlbx = ver('bioinfo');
tlbx = tlbx(1);
mailstr = ['web(''mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20',...
    'MSViewer%20in%20Bioinformatics%20Toolbox%20', tlbx.Version,''')'];
uimenu(helpmenu, 'Label', 'Send Feedback...', 'Separator', 'on',...
    'Callback', mailstr);
uimenu(helpmenu, 'Label', 'About Bioinformatics Toolbox', ...
    'Callback', @bioinfoabout, ...
    'Separator', 'on');



function bioinfoabout(varargin)
%BIOINFOABOUT About the Bioinformatics Toolbox.

tlbx = ver('bioinfo');
tlbx = tlbx(1);
str = sprintf('%s %s\nCopyright 1993-%s The MathWorks, Inc.', ...
    tlbx.Name, tlbx.Version, datestr(tlbx.Date, 10));
msgbox(str,tlbx.Name,'modal');

%=================== Import/Export markers function ======================
function export_successful = exportMarkersToWS(hSrc, event, hfig)
%Exports markers to workspace

export_successful = false;

% Define the prompt
dlg_prompt = {'Markers variable name:'};

% Define the title
dlg_title = 'Export Markers To Workspace';
num_lines_per_prompt = 1;

while ~export_successful
    % Display the input dialog using INPUTDLG from MATLAB
    user_var = inputdlg(dlg_prompt{1}, dlg_title, num_lines_per_prompt);

    user_pressed_cancel = isempty(user_var);
    if user_pressed_cancel
        export_successful = false;
        return;
    end

    %         hfig = gcbf;
    ud = get(hfig, 'userdata');

    marker_data = zeros(size(ud.markers));

    for i =1:length(ud.markers)
        marker_data(i)= ud.markers(i).mzvalue;
    end
    marker_data = {sort(marker_data(:))}; % in a cell array
    export_successful = assignvartobaseworkspace(user_var, marker_data);
end
%------------------------------------------------
function importMarkersFromWS(hSrc, event, hfig)
%Import markers from workspace
%     hfig = gcbf;
getvarfromworkspace(hfig);
ud = get(hfig, 'userdata');

if ud.importcanceled
    return
end
x=[];
marker=[];
markers=sort(ud.importmarkers);
if ~isempty(ud.markers)
    marker=ud.markers(1);
    x = zeros(size(ud.markers));
    for i =1:length(ud.markers)
        x(i)= ud.markers(i).mzvalue;
    end

end
for i=1:length(ud.importmarkers)
    idx = markers(i,:)==x;
    if isempty(idx)
        marker.mzvalue = markers(i,:);
        marker.max_markerline=[];
        marker.oax_markerline=[];
        marker.markerlabel=[];
        marker.cmenu=[];
        ud.markers = [ud.markers marker];
    end
end
set(hfig, 'userdata', ud);
ud = msvmarkers([], [], 'init', hfig);
set(hfig, 'userdata', ud);

%==============================================

function getvarfromworkspace(varargin)
if nargin == 1
    % client needs to be a figure
    hfig = varargin{1};
end

% Output variables for function scope
uData.last_selected_value = [];
uData.hvarlist=[];
% Get workspace variables and store variable names
workspace_vars = evalin('base','whos');
num_of_vars = length(workspace_vars);
[uData.all_var_names{1:num_of_vars}] = deal(workspace_vars.name);

%== Create import figure
hImportFig = figure('Toolbar','none',...
    'Menubar','none',...
    'NumberTitle','off',...
    'IntegerHandle','off',...
    'Tag','importMarkerFromWS',...
    'Visible','off',...
    'HandleVisibility','callback',...
    'Name','Import Markers From Workspace',...
    'WindowStyle','modal',...
    'Resize','off');

% Layout management
fig_height = 300;
fig_width  = 300;
fig_size = [fig_width fig_height];
left_margin = 10;
right_margin = 10;
bottom_margin = 10;
top_margin = 10;
spacing = 5;
button_size = [60 25];
default_panel_width = fig_width - left_margin -right_margin;

screen_size = get(0,'ScreenSize');
lower_left_pos = 0.5 * (screen_size(3:4) - fig_size);
set(hImportFig,'Position',[lower_left_pos fig_size]);

%===Create button panel
btnPanelPos = [left_margin, bottom_margin, default_panel_width, button_size(2)];

hButtonPanel = uipanel('parent',hImportFig,...
    'Tag','buttonPanel',...
    'Units','pixels',...
    'Position',btnPanelPos,...
    'BorderType','none');

set(hButtonPanel,'BackgroundColor',get(hImportFig,'Color'));

% add buttons
button_strs_n_tags = {'OK', 'okButton';...
    'Cancel','cancelButton'};

num_of_buttons = length(button_strs_n_tags);

button_spacing = (btnPanelPos(3)-(num_of_buttons * button_size(1)))/(num_of_buttons+1);
posX = button_spacing;
posY = 0;
buttons = zeros(num_of_buttons,1);

for n = 1:num_of_buttons
    buttons(n) = uicontrol('parent',hButtonPanel,...
        'Style','pushbutton',...
        'String',button_strs_n_tags{n,1},...
        'Tag',button_strs_n_tags{n,2});
    set(buttons(n),'BackgroundColor',get(hButtonPanel,'BackgroundColor'));

    set(buttons(n),'Position',[posX, posY, button_size]);
    set(buttons(n),'Callback',@doButtonPress);
    posX = posX + button_size(1) + button_spacing;
end

%===Create variable panel
height = fig_height - top_margin - (button_size(2)+ bottom_margin + 2*spacing);
panelPos = [left_margin, button_size(2)+ bottom_margin + spacing,...
    default_panel_width, height];

hVarPanel = uipanel('parent',hImportFig,...
    'Tag','markerVarPanel',...
    'Units','pixels',...
    'BorderType','none',...
    'Position',panelPos,...
    'Visible','off');
set(hVarPanel,'BackgroundColor',get(hImportFig,'Color'));

hVarLabel = uicontrol('parent',hVarPanel,...
    'Style','text',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'String','Variables:');
set(hVarLabel,'BackgroundColor',get(hVarPanel,'BackgroundColor'));

label_extent = get(hVarLabel,'Extent');
label_posX = left_margin;
label_posY = panelPos(4) - label_extent(4) - spacing;
label_width = label_extent(3);
label_height = label_extent(4);
label_position = [label_posX label_posY label_width label_height];

set(hVarLabel,'Position',label_position);

hVarList = uicontrol('parent',hVarPanel,...
    'Style','listbox',...
    'fontname','Courier',...
    'Value',1,...
    'Units','pixels',...
    'Tag', 'markerList',...
    'BackgroundColor','white');

list_posX = left_margin;
list_posY = bottom_margin;
list_width = panelPos(3) - 2*list_posX;
list_height = panelPos(4) - list_posY - label_height - spacing;
list_position = [list_posX list_posY list_width list_height];

set(hVarList,'Position',list_position);
set(hVarList,'Callback',@listSelected);

%== Filter the right type of variable
varInd = [];

for n = 1:num_of_vars
    is_double = strcmpi(workspace_vars(n).class,'double');
    if is_double
        is_M_by_1 = (length(workspace_vars(n).size) == 2 && workspace_vars(n).size(end) == 1);
        is_1_by_M = (length(workspace_vars(n).size) == 2 && workspace_vars(n).size(1) == 1);
        if is_M_by_1 || is_1_by_M
            varInd = [varInd,n];  %#ok<AGROW>
        end
    end
end

%== Display the variables
ws_vars = orderfields(workspace_vars(varInd));
num_of_vars = length(ws_vars);

var_str = cell(num_of_vars,1);
var_names = {ws_vars.name};
longest_var_name = max(cellfun('length',var_names));
format1 = sprintf('%%-%ds',longest_var_name+2);
format2 = sprintf('%%-12s %%-6s');
format_all = sprintf('%s%s',format1,format2);

for n = 1:num_of_vars
    if length(ws_vars(n).size) == 3
        sz_str = sprintf('%dx%dx%d',ws_vars(n).size);
        tmp_str= sprintf(format_all,ws_vars(n).name,...
            sz_str,ws_vars(n).class);
    else
        sz_str = sprintf('%dx%d',ws_vars(n).size);
        tmp_str= sprintf(format_all,ws_vars(n).name,...
            sz_str,ws_vars(n).class);
    end
    var_str{n} = sprintf('%s\n',tmp_str);
end

set(hVarList,'String',var_str);
set(hVarList,'HorizontalAlignment','left')
set(hVarPanel,'Visible','on');
set(hImportFig,'Visible','on');
uData.hvarlist = hVarList;
uData.mainfig = hfig;
set(hImportFig, 'userdata', uData);

% This blocks until the user explicitly closes the tool.
uiwait(hImportFig);

%+======================================================================
function listSelected(hSrc,evt) %#ok<INUSD>
% callback for the  list boxes
% we disable the colormap panel controls for an RGB image
himportfig = gcbf;
ud = get(himportfig, 'userdata');
mainud = get(ud.mainfig, 'userdata');

ind = get(hSrc,'Value');
list_str = get(hSrc,'String');

if isempty(list_str)
    return
end

double_click = strcmpi(get(himportfig,'SelectionType'),'open');
[issuccess, ud] =  getVars(ud);
clicked_same_list_item = (ud.last_selected_value == ind);
if double_click && clicked_same_list_item && issuccess
    mainud.importcanceled = false;
    set(ud.mainfig, 'userdata', mainud)
    ud.last_selected_value = ind;
    set(himportfig, 'userdata', ud);
    close(himportfig);
else
    set(himportfig,'SelectionType','normal');
    ud.last_selected_value = ind;
    set(himportfig, 'userdata', ud);
end


%------------------------------------------------
function doButtonPress(hSrc,evt) %#ok<INUSD>
% call back function for the OK and Cancel buttons
tag = get(hSrc,'tag');
himportfig = gcbf;
ud = get(himportfig, 'userdata');

switch tag
    case 'okButton'
        [issucess, ud] = getVars(ud);
        if issucess
            mainud = get(ud.mainfig, 'userdata');
            mainud.importcanceled = false;
            set(ud.mainfig, 'userdata', mainud)
            set(himportfig, 'userdata', ud);
            close(himportfig);
        end

    case 'cancelButton'
        close(himportfig);
end


%------------------------------------------------
function [issuccess, ud] = getVars(ud)
issuccess = true;
list_str = get(ud.hvarlist,'String');

% return if there are no variables listed in current panel
if isempty(list_str)
    error_str = sprintf('There are no variables containing marker');
    error_str = sprintf('%s data in the workspace.',error_str);
    errordlg(error_str);
    issuccess = false;
    return;
end

ind = get(ud.hvarlist,'Value');
var_name = strtok(list_str{ind});
try
    importmarkers = evalin('base',sprintf('%s;',var_name));
    importmarkers=importmarkers(:);
catch theException
    error_str = theException.message;
    errordlg(error_str)
    issuccess = false;
    return;
end

% check that the value is in range
mainud = get(ud.mainfig, 'userdata');
if max(importmarkers)>max(mainud.xdata)...
        || min(importmarkers)<min(mainud.xdata)
    error_str = sprintf('The variable you selected is out of range.');
    errordlg(error_str);
    issuccess = false;
    return;
elseif length(importmarkers) >= length(mainud.xdata)
    error_str = sprintf('The variable you selected is not suitable for markers.');
    errordlg(error_str);
    issuccess = false;
    return;
end
mainud.importmarkers = importmarkers ;
set(ud.mainfig, 'userdata', mainud);


%=================== Utility functions ===================================
function [d1, d2] = msvicondir
%   (D1) and the MATLAB icons (D2).
pathname = fileparts(which(mfilename));
% pathname = fileparts(pathname);
d1 = fullfile(pathname, 'icons');
if nargout > 1
    d2 = fullfile(toolboxdir('matlab'), 'icons'); 
end


function varargout = msvswitch(varargin)
%MSVSWITCH Function switch-yard.
%    MSVSWITCH('FOO',ARG1,ARG2,...) is the same as FOO(ARG1,ARG2,...).  This
%    provides access to private functions for Handle Graphics callbacks.
if (nargout == 0)
    feval(varargin{:});
else
    [varargout{1:nargout}] = feval(varargin{:});
end


function setpdata(hpatch,xlim,ylim)
%setpdata - set xdata and ydata of patch object to rectangle specified by
% xlim and ylim input
set(hpatch,'xdata',[xlim(1) xlim(2) xlim(2) xlim(1) xlim(1)], ...
    'ydata',[ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)]) % thumb patch


function bool = pointsinrect(pts,rect)
%POINTSINRECT Determine if points lie in or on rectangle.
%   Inputs:
%     pts - n-by-2 array of [x,y] data
%     rect - 1-by-4 vector of [xlim ylim] for the rectangle
%   Outputs:
%     bool - length n binary vector
[i,j] = find(isnan(pts));  %#ok<NASGU>
bool = (pts(:,1)<rect(1))|(pts(:,1)>rect(2))|...
    (pts(:,2)<rect(3))|(pts(:,2)>rect(4));
bool = ~bool;
bool(i) = 0;

function x = inbounds(x,lim,logscale)
%INBOUNDS Returns value limited to interval.
%   inbounds(x,lim), where x is a scalar and lim = [lower_bound
%              upper_bound], clips x to the interval lim.
%   inbounds(interval,lim), where interval is an interval with
%          length less than or equal to lim, clips the interval
%          to lim, keeping its width constant.
%   inbounds(interval,lim,1) clips the interval to lim, keeping the log of its
%      width constant.
if nargin < 3
    logscale = 0;
end

if length(x)==1
    if x<lim(1)
        x = lim(1);
    elseif x>lim(2)
        x = lim(2);
    end
else
    if x(1)<lim(1)
        if logscale
            x = 10.^(log10(x) - log10(x(1)/lim(1)));
        else
            x = x - x(1) + lim(1);
        end
    elseif x(2)>lim(2)
        if logscale
            x = 10.^(log10(x) - log10(x(2)/lim(2)));
        else
            x = x - x(2) + lim(2);
        end
    end
    if diff(x)>diff(lim)
        x = lim;
    end
end

%---------------------------------------------------------------------
function cdata = makeToolbarIcon(filename,guessalpha)
% MAKETOOLBARICON read an image file and convert it to CData for a toolbar icon.
%
% CDATA=MAKETOOLBARICON(FILENAME)
%   Read an image file and convert it to CData with automatic transparency
%   handling. If the image has transparency data, PNG files sometimes do,
%   the transparency data is used. If the image has no CData, the top left
%   pixel is treated as the transparent color.
%
% CDATA=MAKETOOLBARICON(FILENAME, FALSE)
%   Same as above but suppress the usage of the top left pixel for images
%   with no transparency data. This may require the caller to handle the
%   transparency explicitly. View the contents of this Matlab file for an
%   example of how to handle transparency.
%
% Adopted from Bill York's iconread.m


if nargin < 2
    guessalpha = true;
end

[p,f,ext] = fileparts(filename); 
% if this is a mat-file, look for the variable cdata (or something like it)
if isequal(lower(ext),'.mat')
    cdata = [];
    s = whos('-file',filename);
    for i=1:length(s)
        if ~isempty(strfind(lower(s(i).name), 'cdata'))
            data = load(filename,s(i).name);
            cdata = data.(s(i).name);
        end
    end
    return
end

[cdata,map,alpha] = imread(filename);
if isempty(cdata)
    return;
end

if isempty(map)
    % need to use doubles cuz nan's only work as doubles
    cdata = double(cdata);
    cdata = cdata/255;
else
    cdata = ind2rgb(cdata,map);
end

if isempty(alpha)
    if ~guessalpha
        return;
    end
    % guess the alpha pixel by using the top left pixel in the icon
    ap1 = cdata(1,1,1);
    ap2 = cdata(1,1,2);
    ap3 = cdata(1,1,3);
    alpha = cdata(:,:,1) == ap1 & cdata(:,:,2) == ap2 & cdata(:,:,3) == ap3;
    alpha = ~alpha;
end

% process alpha data
r = cdata(:,:,1);
r(alpha == 0) = NaN;
g = cdata(:,:,2);
g(alpha == 0) = NaN;
b = cdata(:,:,3);
b(alpha == 0) = NaN;
cdata = cat(3,r,g,b);
%-----------------------------------------------------------------------
function successful = assignvartobaseworkspace(variable_name, variable_data)
% Assigns variable_data to variable_name in the base workspace
% Both arguments are cell arrays of the same length

if isempty(variable_name) || isempty(variable_data)
    str1 = sprintf('You did not specify any variable names.');
    uiwait(errordlg(str1,'Export To Workspace'));
    successful = false;
    return
end

% Grab all the workspace variables and store their names in a cell
% array.
ws_vars = evalin('base','whos');
[ws_var_names{1:length(ws_vars)}] = deal(ws_vars.name);

user_response = '';
valid_var_name = genvarname(variable_name);

for n = 1:length(variable_name)

    var_already_exists = any(strcmp(valid_var_name{n}, ws_var_names),2);
    user_spec_var_name_changed = ~strcmpi(variable_name{n}, valid_var_name{n});
    skip_next_check = false;

    if user_spec_var_name_changed
        str = sprintf('"%s" is an invalid variable name.  "%s" is more appropriate.',...
            variable_name{n}, valid_var_name{n});
        question_str = sprintf('%s  Would you like to use the suggested name?',str);

        user_response = questdlg(question_str,'Name change','Yes');

        if strcmpi(user_response,'no')
            skip_next_check = true;
        end

    end

    if ~isempty(var_already_exists) && var_already_exists && ~skip_next_check

        question_str = sprintf('The variable %s already exists in the workspace.\n',...
            valid_var_name{n});

        question_str = sprintf('%sWould you like to overwrite?',question_str);

        user_response = questdlg(question_str, 'Variable already exists','No');

    end

    %Handle user response appropriately
    switch lower(user_response)
        case {'yes',''} % user answered yes or didn't have to specify any responses
            % good! now do it again
        case 'no'
            successful = false;
            return;
        case 'cancel'
            successful = true;
            return
    end

end

for n = 1:length(variable_name)
    assignin('base',valid_var_name{n}, variable_data{n});
    successful = true;
end

%---------------------------------------------------------------------
function newFigPos = setminfigsize(newFigPos, lastFigPos, minSize)
%setMinFigSize sets the minimum size of a figure
%
%   setMinFigSize(newFigPos, lastFigPos, minXYSize) constrains the width
%   and height of a figure by using its last position and new position.
%   newFigPos returned by setMinFigSize contains the new position of the
%   figure such that only the resizing edge/corner of the figure is changed.

newSize = newFigPos(3:4);
tooSmall = (newSize < minSize);
newSize(tooSmall) = minSize(tooSmall);

if tooSmall(1) % X needs to be resized
    changeX = newSize(1) - newFigPos(3);

    %determine the resize direction
    leftSideMoved = (lastFigPos(1) ~= newFigPos(1));

    if (leftSideMoved) % the left side moved
        newFigPos(1) = newFigPos(1)-changeX;
        newFigPos(3) = newSize(1);
    else % the right side moved
        newFigPos(3) = newFigPos(3)+changeX;
    end
end

if tooSmall(2) % Y needs to be resized
    changeY = newSize(2) - newFigPos(4);

    %determine the resize direction
    bottomMoved = (lastFigPos(2) ~= newFigPos(2));
    if (bottomMoved) % the bottom moved
        newFigPos(2) = newFigPos(2) - changeY;
        newFigPos(4) = newSize(2);
    else % the top moved
        newFigPos(4) = newFigPos(4)+changeY;
    end
end





