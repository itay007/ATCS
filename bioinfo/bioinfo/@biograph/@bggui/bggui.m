function h = bggui(bho)
%BGGUI class constructor.


% Copyright 2003-2011 The MathWorks, Inc.

h = biograph.bggui; % constructor

% make a deep copy of the input object
bh = bho.deepCopy;

% create a listener to close gui when object is deleted
bh.DestroyListener = handle.listener(bh,'ObjectBeingDestroyed',@deletebggui);

% if biogrph has not been arranged, do it now
if ~bh.isLaidout
    bh.dolayout
end

% figure out the figure number
figureName = ['Biograph Viewer ' getbiographtoolnumber];

% create a new figure to contain the bggui
fh = figure('visible','off','Name',figureName,...
     'NumberTitle','off','IntegerHandle','off','tag','BioGraphTool');

% correct the figure position
set(fh,'units','points','position',figurePosition(fh,bh))
set(fh,'PaperPositionMode','auto','Color',[1 1 1]);

% save bggui handle in figure
setappdata(fh,'bggui',h)

% save figure handle into bggui object
h.hgFigure = fh;

% setup UIMenus, context menus and toolbar
makeBiographViewerUIMenus(fh,bh);                    
makeBiographViewerToolbar(fh,bh);
[hcm1,hcm2,hcm3] = makeBiographViewerContextMenus(fh,bh);

% connect biograph to bggui
h.biograph = bh;
h.connect(bh,'down');

% setup layout in HG
bh.hgSetup
set(bh.hgAxes,'UIContextMenu',hcm3)

% draw edges in HG (vectorized)
if ~isempty(bh.Edges)
    bh.Edges.hgDraw(hcm2)
end

% draw nodes in HG (vectorized)
bh.Nodes.hgDraw(hcm1);

% set up listeners to catch zoom
              
% listens when the Ylim of axes has changed
zoomListenerHandle = addlistener(gca,'YLim',...
                     'PostSet',@(lH,lEvent)zoomListener(lH,lEvent,fh,bh));
% store the zoom listeners
setappdata(fh,'bgguiListeners',zoomListenerHandle);

% set callback for resizing and make gui visible
set(fh,'ResizeFcn',{@bgguiResize,bh},'visible','on')

% setup figure callback functions
set(fh,'WindowButtonDownFcn',{@mouseClickOnFigure,bh,h});
set(fh,'WindowButtonUpFcn',{@mouseRelease,bh,h});
set(fh,'WindowButtonMotionFcn',{@localWindowButtonMotion,bh,h});

h.dragBox = plot(1,1,':k','LineWidth',1,'Visible','off'); 
h.ghostNodes = patch(1,1,[.5 .5 .5],'LineStyle','none','Visible','off');
h.dataTip = text(0,1,1,'k','Tag','TreeTag','BackgroundColor',[1 1 .93],...
                'Color', [0 0 0],'EdgeColor', [0.8 0.8 0.8],...
                'VerticalAlignment','Bottom','Clipping','off',...
                'Visible','off','Fontsize',8,'Interpreter','none'); 

set(fh,'HandleVisibility','callback') % after all init, make it invisible

bgpk = findpackage('biograph');
nodec = bgpk.findclass('node');
edgec = bgpk.findclass('edge');
biographc = bgpk.findclass('biograph');
listener1 = handle.listener(bh.nodes, nodec.Properties, 'PropertyPostSet', @(h,e) e.AffectedObject.hgUpdate);
listener2 = handle.listener(bh.edges, edgec.Properties, 'PropertyPostSet', @(h,e) e.AffectedObject.hgUpdate);
listener3 = handle.listener(bh, biographc.Properties, 'PropertyPostSet', @(h,e) e.AffectedObject.hgReDraw(e));
setappdata(fh,'biographListeneres',[listener1,listener2,listener3])
setappdata(fh,'biographContextMenus',[hcm1,hcm2,hcm3])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pfp = figurePosition(fh,bh)
% Helper function to estimate the optimal position for the gui figure given
% the extent of the system monitor and the bounding box required after
% laying out the graph.
%
% fh if the figure handle
% bh is the biograph object handle

% find out the monitor extent 
me = get(0,'ScreenSize') + [4 33 -8 -106];

% maximum allowable size of figure (in points) given the monitor properties
msf = hgconvertunits(fh,me,'pixels','points',0);

% get the default MATLAB figure position (in points)
% so we preserve try to keep the left-up corner
dfp = hgconvertunits(fh,get(fh,'Position'),'pixels','points',0);

% required figure size after doing the layout (add 4 pts at the borders)
rfp = (bh.Scale*bh.BoundingBox)+8;


if any(rfp(3:4)>msf(3:4)) % need to scale down?
    scl= min(msf(3:4)./rfp(3:4));
elseif all(rfp(3:4)<dfp(3:4)) % need to scale up? 
    scl = min(dfp(3:4)./rfp(3:4));
else
    scl = 1;
end

% figure size after scaling
fsas = rfp*scl;
fsas = max([fsas;dfp]);

% final figure position after using default location and new size
pfp = [dfp(1) dfp(2)+dfp(4)-fsas(4) fsas(3:4)];

% now shift if it goes beyond the monitor limits
pfp(2) = pfp(2)-min(0,pfp(2)-msf(2));
pfp(1) = pfp(1)-max(0,pfp(3)-dfp(3))/2;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zoomListener(~,~,fh,bh)
% Helper function to update the rendering scale and change the fontsize
% when the user is Zooming. The function is triggered by a listener
% attached to 'Ylim' property of the axes. 
%
% fh if the figure handle
% bh is the biograph object handle

if ~isequal(class(bh),'biograph.biograph')
    return
end

% shortcut to leave if we arrived here when no Zooming
% ZoomOnState = getappdata(fh,'ZoomOnState')
% if isempty(ZoomOnState) || isequal(ZoomOnState,'off')
%     return
% end

% find handles to listeners
hl = getappdata(fh,'bgguiListeners'); 
hk = getappdata(fh,'biographListeners'); 

% turn 'off' listener so I do not trigger it again when changing the axes
% properties
%set(hl,'Enabled','off')
%set(hk,'Enabled','off')
bioinfoprivate.bioToggleListenerState(hl,'off')
bioinfoprivate.bioToggleListenerState(hk,'off')

% handle to axes
ha = bh.hgAxes;

% compute the a new scale; based on the required zoom and the actual yx
% limits of the axes. Make sure that we keep the same aspect ratio, so we
% do not loose the proportion of the nodes
haxlim = get(ha,'XLim'); haylim = get(ha,'YLim'); hapos = get(ha,'Position');
cf = (max([diff(haxlim) diff(haylim)]./hapos(3:4))*hapos(3:4)-[diff(haxlim) diff(haylim)])/2;
haxlim = haxlim  - [cf(1) -cf(1)];
haylim = haylim  - [cf(2) -cf(2)];
set(ha,'XLim',haxlim);set(ha,'YLim',haylim);
hapos = get(ha,'Position');
scale = hapos(3)/diff(haxlim);
oldScale = getappdata(ha,'Scale');
setappdata(ha,'Scale',scale)

if abs(oldScale - scale) > sqrt(eps)
  % adjusting font size based on the new scale
  bh.hgCorrectFontSize
end

% enable listener before leaving
%set(hl,'Enabled','on')
%set(hk,'Enabled','on')
bioinfoprivate.bioToggleListenerState(hl,'on')
bioinfoprivate.bioToggleListenerState(hk,'on')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bgguiResize(fh,~,bh)
% Helper function to update the rendering scale and change the fontsize
% when the user Resized the figure. The function is triggered by the
% callback 'ResizeFcn' in the figure properties.
%
% fh if the figure handle
% bh is the biograph object handle

if ~isequal(class(bh),'biograph.biograph')
    return
end

% find handle to zoom listener
lh = getappdata(fh,'bgguiListeners'); 
kh = getappdata(fh,'biographListeners'); 

% turn 'off' the zoom listener so I do not trigger it.
%set(lh,'Enabled','off')
%set(kh,'Enabled','off')
bioinfoprivate.bioToggleListenerState(lh,'off')
bioinfoprivate.bioToggleListenerState(kh,'off')

% update layout and then correct the size of texts
bh.hgUpdate

% enable listener
%set(lh,'Enabled','on')
%set(kh,'Enabled','on')
bioinfoprivate.bioToggleListenerState(lh,'on')
bioinfoprivate.bioToggleListenerState(kh,'on')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = getbiographtoolnumber()
% Computes the index number for this particular tool

% first, finds the used numbers so far
allFigs = findall(0,'tag','BioGraphTool');
usedNumbers = zeros(1,numel(allFigs)+1);
for i = 1:numel(allFigs)
    str = get(allFigs(i),'Name');
    % check that no-one  has chenged the mane to the figure window
    if isequal(str(1:min(16,end)),'Biograph Viewer ')
        tmp = str2num(str(17:end)); %#ok<ST2NM>
        if ~isempty(tmp)
            usedNumbers(i) = tmp;
        end
    end
end

% This is how we find the next index.  The rule is that we find the lowest
% integer value (non-zero and positive) not yet prescribed to a phytree
% tool, This is the same way MATLAB figures behave.
n = num2str(min(setdiff(1:(max(usedNumbers)+1),usedNumbers)));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeBiographViewerToolbar(fig,bh) %#ok<INUSD>
% helper function to set the toolbar
oldSH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')
set(fig,'toolbar','figure')  % needs to update because uicontrols turn it off

% Fix toolbar options, we keep: ZoomIn,ZoomOut,Pan
hw = findall(fig,'type','uitoolbar');
hf = get(hw,'Children');
set(hf,'Separator','off')
h1 = findall(hf,'Tag','Exploration.Pan');
h2 = findall(hf,'Tag','Exploration.ZoomOut');
h3 = findall(hf,'Tag','Exploration.ZoomIn');
delete(setxor(hf,[h1,h2,h3]))
set(0,'ShowHiddenHandles',oldSH)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeBiographViewerUIMenus(fig,bh)
% helper function to set UI menus
oldSH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

% delete figure menus not used
%h1 = findall(fig,'Type','uimenu', 'Label','&Edit');
h1 = findall(fig,'Type','uimenu', 'Tag','figMenuEdit');
%h2 = findall(fig,'Type','uimenu', 'Label','&View');
h2 = findall(fig,'Type','uimenu', 'Tag','figMenuView');
%h3 = findall(fig,'Type','uimenu', 'Label','&Insert');
h3 = findall(fig,'Type','uimenu', 'Tag','figMenuInsert');
%h4 = findall(fig,'Type','uimenu', 'Label','&Desktop');
h4 = findall(fig,'Type','uimenu', 'Tag','figMenuDesktop');
delete([h1,h2,h3,h4])

% Repair "File" menu
%hw = findall(fig,'Type','uimenu', 'Label','&File');
hw = findall(fig,'Type','uimenu', 'Tag','figMenuFile');
hf = get(hw,'children');
%h1 = findall(hw,'Label','Expo&rt Setup...');
h1 = findall(hw,'Tag','figMenuFileExportSetup');
%h3 = findall(hw,'Label','Print Pre&view...');
h3 = findall(hw,'Tag','figMenuFilePrintPreview');
%h4 = findall(hw,'Label','&Print...');
h4 = findall(hw,'Tag','printMenu');
delete(setxor(hf,[h1,h3,h4]))

uimenu(hw,'Label','Export to Workspace...','Position',1,'Callback',{@exportToWorkspace,bh})
uimenu(hw,'Label','Print to Figure','Position',2,'Callback',{@printToFigure,bh})
uimenu(hw,'Label','Exit','Separator','on','Position',6,'Callback','close(gcbf)')
set(h1,'Separator','on')
    
% Repair "Tools" menu
%hw = findall(fig,'Type','uimenu','Label','&Tools');
hw = findall(fig,'Type','uimenu','Tag','figMenuTools');
hf = get(hw,'children');
h1 = findall(hw,'Tag','figMenuZoomIn');    
h2 = findall(hw,'Tag','figMenuZoomOut');   
h3 = findall(hw,'Tag','figMenuPan');       
h4 = findall(hw,'Tag','figMenuResetView'); 
delete(setxor(hf,[h1,h2,h3,h4]))
set(h1,'Separator','off')
h5 = uimenu(hw,'Label','Layout Properties...','Separator','on','Position',5,'Callback',{@callPropertyEditor,bh}); %#ok<NASGU>
h6 = uimenu(hw,'Label','Refresh Layout','Position',6,'Callback',{@redoLayout,bh}); %#ok<NASGU>
h7 = uimenu(hw,'Label','Refresh Edges','Position',7,'Callback',{@redoLayout,bh}); %#ok<NASGU>

% Repair "Help" menu
%hw = findall(fig,'Type','uimenu','Label','&Help');
hw = findall(fig,'Type','uimenu','Tag','figMenuHelp');
delete(get(hw,'children'));
uimenu(hw,'Label','Bioinformatics Toolbox Help','Position',1,'Callback',...
       'helpview(fullfile(docroot,''toolbox'',''bioinfo'',''bioinfo.map''),''bioinfo_product_page'')')
uimenu(hw,'Label','Examples','Position',2,'Separator','on',...
       'Callback','demo(''toolbox'',''bioinfo'')')   
tlbx = ver('bioinfo');
mailstr = ['web(''mailto:bioinfofeedback@mathworks.com?subject=',...
           'Feedback%20for%20Biograph%20Viewer%20in%20Bioinformatics',...
           '%20Toolbox%20',tlbx(1).Version,''')'];
uimenu(hw,'Label','Send Feedback','Position',3,'Separator','on',...
       'Callback',mailstr);

set(0,'ShowHiddenHandles',oldSH)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hcm1,hcm2,hcm3] = makeBiographViewerContextMenus(fig,bh) 
% helper function to set context menus

hcm1 = uicontextmenu('Tag','Node','Parent',fig);
set(hcm1,'Callback',{@updateNodeContextMenus,bh})
h1 = uimenu(hcm1,'Label','Ancestors',  'Callback',{@findRelatives,bh});
h2 = uimenu(hcm1,'Label','Descendants','Callback',{@findRelatives,bh});
h3 = uimenu(hcm1,'Label','Both',  'Callback',{@findRelatives,bh});
h4 = uimenu(hcm1,'Label','Node Properties...', 'Separator','on',...
                                'Callback',{@callPropertyEditor,bh});
setappdata(hcm1,'fixedMenuItems',[h1,h2,h3,h4])                         

hcm2 = uicontextmenu('Tag','Edge','Parent',fig);  
set(hcm2,'Callback',{@updateEdgeContextMenus,bh})
h6 = uimenu(hcm2,'Label','Edge Properties...',...
                                'Callback',{@callPropertyEditor,bh});
setappdata(hcm2,'fixedMenuItems',h6)

hcm3 = uicontextmenu('Tag','Layout','Parent',fig);                                    
h8 = uimenu(hcm3,'Label','Layout Properties...',...
                                  'Callback',{@callPropertyEditor,bh}); %#ok<NASGU>
h9 = uimenu(hcm3,'Label','Refresh Layout','Callback',{@redoLayout,bh}); %#ok<NASGU>
h10 = uimenu(hcm3,'Label','Refresh Edges','Callback',{@redoLayout,bh}); %#ok<NASGU>

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateNodeContextMenus(hcm1,~,bh)
setappdata(gcbo,'node',getappdata(gco,'node'))
delete(setxor(get(hcm1,'Children'),getappdata(hcm1,'fixedMenuItems')))
if iscell(bh.NodeCallbacks)
    for i = 1:numel(bh.NodeCallbacks)
        uimenu(hcm1,'Label',sprintf('User Callback %d',i),...
               'Callback',{@callUserCallback,bh,'Node',i});
    end
else
    uimenu(hcm1,'Label','User Callback',...
               'Callback',{@callUserCallback,bh,'Node',0});
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateEdgeContextMenus(hcm2,~,bh)
setappdata(gcbo,'edge',getappdata(gco,'edge'))
delete(setxor(get(hcm2,'Children'),getappdata(hcm2,'fixedMenuItems')))
if iscell(bh.EdgeCallbacks)
    for i = 1:numel(bh.EdgeCallbacks)
        uimenu(hcm2,'Label',sprintf('User Callback %d',i),...
               'Callback',{@callUserCallback,bh,'Edge',i});
    end
else
    uimenu(hcm2,'Label','User Callback',...
               'Callback',{@callUserCallback,bh,'Edge',0});
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function callUserCallback(h,~,bh,caller,multicallbks)

switch caller
    case 'Node'
        node = getappdata(get(h,'Parent'),'node');
        if multicallbks
           cbck = bh.nodeCallback{multicallbks};
        else
           if iscell(bh.nodeCallback)
           cbck = bh.nodeCallback{1};
           else
               cbck = bh.nodeCallback;
           end
        end
        if ischar(cbck)
            eval([cbck '(node);'])
        else
            cbck(node)
        end
    case 'Edge'
        edge = getappdata(get(h,'Parent'),'edge');
        if multicallbks
           cbck = bh.edgeCallback{multicallbks};
        else
           if iscell(bh.edgeCallback)
               cbck = bh.edgeCallback{1};
           else
               cbck = bh.edgeCallback;
           end
        end
        if ischar(cbck)
            eval([cbck '(edge);'])
        else
            cbck(edge)
        end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function callPropertyEditor(h,~,bh)

switch get(h,'Label')
    case 'Node Properties...'
        node = getappdata(get(h,'Parent'),'node');
        inspect(node)
    case 'Edge Properties...'
        edge = getappdata(get(h,'Parent'),'edge');
        inspect(edge)
    case 'Layout Properties...'
        inspect(bh)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function findRelatives(h,~,bh)

% find selected Nodes
selNodes = bh.nodes(mycell2mat(get(bh.nodes,'isSelected')));

% get the calling node
node = getappdata(get(h,'Parent'),'node');

if ~ismember(node,selNodes)
    if ~isempty(selNodes)
       set(selNodes,'isSelected',false)
       selNodes.hgUpdate
    end
    selNodes = node;
    selNodes.isSelected = true;
    selNodes.hgUpdate
end

switch get(h,'Label')
    case 'Ancestors'
        selNodes = setdiff(getancestors(selNodes),selNodes);
    case 'Descendants'
        selNodes = setdiff(getdescendants(selNodes),selNodes);
    case 'Both'
        selNodes = setdiff(getrelatives(selNodes),selNodes);
end

if ~isempty(selNodes)
    set(selNodes,'isSelected',true)
    selNodes.hgUpdate
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function redoLayout(h,~,bh)
set(gcbf,'Pointer','watch')
switch get(h,'Label')
   case 'Refresh Layout'
       bh.dolayout
   case 'Refresh Edges'
       bh.dolayout('pathsonly',true)
end
bh.hgUpdate
bh.hgReDraw
set(gcbf,'Pointer','arrow')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function exportToWorkspace(~,~,bh)
global private_pointer_to_BIOGRAPH
s = inputdlg('Workspace variable name ?','Export to Workspace',1);
while ~(isempty(s) || isvarname(s{1}) || isempty(s{1}))
      s = inputdlg('Not a valid variable name, type a MATLAB variable name ?','Export to Workspace',1);
end
if ~(isempty(s) || isempty(s{1}))
    evalin('base','global private_pointer_to_BIOGRAPH');
    private_pointer_to_BIOGRAPH = bh;
    evalin('base',[s{1} '= private_pointer_to_BIOGRAPH.deepCopy;'])
    evalin('base','clear private_pointer_to_BIOGRAPH')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printToFigure(~,~,bh)
g = figure('Units','Points','Position',get(get(bh.hgAxes,'Parent'),'Position')+[20 -20 0 0]);
copyobj(bh.hgAxes,g);
set(g,'Units','Pixels');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseClickOnFigure(fh,~,bh,gh)
% helper for mouse click

% in case the contained BIOGRAPH was destroyed this control is deactivated
if ~isequal(class(bh),'biograph.biograph')
    return
end

[node,cp] = getNodeAtMousePosition(bh);
innode = ~isempty(node);

if innode && node.isSelected %start moving nodes
    switch get(fh,'SelectionType')
        case {'normal','extend'}
            gh.moving = true;
            gh.innode = true;
            gh.currentNode = node;
            gh.movNodesPos = cp; 
            selNodes = find(mycell2mat(get(bh.nodes,'isSelected')));
            patchLengths = zeros(numel(selNodes),1);
            for k=1:numel(selNodes)
                patchLengths(k) = numel(get(bh.Nodes(selNodes(k)).hgPatch,'Xdata'));
            end
            patchXData = zeros(sum(patchLengths)+numel(selNodes),1);
            patchYData = zeros(sum(patchLengths)+numel(selNodes),1);
            j = 0;
            for k=1:numel(selNodes)
                v = j+(1:patchLengths(k));
                patchXData(v) = get(bh.Nodes(selNodes(k)).hgPatch,'Xdata');
                patchYData(v) = get(bh.Nodes(selNodes(k)).hgPatch,'Ydata');
                j = j + patchLengths(k)+1;
                patchXData(j) = patchYData(1);
                patchYData(j) = patchYData(1);    
            end    
            gh.movPatchXData = patchXData;   
            gh.movPatchYData = patchYData;
            set(gh.dataTip,'Visible','off');
        case 'open'
            gh.moving = false;
            gh.opening = true;
            gh.innode = true;
            gh.currentNode = node;
        case 'alt'
            % alt key unavailable
    end %switch
else %start selecting nodes
    switch get(fh,'SelectionType')
        case {'normal','extend'}
            gh.selecting = true;
            gh.innode = innode;
            gh.currentNode = node;
            set(gh.dragBox,'XData',repmat(cp(1),5,1));
            set(gh.dragBox,'YData',repmat(cp(3),5,1));
            set(gh.dataTip,'Visible','off');
        case 'open'
            gh.selecting = false;
            gh.opening = true;
            gh.innode = innode;
            gh.currentNode = node;
        case 'alt'
            % alt key unavailable
    end %switch
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localWindowButtonMotion(fh,event,bh,gh) %#ok
% Callback function activated when moving over the axes, checks location of
% the mouse and puts datatip if over an active element.

% in case the contained BIOGRAPH was destroyed this control is deactivated
if ~isequal(class(bh),'biograph.biograph')
    return
end

if ~gh.selecting && ~ gh.moving   % no click, just show data tips if needed
    [sn,cp] = getNodeAtMousePosition(bh);
    if isempty(sn)
        set(gh.dataTip,'Visible','off');
    else
        strc{1} = sprintf('ID: %s',sn.ID);
        if ~isempty(sn.label)
            strc{2} = sn.label;
        end
        if ~isempty(sn.description)
            strc{end+1} = sn.description;
        end
        str = char(strc);
        axs = [get(bh.hgAxes,'Xlim');get(bh.hgAxes,'Ylim')];
        set(gh.dataTip,'String',str);
        extmp = get(gh.dataTip,'Extent');
        ex = extmp(3:4)';
        p = cp([1;3])+diff(axs,[],2)/100;
        set(gh.dataTip,'Position',[p' 1]);
        extmp = get(gh.dataTip,'Extent'); 
        t = sum(reshape(extmp,2,2),2)>max(axs,[],2);
        if any(t)
            if t(2)
                p = cp([1;3]) - diff(axs,[],2)/100 - ex;
            else
                p(1) = cp(1) - diff(axs(1,:))/100 - ex(1);
            end
            set(gh.dataTip,'Position',[p' 1]);
        end
        set(gh.dataTip,'Visible','on');
    end
elseif gh.selecting               % dragging the selecting box
    gh.dragging = true;
    cp = get(bh.hgAxes,'CurrentPoint');
    Xda = get(gh.dragBox,'XData');
    Yda = get(gh.dragBox,'YData');
    Xda([3,4]) = cp(1);
    Yda([2,3]) = cp(3);
    set(gh.dragBox,'XData',Xda)
    set(gh.dragBox,'YData',Yda)
    set(gh.dragBox,'Visible','on')
elseif gh.moving                  % dragging to move selected nodes
    gh.dragging = true;  
    cp = get(bh.hgAxes,'CurrentPoint');
    set(gh.ghostNodes,'Xdata',gh.movPatchXData + cp(1) - gh.movNodesPos(1));
    set(gh.ghostNodes,'Ydata',gh.movPatchYData + cp(3) - gh.movNodesPos(3));
    set(gh.ghostNodes,'Visible','on');
else
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseRelease(fh,event,bh,gh) %#ok
% helper for mouse release


% in case the contained BIOGRAPH was destroyed this control is deactivated
if ~isequal(class(bh),'biograph.biograph')
    return
end

% deactivate property listeners
%set(getappdata(fh,'biographListeneres'),'enable','off')
bioinfoprivate.bioToggleListenerState(getappdata(fh,'biographListeneres'),'off')


if ~gh.dragging && gh.innode && ~gh.opening % select/unselect without dragging
    switch get(fh,'SelectionType')
        case 'normal'
            previousState = gh.currentNode.isSelected;
            selNodes = find(mycell2mat(get(bh.nodes,'isSelected')));
            if ~isempty(selNodes)
                set(bh.Nodes(selNodes),'isSelected',false)
            end
            gh.currentNode.isSelected = ~previousState;
            nodesToUpdate =[selNodes;gh.currentNode.idx];
            bh.Nodes(nodesToUpdate).hgUpdate
        case 'extend'    
            gh.currentNode.isSelected = ~gh.currentNode.isSelected;
            gh.currentNode.hgUpdate
    end
elseif gh.dragging && gh.moving && ~gh.opening % moving nodes
    selNodes = find(mycell2mat(get(bh.nodes,'isSelected')));
    cp = get(bh.hgAxes,'CurrentPoint');
    vector = cp([1,3])-gh.movNodesPos([1,3]);
    vector = vector/bh.Scale;
    if ~isempty(selNodes)
        for k = 1:numel(selNodes)
            bh.Nodes(selNodes(k)).Position = ...
                           bh.Nodes(selNodes(k)).Position + vector;
        end
        bh.Nodes(selNodes).hgUpdate
        
    end
    edges = bh.edges(union(nonzeros(bh.to(:,selNodes)),nonzeros(bh.from(:,selNodes))));
    if ~isempty(edges)
        edges = unique(edges);
        for k = 1:numel(edges)
            if edges(k).ToNode.isSelected && edges(k).FromNode.isSelected
                edges(k).ControlPoints = edges(k).ControlPoints + ...
                              repmat(vector',1,size(edges(k).ControlPoints,2));
            elseif edges(k).ToNode.isSelected
                segm = sqrt(sum(diff(edges(k).ControlPoints,[],2).^2));
                if ~any(segm) 
                    segm = (1:numel(segm))/numel(segm); 
                end
                edges(k).ControlPoints = edges(k).ControlPoints + ...
                                         vector'*[0 cumsum(segm)]/sum(segm);
            else
                segm = sqrt(sum(diff(edges(k).ControlPoints,[],2).^2));
                if ~any(segm)
                    segm = (1:numel(segm))/numel(segm);
                end
                edges(k).ControlPoints = edges(k).ControlPoints + ...
                           vector'*[fliplr(cumsum(segm(end:-1:1))) 0]/sum(segm);
            end
        end
    edges.hgUpdate
    end
    set(gh.ghostNodes,'Visible','off');
elseif gh.dragging && gh.selecting % selecting nodes with dragging box
    if numel(bh.Nodes)==1
        nodePos = get(bh.Nodes,'Position')*bh.Scale;
        nodeVis = get(bh.Nodes,'Visible');
    else
        nodePos = cell2mat(get(bh.Nodes,'Position'))*bh.Scale;
        nodeVis = cell2mat(get(bh.Nodes,'Visible'));
    end
    nodesPicked = find(nodePos(:,1)>min(get(gh.dragBox,'XData')) & ...
                       nodePos(:,1)<max(get(gh.dragBox,'XData')) & ...
                       nodePos(:,2)>min(get(gh.dragBox,'YData')) & ...
                       nodePos(:,2)<max(get(gh.dragBox,'YData')) & ...
                       nodeVis);
    switch get(fh,'SelectionType')
        case 'normal'
            selNodes = find(mycell2mat(get(bh.nodes,'isSelected')));
            if ~isempty(selNodes)
                set(bh.Nodes(selNodes),'isSelected',false)
            end
            if ~isempty(nodesPicked)
                set(bh.Nodes(nodesPicked),'isSelected',true)
            end
            nodesToUpdate =[selNodes;nodesPicked];
            if ~isempty(nodesToUpdate)
                bh.Nodes(nodesToUpdate).hgUpdate
            end
        case 'extend'
            if ~isempty(nodesPicked)
                set(bh.Nodes(nodesPicked),'isSelected',true)
                bh.Nodes(nodesPicked).hgUpdate
            end
    end
    set(gh.dragBox,'Visible','off');
elseif ~gh.dragging && gh.selecting % single click on axes
     selNodes = find(mycell2mat(get(bh.nodes,'isSelected')));
     if ~isempty(selNodes)
         set(bh.Nodes(selNodes),'isSelected',false)
         bh.Nodes(selNodes).hgUpdate
     end
elseif gh.opening && gh.innode % opening a node executes the user nodeCallback
     if iscell(bh.nodeCallback)
        cbck = bh.nodeCallback{1};
     else 
        cbck = bh.nodeCallback;
     end
     if ischar(cbck)
        eval([cbck '(gh.currentNode);'])
    else 
        cbck(gh.currentNode)
    end
elseif gh.opening && ~gh.innode  % double click on axes
    bh.hgUpdate
    %bh.hgReDraw  %expensive, only do it if necessary  
elseif ~gh.dragging && ~gh.innode && ~gh.opening && ~gh.moving && ~gh.selecting   
    % no action 
else
    % no action - error ?
end

% done with action, resetting status to inactive
gh.selecting = false; gh.dragging = false; gh.opening = false;
gh.moving = false;    gh.innode = false;
 
%set(getappdata(fh,'biographListeneres'),'enable','on')
bioinfoprivate.bioToggleListenerState(getappdata(fh,'biographListeneres'),'on')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [node,cp] = getNodeAtMousePosition(bh)
% Get the node ID at the current mouse position (node must be visible)
cp = get(bh.hgAxes,'CurrentPoint');
node = [];
if numel(bh.Nodes)==1
    pos = cp([1 3]) - bh.Scale*get(bh.Nodes,'Position');
    pos(~get(bh.Nodes,'Visible'),:) = inf;
else
    pos = bh.Scale*cell2mat(get(bh.Nodes,'Position'));
    pos(:,1) = cp(1)-pos(:,1);
    pos(:,2) = cp(3)-pos(:,2);
    pos(~cell2mat(get(bh.Nodes,'Visible')),:) = inf;
end
[pos,h] = sort(sum(pos.^2,2));
for k = 1:min(3,numel(h)) %look within the three closest nodes
    if bh.nodes(h(k)).visible && inpolygon(cp(1),cp(3),...
                 get(bh.nodes(h(k)).hgPatch,'xdata'),...
                 get(bh.nodes(h(k)).hgPatch,'ydata'))
        node = bh.nodes(h(k)); % got this node
        break % if found, break
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deletebggui(hSrc,event)  %#ok<INUSD>
% delete figure if object has been deleted, as HG does

if ishandle(hSrc) && ~isempty(hSrc.up) && ~isempty(hSrc.up.hgFigure) ...
   && ishandle(hSrc.up.hgFigure)
       delete(hSrc.up.hgFigure);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = mycell2mat(c)
if numel(c)==1
    m = c;
else
    m = cell2mat(c);
end
