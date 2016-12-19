function view(tr,sel,propsForFigure)
%VIEW views a phylogenetic tree in phytreeviewer.
%
%   VIEW(TREE) shows a phylogenetic tree object. The significant distances
%   between branches and nodes are in horizontal direction, vertical
%   coordinates are accommodated only for display purposes. Tree
%   Edit/Analysis tools are accessible through the mouse left/right buttons
%   and also using the 'Tree' menu.  
%
%   VIEW(TREE,SEL) starts the viewer with an initial selection of nodes
%   specified by SEL. SEL can be a logical array of any of the following
%   sizes: [NUMLEAVES+NUMBRANCHES x 1], [NUMLEAVES x 1], or [NUMBRANCHES x
%   1]. SEL may also be a list of indices. 
% 
%   Examples:
%      
%       tr = phytreeread('pf00002.tree')
%       view(tr)
%
%   See also PHYTREE, PHYTREE/PLOT, PHYTREEREAD, PHYTREEVIEWER, SEQLINKAGE,
%   SEQNEIGHJOIN.

% Copyright 2003-2012 The MathWorks, Inc.


if numel(tr)~=1
     error(message('bioinfo:phytree:view:NoMultielementArrays'));
end

tr = doBasicCalculations(tr);

nodeIndex   = 1:tr.numLabels;
leafIndex   = 1:tr.numLeaves;
branchIndex = tr.numLeaves+1:tr.numLabels;

% check empty names
for ind = nodeIndex
    if isempty(tr.names{ind}) 
        if ind > tr.numLeaves
            tr.names{ind} = ['Branch ' num2str(ind-tr.numLeaves)]; 
        else
            tr.names{ind} = ['Leaf ' num2str(ind)]; 
        end
    end
end

% initial drawing
if nargin<3
    propsForFigure.Name = ['Phylogenetic Tree ' getphytreetoolnumber];
end
propsForFigure.PruneWarning = getacceptedwarningfromothertools;
fig = figure('Renderer','ZBuffer','Name',propsForFigure.Name,...
           'NumberTitle','off','IntegerHandle','off','tag','PhyTreeTool');
setappdata(fig,'propsForFigure',propsForFigure)  
setappdata(fig,'backupTree',tr) 
tr.ha = axes; hold on;
set(tr.ha,'Position',[.05 .05 .7 .9],'YTick',leafIndex,'FontSize',9,'Ydir','reverse',...
          'YAxisLocation','Right','YTickLabel',char(tr.names{leafIndex}))
tr.hlines = plot( ...
  tr.x([nodeIndex;repmat([tr.par(1:tr.numLabels-1) tr.numLabels],2,1)]),...
  tr.y([repmat(nodeIndex,2,1);[tr.par(1:tr.numLabels-1) tr.numLabels]]),...
                '-k');
tr.hpathline = plot(1,1,'--r','LineWidth',2,'Visible','off');  
tr.hdragbox = plot(1,1,':k','LineWidth',1,'Visible','off'); 
tr.hdots(1,1) = plot(tr.x(branchIndex),tr.y(branchIndex),'o',...
               'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','b');
tr.hdots(1,2) = plot(tr.x(leafIndex),tr.y(leafIndex),'square',...
               'MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor','w');
tr.hseldots(1,1) = plot(tr.x(branchIndex),tr.y(branchIndex),'o',...
               'MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r');
tr.hseldots(1,2) = plot(tr.x(leafIndex),tr.y(leafIndex),'square',...
               'MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r');
tr.hldots(1,1) = plot(tr.x(branchIndex),tr.y(branchIndex),'o',...
               'MarkerSize',5,'MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor',[.6 .6 1]);
tr.hldots(1,2) = plot(tr.x(leafIndex),tr.y(leafIndex),'square',...
               'MarkerSize',4,'MarkerEdgeColor',[.5 .5 .5],...
               'MarkerFaceColor','w');     
set(tr.hldots(1),'Xdata',[],'Ydata',[])
set(tr.hldots(2),'Xdata',[],'Ydata',[])
tr.axhold = plot([-eps -eps],[0 0],'.','MarkerSize',eps,'Color','w');
tr.datatip = text(0,1,1,'k','Tag','TreeTag','BackgroundColor',[1 1 .93],...
               'Color', [0 0 0],'EdgeColor', [0.8 0.8 0.8],...
               'VerticalAlignment','Top','Clipping','off',...
               'Visible','off','Fontsize',8,'Interpreter','none');         

if nargin == 1 || isempty(sel)
   tr.selected = false(tr.numLabels,1);        % selected nodes
else
    % validate sel
    if islogical(sel)
        if numel(sel)==tr.numLabels 
            sel = sel(:)==true;
        elseif numel(sel)==tr.numLeaves
            sel = [sel(:);false(tr.numBranches,1)];
        elseif numel(sel)==tr.numBranches
            sel = [false(tr.numLeaves,1);sel(:)];
        else
        close(fig)    
        error(message('bioinfo:phytree:view:IncorrectLogical'));
        end
    elseif isnumeric(sel) && isreal(sel) && all(sel>=1) && all(sel<=tr.numLabels)
        tem(tr.numLabels)=false;
        tem(floor(sel))=true;
        sel=tem(:);
    else
        close(fig)    
        error(message('bioinfo:phytree:view:IncorrectTypeofArguments'));
    end
    tr.selected =sel;
end              
               
% save more figure data needed for the gui functionality
tr.activeNodes    = true(tr.numLabels,1);   % active nodes
tr.activeBranches = true(tr.numBranches,1); % active Branches
tr.sel2root = false(tr.numLabels,1);        % path sel-node to root
tr.editMode = 'Select';                            % initial edit mode
tr.indicativeMode = false;                         % data-tip flag
tr.lastThresholdValue = [];                        % remembers last cut

% create uicontrols (will appear as needed, initially invisible)
tr.editBox =  uicontrol(fig,'Background',[1 1 1],'style','edit',...
                        'visible','off','callback',@doneRenaming);                     
tr.slider  =  uicontrol(fig,'style','slider','SliderStep',[.1 .1],...
                        'visible','off','callback',@sliderCallback);
tr.slidertx = uicontrol(fig,'style','text','visible','off');
tr.sliderok = uicontrol(fig,'style','pushbutton','visible','off',...
                        'string','OK','callback',@doThresholdCut);

tr.slider2  =  uicontrol(fig,'style','slider','SliderStep',[.1 .1],...
                        'visible','off','callback',@slider2Callback);
tr.slider2tx = uicontrol(fig,'style','edit','visible','off',...
                        'Background',[1 1 1],'callback',@slider2txCallback);
tr.slider2ok = uicontrol(fig,'style','pushbutton','visible','off',...
                        'string','OK','callback',@doSelectTool);
                    
                    
% setup callback for click over nodes                    
set([tr.hseldots,tr.hdots,tr.hldots],'ButtonDownFcn',@toggleNode)
% setup figure callback functions
set(fig,'WindowButtonDownFcn',@mouseClickOnFigure);
set(fig,'WindowButtonUpFcn',@mouseRelease);
set(fig,'WindowButtonMotionFcn',@localWindowButtonMotion);

% setup UIMenus, context menus and toolbar
tr.hToggleUIMenu  = makePhyTreeViewerUIMenus(fig);                    
tr.hToggleToolbar = makePhyTreeViewerToolbar(fig);
[tr.hToggleContextMenu,tr.hAxisContextMenu,tr.hDotsContextMenu] = ...
                                        makePhyTreeViewerContextMenus(fig);
% activate Context Menus
set(tr.ha,'UIContextMenu',tr.hAxisContextMenu);
set([tr.hdots tr.hldots tr.hseldots],'UIContextMenu',tr.hDotsContextMenu);

% figures out the label extents once to speedup the automatic shifting of
% the right edge of the axis
tr.anytext = text(1,1,repmat('A',1,1000),'Visible','off','fontsize',9,'units','pixels','interpreter','none');
tr.labelExtents = zeros(tr.numLabels,1);
for i = 1:tr.numLabels
   set(tr.anytext,'string',tr.names{i})    
   tr.labelExtents(i) = get(tr.anytext,'Extent')*[0;0;1;0];
end
tr.maxLabelExtent = max(tr.labelExtents(1:tr.numLeaves));

set(fig,'UserData',tr)           % save figure data

correctFigureSize(fig, 15 * tr.numLeaves);        % resize figure if needed
setupYLabelsListeners(fig,tr.ha);           % listeners for YLabels
updateTree(fig,[])            % updates figure after all initializations 
set(tr.ha,'Ylim',[min(tr.y(tr.activeNodes))-1,max(tr.y(tr.activeNodes))+1]);
set(tr.ha,'Xlim',[0  max(tr.x)] + max(tr.x) * [-.1 .05]);  
tr.yLim = get(tr.ha,'Ylim');
tr.xLim = get(tr.ha,'Xlim');
set(fig,'UserData',tr)           % save figure data
toolsmenufcn(fig,'PanY')         % set zoom mode to vertical constraining
toolsmenufcn(fig,'ZoomY')        % set pan  mode to vertical constraining
set(fig,'HandleVisibility','callback','Visible','on') % after all init, 
%                               hide the handles and force to the front

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ylabelsListener(~,~,hf,ha) 
% Auto sizes the ylabels
ratio = max(get(hf,'Position').*[0 0 0 1])/diff(get(ha,'YLim'));
set(ha,'Fontsize',max(2,min(9,ceil(ratio/1.7))));    % the gold formula

tr=get(hf,'Userdata');
wS = get(hf,'Position');  % window dimensions
aP = get(ha,'Position'); % axes position

% check if right edge of the axis needs to be moved
ff =[1 2 3 4 5 5 6 7 8 9 9 11];
fe = get(hf,'Position')*[0;0;1;0]; % figure extent
le = ff(get(ha,'Fontsize'))/7*tr.maxLabelExtent; % label extent (in pixels)
le = min(max(le,fe*0.04),fe*0.65);
aP(3) = 0.94-le/fe;             % figure out the new axis extent (normalized)
set(ha,'position',aP);

% Also verify if we need to re-position the slidebar of threshold cut
if any(strcmp(tr.editMode,{'Distance to Leaves','Distance to Root'}))  
    set(tr.slider,  'Position',[aP(1)*wS(3) wS(4)-30 max(1,aP(3)*wS(3)-80) 20])
    set(tr.slidertx,'Position',[sum(aP([1 3]))*wS(3)-80 wS(4)-30 60 20])
    set(tr.sliderok,'Position',[sum(aP([1 3]))*wS(3)-20 wS(4)-30 30 20])
end
% Also verify if we need to re-position the slidebar of threshold cut
if strcmp(tr.editMode,'Select Tool')  
    set(tr.slider2,  'Position',[aP(1)*wS(3) wS(4)-30 max(1,aP(3)*wS(3)-80) 20])
    set(tr.slider2tx,'Position',[sum(aP([1 3]))*wS(3)-80 wS(4)-30 60 20])
    set(tr.slider2ok,'Position',[sum(aP([1 3]))*wS(3)-20 wS(4)-30 30 20])
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseClickOnFigure(h,varargin)
% This callback function is activated when a mouse button is pressed in any
% location of the figure and under any of my edit modes
tr = get(gcbf,'Userdata');
switch tr.editMode
    case 'Renaming';           doneRenaming(h,varargin);
    case 'Distance to Leaves'; cancelThresholdCut(h,varargin);
    case 'Distance to Root';   cancelThresholdCut(h,varargin);
    case 'Select Tool';        cancelSelectTool(h,varargin);
    case 'Select';       
        switch get(gcbf,'SelectionType')
            case {'normal','extend'}
                tr = get(gcbf,'userdata');
                cp = get(tr.ha,'CurrentPoint');
                xPos = cp(1,1); yPos = cp(1,2); 
                set(tr.hdragbox,'Visible','on',...
                    'Xdata',repmat(xPos,5,1),'Ydata',repmat(yPos,5,1))
            case 'open'
                autoFit
        end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggleNode(h,varargin)
% This callback function is activated when a mouse button is pressed over
% any of the displayed nodes under any of my edit modes
hideActiveIndicators(h,varargin)
tr = get(gcbf,'Userdata');
switch get(gcbf,'SelectionType')
    case 'normal'
        switch tr.editMode
            case 'Select';             selectNode(h,varargin);
            case 'Inspect';            inspectNode(h,varargin);
            case 'Collapse/Expand';    collapseExpand(h,varargin);
            case 'Rotate Branch';      rotateBranch(h,varargin);
            case 'Rename';             renameNode(h,varargin);
            case 'Renaming';           doneRenaming(h,varargin);
            case 'Prune';              pruneTree(h,varargin);
            case 'Distance to Leaves'; cancelThresholdCut(h,varargin);    
            case 'Distance to Root';   cancelThresholdCut(h,varargin);
            case 'Select Tool';        cancelSelectTool(h,varargin);
        end
    case 'extend'
        switch tr.editMode
            case 'Select';             selectNode(h,varargin);
        end
    case 'alt'
    case 'open'    
 end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeEditMode(h,varargin)
% Callback function to change the edit mode, this function is
% called from the toolbar, the context menu or the uimenu.
tr = get(gcbf,'Userdata');
myModes = {'Inspect','Collapse/Expand','Rotate Branch','Rename','Prune'};

% first, disable any present edit mode
switch tr.editMode
    case [myModes,{'Select'}];  
                      disableMyContextMenus(h)
                      disableMyWindowButtonActions(h)
                      ind = strmatch(tr.editMode,myModes);
                      set(tr.hToggleToolbar(ind),    'State','off')
                      set(tr.hToggleUIMenu(ind),     'Checked','off')
                      %set(tr.hToggleContextMenu(ind),'Checked','off')
    %case '&Zoom In';  toolsmenufcn(gcbf,'ZoomIn');
    case '&Zoom In';  zoom(gcbf,'off')
    %case 'Zoom &Out'; toolsmenufcn(gcbf,'ZoomOut');
    case 'Zoom &Out'; zoom(gcbf,'off')
    %case '&Pan';      toolsmenufcn(gcbf,'Pan');
    case '&Pan';      pan(gcbf,'off');
    case {'Distance to Leaves','Distance to Root'}
                      enableAllUI(h)
                      disableMyWindowButtonActions(h)
    case 'Select Tool'
                      enableAllUI(h)
                      disableMyWindowButtonActions(h)                      
end
 
% depending on the caller instance, determine the new edit mode
switch get(h,'Type')
    case 'uimenu'; newEditMode = get(h,'Label');
    case 'uitoggletool'
        newEditMode = get(h,'Tag');
        switch newEditMode
            case 'Exploration.ZoomIn';  newEditMode = '&Zoom In';
            case 'Exploration.ZoomOut'; newEditMode = 'Zoom &Out';
            case 'Exploration.Pan';     newEditMode = '&Pan';
        end
    otherwise; newEditMode = 'Select';
end
%disp( [tr.editMode ' -->  ' newEditMode]  )
% if new mode is the same then we are toggling off
if strcmp(newEditMode,tr.editMode) 
    newEditMode = 'Select'; 
end

% if changing to Prune, verify the warnign has been accepted
if strcmp(newEditMode,'Prune') 
    propsForFigure = getappdata(gcbf,'propsForFigure');
    if isequal(propsForFigure.PruneWarning,'NotDone')
        warndlg(['Pruning nodes cannot be undone. Before continuing,',...
             ' you may want to export the current tree to a new viewer.'],...
             'Warning','modal')
        setacceptedwarningtoothertools
    end
end

switch newEditMode
    case '&Zoom In';  toolsmenufcn(gcbf,'ZoomIn'); 
    case 'Zoom &Out'; toolsmenufcn(gcbf,'ZoomOut');
    %case '&Pan';      toolsmenufcn(gcbf,'Pan'); 
    case '&Pan';      pan(gcbf,'on');
    case myModes;     enableMyContextMenus(h)
                      enableMyWindowButtonActions(h)
                      ind = strmatch(newEditMode,myModes);
                      set(tr.hToggleToolbar(ind),    'State','on')
                      set(tr.hToggleUIMenu(ind),     'Checked','on')
                      %set(tr.hToggleContextMenu(ind),'Checked','on')
    case 'Select';    enableMyContextMenus(h)
                      enableMyWindowButtonActions(h)
    case {'Distance to Leaves','Distance to Root'}
                      disableAllUI(h)
                      enableMyWindowButtonActions(h)
                      set(gcbf,'WindowButtonMotionFcn',[])
                      set(gcbf,'WindowButtonUpFcn',[])
                      set([tr.hseldots,tr.hdots,tr.hldots],'ButtonDownFcn',[])
    case 'Select Tool'
                      disableAllUI(h)
                      enableMyWindowButtonActions(h)
                      set(gcbf,'WindowButtonMotionFcn',[])
                      set(gcbf,'WindowButtonUpFcn',[])
                      set([tr.hseldots,tr.hdots,tr.hldots],'ButtonDownFcn',[])                      
end

switch newEditMode
    case 'Inspect';  if sum(tr.selected(:)) ~= 1
                         tr.selected(:) = false;
                         tr.selected(end) = true;
                     end
    otherwise
end    

tr.sel2root = path2root(tr, tr.selected);
tr.editMode = newEditMode;
set(gcbf,'userdata',tr)
updateTree(gcbf,[])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hideActiveIndicators(~,varargin) 
tr = get(gcbf,'userdata');
set([tr.hpathline,tr.datatip],'visible','off')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disableAllUI(~,varargin) 
hw = findall(gcbf,'Type','uimenu','Parent',gcbf);
gw = findall(gcbf,'Type','UIToggleTool');
set([hw;gw],'Enable','off')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enableAllUI(~,varargin) 
hw = findall(gcbf,'Type','uimenu','Parent',gcbf);
gw = findall(gcbf,'Type','UIToggleTool');
set([hw;gw],'Enable','on')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disableMyContextMenus(~,varargin) 
tr = get(gcbf,'userdata');
set(tr.ha,'UIContextMenu',[]);
set([tr.hdots tr.hseldots],'UIContextMenu',[]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enableMyContextMenus(~,varargin) 
tr = get(gcbf,'userdata');
set(tr.ha,'UIContextMenu',tr.hAxisContextMenu);
set([tr.hdots tr.hldots tr.hseldots],'UIContextMenu',tr.hDotsContextMenu);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disableMyWindowButtonActions(~,varargin) 
set(gcbf,'WindowButtonDownFcn',[]);
set(gcbf,'WindowButtonUpFcn',[])
set(gcbf,'WindowButtonMotionFcn',[]);
tr = get(gcbf,'userdata');
set([tr.hseldots,tr.hdots,tr.hldots],'ButtonDownFcn',[])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enableMyWindowButtonActions(~,varargin) 
set(gcbf,'WindowButtonDownFcn',@mouseClickOnFigure);
set(gcbf,'WindowButtonUpFcn',@mouseRelease);
set(gcbf,'WindowButtonMotionFcn',@localWindowButtonMotion);
tr = get(gcbf,'userdata');
set([tr.hseldots,tr.hdots,tr.hldots],'ButtonDownFcn',@toggleNode)
set(gcbf,'KeyPressFcn',[]);
set(gcbf,'Pointer','arrow');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localWindowButtonMotion(h,varargin) 
% Callback function activated when moving over the axes, checks location of
% the mouse and puts datatip if over an active node.

tr = get(h,'userdata');
% set a virtual grid to get the point
xThres=diff(get(tr.ha,'Xlim'))/100;
yThres=diff(get(tr.ha,'Ylim'))/100;
cp = get(tr.ha,'CurrentPoint');
xPos = cp(1,1); yPos = cp(1,2);
hp = tr.x<(xPos+xThres) & tr.x>(xPos-xThres) & ...
     tr.y<(yPos+yThres) & tr.y>(yPos-yThres);
hp = find (hp & tr.activeNodes);

% shortcut out when dragging a box in select mode
if strcmp(get(tr.hdragbox,'Visible'),'on')
    xdata = get(tr.hdragbox,'XData');xdata([3,4]) = xPos;
    ydata = get(tr.hdragbox,'YData');ydata([2,3]) = yPos;
    set(tr.hdragbox,'XData',xdata,'YData',ydata)
% shortcut out when turning off 'indicative' mode
elseif tr.indicativeMode && isempty(hp) %&& isempty(tr.highligth)
    set([tr.datatip tr.hpathline],'visible','off')
    set(tr.hlines,'color','black')
    set(tr.hldots(1),'Xdata',[],'Ydata',[])
    set(tr.hldots(2),'Xdata',[],'Ydata',[])
% turn on or update 'indicative' mode
elseif numel(hp) % && isempty(tr.highligth)
   % find leaves (children) below this branch
   children = false(1,tr.numLabels);
   children(hp(1)) = true;
   for ind = hp(1)-tr.numLeaves:-1:1
       if children(ind+tr.numLeaves)
           children(tr.tree(ind,:))=true;
       end
   end 
   
   % find and draw path to selected
   if strcmp(tr.editMode,'Inspect')
       [pathA,pathB] = path2sel(tr,hp(1)); 
       dis2sel = tr.x(find(pathA,1))+tr.x(find(pathB,1))...
                 -2*tr.x(find(pathA,1,'last'));
       if any(pathB)
          xx = [tr.x(pathA);NaN;tr.x(pathB)];
          yy = [tr.y(pathA);NaN;tr.y(pathB)];
          hh=zeros(2*numel(xx),1); hh(1:2:end)=1; hh=cumsum(hh);
          set(tr.hpathline,'XData',xx(hh(2:end)),...
             'YData',yy(hh(1:end-1)),'Visible','on');
       end
   end
   
   % place text to datatips
   name = [tr.names{hp(1)} ' '];
   %name(name=='_')=' ';
   children(hp(1)) = false;
   numChil = sum(children(1:tr.numLeaves));
   childrenNames = char(tr.names(children(1:tr.numLeaves)));
   %childrenNames(childrenNames=='_')=' ';
   childrenNames=[repmat('   ',size(childrenNames,1),1) childrenNames];
   switch tr.editMode
       case 'Inspect'
             set(tr.datatip,'string',char({...
                 [' Path length: '  num2str(dis2sel)];...
                 [' Selected node: ' tr.names{find(tr.selected,1)}];...
                 [' Current node:   ' name]}))
             extraLines = 3;
             numChil = 0;
       case {'Collapse/Expand','Rotate Branch','Rename','Prune','Select'}
           if numChil
               set(tr.datatip,'string',char([ 
                   {[name '  (' num2str(numChil) ' samples)']};...
                   mat2cell(childrenNames,ones(size(childrenNames,1),1),...
                   size(childrenNames,2))]))
           else
               set(tr.datatip,'string', name)
           end
           extraLines = 1;
       otherwise % all other modes
   end

   %compute some values before adjusting data tip
   fp = get(gcbf,'Position'); % fig position in points
   fh = fp(4);%fw = fp(3);     % fig size (height & width) in points
   ap  = get(tr.ha,'Position');          % axis position normalized
   yl  = ylim(tr.ha); yl  = yl - ...
         [ap(2) ap(2)+ap(4)-1]*diff(yl)/ap(4); % fig height limits in axis units
   xl  = xlim(tr.ha); xl  = xl - ...
         [ap(1) ap(1)+ap(3)-1]*diff(xl)/ap(3); % fig width limits in axis units
   yPosPt = (-yPos -4*yThres + yl(2))*fh/diff(yl); % datatip position in pts
   reqPt  = (numChil+extraLines)*14+2;  % required datatip height in pts
                                        % adjust if other fontsize is used 
                                         
   %adjust string of datatip if it will not fit (i.e. remove names)
   if reqPt > fh
       str = get(tr.datatip,'String');
       set(tr.datatip,'String',str(1:extraLines,:));
       reqPt = extraLines*14+2;
   end
   
   %adjust vertical position of datatip just below cp
   topEdge = yl(2)-min(fh,max(yPosPt,reqPt))*diff(yl)/fh;
   switch tr.editMode
       case {'Collapse/Expand','Rotate Branch','Prune'} 
       % datatip usually to the left of cp to see shadowing of branches    
           datatipExtent = get(tr.datatip,'Extent');
           datatipWidth = datatipExtent(3);
           rightEdge = max(xPos-3*xThres,xl(1)+datatipWidth);
           % is the datatip over cp ?
           if rightEdge>xPos && topEdge<yPos
               % then try to put it above cp
                topEdge = yPos - 3 * yThres - reqPt*diff(yl)/fh; 
               % does datatip fit above cp ?
               if topEdge<yl(1)
                   % then adjust string by removing names of species
                   str = get(tr.datatip,'String');
                   set(tr.datatip,'String',str(1:extraLines,:));
                   reqPt = extraLines*14+2;
                   topEdge = yl(2)-min(fh,max(yPosPt,reqPt))*diff(yl)/fh;
               end
           end
           set(tr.datatip,'Position',[rightEdge-datatipWidth,topEdge,1])
           set(tr.datatip,'horizontalalignment','left')
       case {'Inspect','Rename','Select'}
       % datatip usually to the right of cp to minimize problems on the
       % left edge
           set(tr.datatip,'Position',[xPos+3*xThres,topEdge,1])
           set(tr.datatip,'horizontalalignment','left')
       otherwise
   end

   switch tr.editMode
       case {'Collapse/Expand','Rotate Branch','Prune'}
           % de-color branches to rotate or collapse
           uncoloredNodes = false(1,tr.numLabels);
           uncoloredNodes(hp(1)) = true;
           for ind =  hp(1)-tr.numLeaves:-1:1
               if uncoloredNodes(ind+tr.numLeaves)
                   uncoloredNodes(tr.tree(ind,:))=true;
               end
           end
           if ~strcmp(tr.editMode,'Prune')
               uncoloredNodes(hp(1)) = false;
           end
           uncoloredNodes = uncoloredNodes & tr.activeNodes';
           set(tr.hlines(uncoloredNodes),'color',[.87 .87 .87])
           ind=find(uncoloredNodes(tr.numLeaves+1:tr.numLabels))+tr.numLeaves;
           set(tr.hldots(1),'Xdata',tr.x(ind),'Ydata',tr.y(ind))
           ind=find(uncoloredNodes(1:tr.numLeaves));
           set(tr.hldots(2),'Xdata',tr.x(ind),'Ydata',tr.y(ind))
       otherwise
   end
  
   if ismac
      set(tr.datatip,'FontSize',10) 
   end
   set(tr.datatip,'Visible','on')
   tr.indicativeMode = true;
   set(gcbf,'userdata',tr)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseRelease(~,varargin) 
tr = get(gcbf,'userdata');
if ~(strcmp(tr.editMode,'Select') && any(strcmp(get(gcbf,'SelectionType'),{'normal','extend'})))
    return
end

xdata = get(tr.hdragbox,'Xdata');
ydata = get(tr.hdragbox,'Ydata');
hp = (tr.x<max(xdata) & tr.x>min(xdata) & ...
          tr.y<max(ydata) & tr.y>min(ydata) & tr.activeNodes) ;
if (strcmp(get(gcbf,'SelectionType'),'normal') && ...
    strcmp(get(tr.hdragbox,'visible'),'on'))
    tr.selected(:) = false;
end
tr.selected(hp) = true;
tr.sel2root = path2root(tr,tr.selected);
set(tr.hdragbox,'Visible','off')
set(gcbf,'userdata',tr)
updateTree(gcbf,[])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectNode(~,varargin) 
% Callback function to select a Node. 
% Entry points: from 1) the dots context menu or 2) toggle node
tr = get(gcbf,'userdata');
% set a virtual grid to get the point
xThres=diff(get(tr.ha,'Xlim'))/100;
yThres=diff(get(tr.ha,'Ylim'))/100;
cp = get(tr.ha,'CurrentPoint');
xPos = cp(1,1); yPos = cp(1,2);
hp = find(tr.x<(xPos+xThres) & tr.x>(xPos-xThres) & ...
          tr.y<(yPos+yThres) & tr.y>(yPos-yThres)) ;

if numel(hp)
     set(tr.hdragbox,'visible','off')
     temp = tr.selected(hp(1));
     switch get(gcbf,'SelectionType')
         case 'normal'; tr.selected(:) = false;
         case 'alt'; if ~temp 
                        tr.selected(:) = false; 
                     end
                     temp=false;
     end
     tr.selected(hp(1)) = ~temp;
     tr.sel2root = path2root(tr,tr.selected);
     set(gcbf,'userdata',tr)
     updateTree(gcbf,[])
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inspectNode(~,varargin) 
% Callback function to inspect the reference Node. 
% Entry points: from 1) the dots context menu or 2) toggle node
tr = get(gcbf,'userdata');
% set a virtual grid to get the point
xThres=diff(get(tr.ha,'Xlim'))/100;
yThres=diff(get(tr.ha,'Ylim'))/100;
cp = get(tr.ha,'CurrentPoint');
xPos = cp(1,1); yPos = cp(1,2);
hp = find(tr.x<(xPos+xThres) & tr.x>(xPos-xThres) & ...
          tr.y<(yPos+yThres) & tr.y>(yPos-yThres)) ;
      
if numel(hp)
     temp = tr.selected(hp(1));
     tr.selected(:) = false;
     tr.selected(hp(1)) = ~temp;
     tr.sel2root = path2root(tr,tr.selected);
     set(gcbf,'userdata',tr)
     updateTree(gcbf,[])
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function collapseExpand(h,varargin) 
% Callback function to Collapse/Expand a branch. 
% Entry points: from 1) the dots context menu or 2) toggle node

tr = get(gcbf,'userdata');
if strcmp(get(h,'Type'),'uimenu') % come from a context menu
    hp = find(tr.selected(tr.numLeaves+1:tr.numLabels));
else
    % set a virtual grid to get the point
    xThres=diff(get(tr.ha,'Xlim'))/100;
    yThres=diff(get(tr.ha,'Ylim'))/100;
    cp = get(tr.ha,'CurrentPoint');
    xPos = cp(1,1); yPos = cp(1,2);
    hp = find(tr.x<(xPos+xThres) & tr.x>(xPos-xThres) & ...
              tr.y<(yPos+yThres) & tr.y>(yPos-yThres)) ;
    hp=hp(hp>tr.numLeaves)-tr.numLeaves;
    if numel(hp)
        hp=hp(1); %just in case it picked two points
    end
end

if numel(hp)
    for ind = 1:numel(hp)
        tr.activeBranches(hp(ind))=~tr.activeBranches(hp(ind));
        activeBranches=find(tr.activeBranches)';
        % find active nodes by expanding active Branches
        tr.activeNodes(:)=false;
        tr.activeNodes(tr.numLabels,1)=true;
        for k = activeBranches(end:-1:1)
            tr.activeNodes(tr.tree(k,:))=tr.activeNodes(k+tr.numLeaves);
        end
        tr.selected(:) = false;
        tr.selected(hp(ind)+tr.numLeaves) = true;
    end
    tr.sel2root = path2root(tr,tr.selected);
    set(gcbf,'userdata',tr)
    updateTree(gcbf,hp(end)+tr.numLeaves)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotateBranch(h,varargin) 
% Callback function to rotate a branch reordering the leaves. 
% Entry points: from 1) the dots context menu or 2) toggle node
tr = get(gcbf,'userdata');
if strcmp(get(h,'Type'),'uimenu') % come from a context menu
    hp = find(tr.selected(tr.numLeaves+1:tr.numLabels));
else
    % set a virtual grid to get the point
    xThres=diff(get(tr.ha,'Xlim'))/100;
    yThres=diff(get(tr.ha,'Ylim'))/100;
    cp = get(tr.ha,'CurrentPoint');
    xPos = cp(1,1); yPos = cp(1,2);
    hp = find(tr.x<(xPos+xThres) & tr.x>(xPos-xThres) & ...
        tr.y<(yPos+yThres) & tr.y>(yPos-yThres)) ;
    hp=hp(hp>tr.numLeaves)-tr.numLeaves;
    if numel(hp)
        hp=hp(1); %just in case it picked two points
    end
end
if numel(hp)
    for ind = 1:numel(hp)
        %find Leaves for every child
        childrenA = false(1,tr.numLabels);
        childrenA(tr.tree(hp(ind),1)) = true;
        for k = tr.tree(hp(ind),1)-tr.numLeaves:-1:1
            if childrenA(k+tr.numLeaves)
                childrenA(tr.tree(k,:))=true;
            end
        end
        childrenB = false(1,tr.numLabels);
        childrenB(tr.tree(hp(ind),2)) = true;
        for k = tr.tree(hp(ind),2)-tr.numLeaves:-1:1
            if childrenB(k+tr.numLeaves)
                childrenB(tr.tree(k,:))=true;
            end
        end
        permuta = 1:tr.numLabels;
        chA = find(childrenA(1:tr.numLeaves));
        chB = find(childrenB(1:tr.numLeaves));
        if chA(1)<chB(1)
            permuta([chA chB])=[chB chA];
        else
            permuta([chB chA])=[chA chB];
        end
        ipermuta = zeros(1,tr.numLabels);
        ipermuta(permuta)=1:tr.numLabels;
        tr.names = tr.names(permuta);
        tr.dist = tr.dist(permuta);
        tr.tree = ipermuta(tr.tree);
        tr.par = tr.par(permuta(1:end-1));
        tr.selected = tr.selected(permuta);
        tr.activeNodes = tr.activeNodes(permuta);
        tr.sel2root = tr.sel2root(permuta);
    end
    set(gcbf,'userdata',tr)
    updateTree(gcbf,[])
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function renameNode(h,varargin) 
% Renames a Node. Puts an uicontrol to input the new name.
% Entry points: from 1) the dots context menu or 2) toggle node
tr = get(gcbf,'userdata');
Xlim=get(tr.ha,'Xlim');Ylim=get(tr.ha,'Ylim');
aPos=get(tr.ha,'Position');
% set a virtual grid to get the point
xThres=diff(Xlim)/100;
yThres=diff(Ylim)/100;
cp = get(tr.ha,'CurrentPoint');
xPos = cp(1,1); yPos = cp(1,2);
hp = find(tr.x<(xPos+xThres) & tr.x>(xPos-xThres) & ...
          tr.y<(yPos+yThres) & tr.y>(yPos-yThres)) ;
if numel(hp)
    tr.previousMode = tr.editMode; tr.editMode = 'Renaming';
    xBoxPos = (.02+(xPos-Xlim(1))/diff(Xlim))*aPos(3)+aPos(1);
    yBoxPos = (.02+(Ylim(2)-yPos)/diff(Ylim))*aPos(4)+aPos(2);
    position=get(gcbf,'position');
    position=[position(3)*xBoxPos position(4)*yBoxPos 150 20];
    set(tr.editBox,'position',position);
    set(tr.editBox,'Visible','on','string',tr.names{hp(1)},'Value',hp(1))
    disableAllUI(h)
    disableMyContextMenus(h)
    set(gcbf,'WindowButtonMotionFcn',[]); % disable windows mouse motion
    set(gcbf,'userdata',tr)
 end
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doneRenaming(h,varargin)
% Output helper function to abandon the "Renaming" mode
tr = get(gcbf,'userdata');
tr.editMode = tr.previousMode;
g = get(tr.editBox,'Value');
proposedName = get(tr.editBox,'String');
if any(strcmp(tr.names([1:g-1 g+1:end]),proposedName))
    proposedName = tr.names{g};
    warndlg('Names for leaves and branches must be unique.','Warning','modal')
end
tr.names{g} = proposedName;
set(tr.anytext,'string',tr.names{g})    
tr.labelExtents(g) = get(tr.anytext,'Extent')*[0;0;1;0];
set(tr.editBox,'Visible','off')
set(gcbf,'userdata',tr)
updateTree(gcbf,[]);
set(tr.ha,'Ylim',get(tr.ha,'Ylim')) % fire sxis listener
enableAllUI(h)
enableMyContextMenus(h)
enableMyWindowButtonActions(h)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pruneTree(h,varargin) 
% Callback function to prune the tree, this function is complex because not
% only the basic structure is updated but all the other handles are
% also updated to contain the new tree.

% if changing to Prune, verify the warning has been accepted
propsForFigure = getappdata(gcbf,'propsForFigure');
if isequal(propsForFigure.PruneWarning,'NotDone')
     warndlg(['Pruning nodes cannot be undone. Before continuing,',...
            ' you may want to export the current tree to a new viewer.'],...
            'Warning','modal')
     setacceptedwarningtoothertools
     return % do not do this pruning
end

tr = get(gcbf,'userdata');
if strcmp(get(h,'Type'),'uimenu') % comes from a context menu
    hp = find(tr.selected);
else
    tr = get(gcbf,'userdata');
    Xlim=get(tr.ha,'Xlim');Ylim=get(tr.ha,'Ylim');
    % aPos=get(tr.ha,'Position');
    % set a virtual grid to get the point
    xThres=diff(Xlim)/100;
    yThres=diff(Ylim)/100;
    cp = get(tr.ha,'CurrentPoint');
    xPos = cp(1,1); yPos = cp(1,2);
    hp = find(tr.x<(xPos+xThres) & tr.x>(xPos-xThres) & ...
        tr.y<(yPos+yThres) & tr.y>(yPos-yThres));
    hp=hp(1); %just in case it picked two points
end

hp(hp==tr.numLabels)=[]; %cannot delete the root

while numel(hp)
    %find all nodes to purge (i.e. all descendants)
    children = false(1,tr.numLabels);
    children(hp(1)) = true;
    for k = hp(1)-tr.numLeaves:-1:1
        if children(k+tr.numLeaves)
            children(tr.tree(k,:))=true;
        end
    end
    mypar = tr.par(hp(1));                                     % parent
    if mypar < tr.numLabels  % my parent is NOT the root
        % connect brother to granparent
        mygrpar = tr.par(mypar);                                 % grandparent
        myuncle = setxor(tr.tree(mygrpar-tr.numLeaves,:),mypar); % uncle
        mybro = setxor(tr.tree(mypar-tr.numLeaves,:),hp(1));   % brother
        tr.tree(mygrpar-tr.numLeaves,:) = [myuncle mybro];
        tr.dist(mybro) = tr.dist(mybro) + tr.dist(mypar);
        temp = get(tr.hlines(mygrpar-tr.numLeaves),'Xdata');
        temp([1 4])=tr.x(tr.tree(mygrpar-tr.numLeaves,:));
        set(tr.hlines(mygrpar-tr.numLeaves),'Xdata',temp);
        highlight = [mybro,mygrpar];
    else % if my parent is the root, now I am the new root
        temp=cell2mat(get(tr.hlines,'Xdata'))-tr.dist(end);
        for k = 1:tr.numBranches
            set(tr.hlines(k),'Xdata',temp(k,:));
        end
        highlight = setxor(tr.tree(mypar-tr.numLeaves,:),hp(1));
    end
    children(mypar) = true; %also delete my par
    % find indexes to change tree
    permuta = 1:tr.numLabels;
    permuta(children) = [];
    ipermuta = zeros(1,tr.numLabels);
    ipermuta(permuta) = 1:length(permuta);
    permutaBranches = permuta(permuta>tr.numLeaves)-tr.numLeaves;
    % update all tree structure fields
    tr.names = tr.names(permuta);
    tr.dist = tr.dist(permuta);
    tr.tree = tr.tree(permutaBranches,:);
    tr.tree = ipermuta(tr.tree);
    if isempty(tr.tree) 
        return; 
    end % one leaf, no branches !
    tr = doBasicCalculations(tr);
    hlines = tr.hlines;
    tr.hlines = tr.hlines(permuta);
    delete(setxor(hlines,tr.hlines));
    tr.activeNodes = tr.activeNodes(permuta);
    tr.activeBranches = tr.activeBranches(permutaBranches);
    tr.selected = tr.selected(permuta);
    tr.selected(:) = false;
    tr.selected(ipermuta(highlight)) = true;
    tr.sel2root = path2root(tr,tr.selected);
    % update the vector with nodes to prune (node index has changed)
    hp=ipermuta(hp);
    hp(1)=[];
    hp(hp==0)=[];
    set(gcbf,'userdata',tr)
    updateTree(gcbf,[])
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectCommonAncestor(~,varargin) 
tr = get(gcbf,'userdata'); 
hl = tr.selected;
if ~any(hl)
    hl(end) = true;
end
for i = tr.numBranches:-1:1
    hl(tr.tree(i,:)) = hl(tr.tree(i,:)) | hl(tr.numLeaves+i);
end
hl(tr.numLeaves+1:end) = false;
hl = double(hl);
nl = sum(hl);
for i = 1:tr.numBranches
    hl(tr.numLeaves+i) = sum(hl(tr.tree(i,:)));
end
tr.selected = hl == nl;
tr.selected(find(tr.selected,1)+1:end) = false;
tr.sel2root = path2root(tr,tr.selected);
set(gcbf,'userdata',tr)
updateTree(gcbf,[])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectLeaves(~,varargin) 
tr = get(gcbf,'userdata'); 
hl = tr.selected;
if ~any(hl)
    hl(:) = true;
end
hl(tr.numLeaves+1:end) = false;
tr.selected = hl;
tr.sel2root = path2root(tr,tr.selected);
set(gcbf,'userdata',tr)
updateTree(gcbf,[])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function propagateSelection(~,varargin) 
tr = get(gcbf,'userdata'); 
hl = tr.selected;
% if ~any(hl)
%     hl(end) = true;
% end
for i = tr.numBranches:-1:1
    hl(tr.tree(i,:)) = hl(tr.tree(i,:)) | hl(tr.numLeaves+i);
end
tr.selected = hl;
tr.sel2root = path2root(tr,tr.selected);
set(gcbf,'userdata',tr)
updateTree(gcbf,[])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function swapSelection(~,varargin) 
tr = get(gcbf,'userdata'); 
tr.selected = ~ tr.selected;
tr.sel2root = path2root(tr,tr.selected);
set(gcbf,'userdata',tr)
updateTree(gcbf,[])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function findNode(~,varargin) 
treefig = gcbf;
tr = get(treefig,'userdata');
s = inputdlg('Regular Expression to match ?','Find Leaf/Branch',1);
if ~isempty(s)
    hc=regexpi(regexprep(tr.names,'_',' '),s);
    h = false(1,tr.numLabels);
    for ind = 1:tr.numLabels
        if ~isempty(hc{ind})
            h(ind)=true;
        end
    end
    hf = find(h);
    for ind = 1:length(hf)
        while ~tr.activeNodes(hf(ind))
            hf(ind)=tr.par(hf(ind));
        end
    end
    tr.selected(:) = false;
    tr.selected(hf) = true;
    tr.sel2root = path2root(tr,tr.selected); % update path to root
    set(treefig,'Userdata',tr);
    updateTree(treefig,[]) 

    % if selected are out of current view then fit the tree
    if (any(min(ylim(tr.ha))>tr.y(tr.selected)) || ...
        any(max(ylim(tr.ha))<tr.y(tr.selected)))
        autoFit
    end
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectTool(h,varargin) 
% sets the slider uicontrol and wait for the selection to be entered.
expandAll
changeEditMode(h)           % turns off any mode
tr = get(gcbf,'userdata'); 
if ~any(tr.selected) || all(tr.selected) % if none or all are selected use the root as anchor
    tr.selected(:) = false;
    tr.selected(end) = true;
end
tr.sel2root = path2root(tr,tr.selected);

D = pdist(phytree(tr),'squareform',1,'nodes','all');
tr.minDist2Sel = min(D(:,tr.selected),[],2);

set(gcbf,'userdata',tr)
updateTree(gcbf,[])

wS = get(gcbf,'Position');  % window dimensions
aP = get(tr.ha,'Position'); % axes position
set(tr.slider2,  'UserData',tr.selected); % save it in case of cancelling
set(tr.slider2,  'Position',[aP(1)*wS(3) wS(4)-30 max(1,aP(3)*wS(3)-80) 20])
set(tr.slider2tx,'Position',[sum(aP([1 3]))*wS(3)-80 wS(4)-30 60 20])
set(tr.slider2ok,'Position',[sum(aP([1 3]))*wS(3)-20 wS(4)-30 30 20])  
set([tr.slider2 tr.slider2tx tr.slider2ok],'visible','on')
set(tr.slider2,'max',max(tr.minDist2Sel),'value',0)
set(tr.slider2,'TooltipString',get(h,'Label'))
slider2Callback(h,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider2Callback(~,varargin) 
% this helper function serves the selectTool slider, it selects/unselects
% nodes based on the distance to the initially selected nodes 
tr = get(gcbf,'userdata');
Value = get(tr.slider2,'Value');
tr.selected = tr.minDist2Sel <= Value;
tr.sel2root = path2root(tr,tr.selected);
set(tr.slider2tx,'String',num2str(Value))
set(gcbf,'Userdata',tr);
updateTree(gcbf,[])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider2txCallback(~,varargin) 
% this helper function serves the edit text next to the slider, it
% selects/unselects nodes based on the distance to the initially selected nodes 
tr = get(gcbf,'userdata');
Value = str2double(get(tr.slider2tx,'String'));
if isnan(Value)
    Value = get(tr.slider2,'Value');
    set(tr.slider2tx,'String',num2str(Value))
end
Value = min(get(tr.slider2,'Max'),Value);
Value = max(Value,0);
set(tr.slider2,'Value',Value)
set(tr.slider2tx,'String',num2str(Value))
tr.selected = tr.minDist2Sel <= Value;
tr.sel2root = path2root(tr,tr.selected);
set(gcbf,'Userdata',tr);
updateTree(gcbf,[])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doSelectTool(h,varargin)
% this helper function inactivates the slider after selectTool and updated
% the tree. Entry point: the only way to get into this function is by the
% 'OK' uicontrol next to the slider.
tr = get(gcbf,'userdata');
set([tr.slider2,tr.slider2tx,tr.slider2ok],'Visible','off')
changeEditMode(h);
updateTree(gcbf,[]);
autoFit
enableAllUI(h)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cancelSelectTool(h,varargin) 
% this helper function cancels the threshold cut mode and returns to the 
% 'select' mode. Entry point: the ways to get into this function is
% by the 'CANCEL' uicontrol next to the slider (does not exist yet) or by
% mouse click over the axes during the slider mode.
tr = get(gcbf,'userdata');
set([tr.slider2,tr.slider2tx,tr.slider2ok],'Visible','off')
tr.selected = get(tr.slider2,'UserData');
tr.sel2root = path2root(tr,tr.selected);
set(gcbf,'Userdata',tr);
changeEditMode(h);
updateTree(gcbf,[]);
enableAllUI(h)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function thresholdCut(h,varargin) %#ok mLint is off because the implementation 
% of thresholdCut has been deactivated in R2009b, however we left the code
% in in case a user might want to roll it back. 
% thresholdCut sets the slider uicontrol and wait for the threshold cut to
% be entered. 
expandAll
changeEditMode(h)           % turns off any mode
tr = get(gcbf,'userdata');  
wS = get(gcbf,'Position');  % window dimensions
aP = get(tr.ha,'Position'); % axes position
set(tr.slider,  'Position',[aP(1)*wS(3) wS(4)-30 max(1,aP(3)*wS(3)-80) 20])
set(tr.slidertx,'Position',[sum(aP([1 3]))*wS(3)-80 wS(4)-30 60 20])
set(tr.sliderok,'Position',[sum(aP([1 3]))*wS(3)-20 wS(4)-30 30 20])  
set([tr.slider tr.slidertx tr.sliderok],'visible','on')
if isempty(tr.lastThresholdValue)
    set(tr.slider,'max',max(tr.x),'value',max(tr.x)*.75)
else
    set(tr.slider,'max',max(tr.x),'value',tr.lastThresholdValue)
end
set(tr.slider,'TooltipString',get(h,'Label'))
sliderCallback(h,varargin)
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sliderCallback(~,varargin) 
% this helper function serves both sliders (from root and from leaves), it
% hides selected nodes based on the current threshold cut
tr = get(gcbf,'userdata');
Value = get(tr.slider,'Value');
switch get(tr.slider,'TooltipString')
    case 'Distance to Leaves'
        tr.tocollapse = tr.dist2Leaf<(max(tr.x)-Value);
        set(tr.slidertx,'String',num2str(max(tr.x)-Value))
    case 'Distance to Root'
        tr.tocollapse = tr.x >= Value;
        set(tr.slidertx,'String',num2str(Value))
end
set(gcbf,'Userdata',tr);
toshow = [tr.tocollapse(tr.par(1:tr.numLabels-1));0]&tr.tocollapse;
mask = (1:tr.numLabels)'>tr.numLeaves;
% update light lines
set(tr.hlines(~toshow),'color','k')
set(tr.hlines(toshow),'color',[.87 .87 .87])
% update light dots
set(tr.hldots(1),'Ydata',tr.y(toshow&mask&tr.activeNodes),...
                 'Xdata',tr.x(toshow&mask&tr.activeNodes))
set(tr.hldots(2),'Ydata',tr.y(toshow&~mask&tr.activeNodes),...
                 'Xdata',tr.x(toshow&~mask&tr.activeNodes))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doThresholdCut(~,varargin) 
% this helper function inactivates nodes based on the threshold cut selected
% with the slider. Entry point: the only way to get into this function is
% by the 'OK' uicontrol next to the slider.
tr = get(gcbf,'userdata');
tr.activeBranches = tr.activeBranches & ~tr.tocollapse(tr.numLeaves+1:tr.numLabels);
% find active nodes by expanding active Branches
activeBranches=find(tr.activeBranches)';
tr.activeNodes(:)=false;
tr.activeNodes(tr.numLabels,1)=true;
for ind = activeBranches(end:-1:1)
    tr.activeNodes(tr.tree(ind,:))=tr.activeNodes(ind+tr.numLeaves);
end
tr.lastThresholdValue = get(tr.slider,'Value');
set([tr.slider,tr.slidertx,tr.sliderok],'Visible','off')
tr.selected(:) = false;
tr.sel2root = path2root(tr,tr.selected);
set(gcbf,'userdata',tr)
changeEditMode(h);
updateTree(gcbf,[]);
autoFit
enableAllUI(h)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cancelThresholdCut(h,varargin) 
% this helper function cancels the threshold cut mode and returns to the 
% 'select' mode. Entry point: the ways to get into this function is
% by the 'CANCEL' uicontrol next to the slider (does not exist yet) or by
% mouse click over the axes diring the slider mode.
tr = get(gcbf,'userdata');
set([tr.slider,tr.slidertx,tr.sliderok],'Visible','off')
changeEditMode(h);
updateTree(gcbf,[]);
enableAllUI(h)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function collapseSelected(~,varargin)
% Callback function to collapse selected nodes. 
% Entry points: from 1) the dots context menu or 2) toggle node

tr = get(gcbf,'userdata');
keepSelected = tr.selected;
hl = tr.selected;
if ~any(hl)
    hl(end) = true;
end
for i = tr.numBranches:-1:1
    hl(tr.tree(i,:)) = hl(tr.tree(i,:)) | hl(tr.numLeaves+i);
end
tr.selected = hl;
hp = find(tr.selected(tr.numLeaves+1:tr.numLabels));
if numel(hp)
    for ind = 1:numel(hp)
        tr.activeBranches(hp(ind))=false;
        activeBranches=find(tr.activeBranches)';
        % find active nodes by expanding active Branches
        tr.activeNodes(:)=false;
        tr.activeNodes(tr.numLabels,1)=true;
        for k = activeBranches(end:-1:1)
            tr.activeNodes(tr.tree(k,:))=tr.activeNodes(k+tr.numLeaves);
        end
        tr.selected(:) = false;
        tr.selected(hp(ind)+tr.numLeaves) = true;
    end
    tr.sel2root = path2root(tr,tr.selected);
    set(gcbf,'userdata',tr)
    updateTree(gcbf,hp(end)+tr.numLeaves)
end
tr.selected = keepSelected & tr.activeNodes;
set(gcbf,'userdata',tr)
updateTree(gcbf,[])
autoFit

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function expandSelected(~,varargin)
% Callback function to expand selected nodes. 
% Entry points: from 1) the dots context menu or 2) toggle node

tr = get(gcbf,'userdata');
keepSelected = tr.selected;
hl = tr.selected;
if ~any(hl)
    hl(end) = true;
end
for i = tr.numBranches:-1:1
    hl(tr.tree(i,:)) = hl(tr.tree(i,:)) | hl(tr.numLeaves+i);
end
tr.selected = hl;
hp = find(tr.selected(tr.numLeaves+1:tr.numLabels));
if numel(hp)
    for ind = 1:numel(hp)
        tr.activeBranches(hp(ind))=true;
        activeBranches=find(tr.activeBranches)';
        % find active nodes by expanding active Branches
        tr.activeNodes(:)=false;
        tr.activeNodes(tr.numLabels,1)=true;
        for k = activeBranches(end:-1:1)
            tr.activeNodes(tr.tree(k,:))=tr.activeNodes(k+tr.numLeaves);
        end
        tr.selected(:) = false;
        tr.selected(hp(ind)+tr.numLeaves) = true;
    end
    tr.sel2root = path2root(tr,tr.selected);
    set(gcbf,'userdata',tr)
    updateTree(gcbf,hp(end)+tr.numLeaves)
end
tr.selected = keepSelected & tr.activeNodes;
set(gcbf,'userdata',tr)
updateTree(gcbf,[])
autoFit

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function expandAll(~,varargin) 
% Callback function to expand all hidden nodes
tr = get(gcbf,'userdata');
x=[0 inf];
[dump,anchor]=min(abs((mean(ylim)-tr.y))+x(1+tr.activeNodes)'); %#ok
tr.activeBranches(:) = true;
tr.activeNodes(:) = true;
set(gcbf,'Userdata',tr);
updateTree(gcbf,anchor)
autoFit

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback function to save tree
function saveNewick(~,varargin)
tr = get(gcbf,'userdata');
newtr.tree = tr.tree;
newtr.dist = tr.dist;
newtr.names = tr.names;
phytreewrite(phytree(newtr),'GUI',true);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback function to restore the original tree
function restoreTree(~,varargin) 
tr = getappdata(gcbf,'backupTree');
view(phytree(tr))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback function to load tree
function loadNewick(h,varargin)
if strcmp(get(h,'Type'),'uimenu') % if caller is the uimenu then needs 
    figtoclose = gcbf;
    tr = phytreeread;             % to pick a file 
    
else                              % if not, caller is the callback from get workspace var
    tr=[];
    pfig = get(h,'Parent');
    if strcmp(get(h,'string'),'Import') || ...
      (strcmp(get(h,'style'),'listbox') && strcmp(get(gcbf,'SelectionType'),'open'))     
        hp = get(pfig,'Userdata'); 
        figtoclose = hp(4);
        ops = get(hp(1),'string');
        if ~isempty(ops)
            tr = evalin('base',ops(get(hp(1),'value'),:));
        end
        close(pfig);
    elseif strcmp(get(h,'string'),'Cancel')
        close(pfig);
    end   
end
if ~isempty(tr)
    propsForFigure = getappdata(figtoclose,'propsForFigure');
    view(tr,[],propsForFigure);
    close(figtoclose)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback function to copy the tree to a figure
function preparePublishFigure(h,varargin,AllNodes,callerIsContextMenu)
tr = get(gcbf,'userdata');
% if called from the uimenu and nothing is selected pick the root
if ~callerIsContextMenu & ~tr.selected %#ok
    tr.selected(end) = true; 
end

selected = find(tr.selected);
commonpath = true(tr.numLabels,1);
for ind = 1:numel(selected)
    commonpath = commonpath & path2root(tr,selected(ind));
end
branchtoexp = find(commonpath,1);
tr.selected(:) = false;
tr.selected(branchtoexp) = true;

set(gcbf,'userdata',tr)
updateTree(gcbf,[]);

hp = branchtoexp;
%find all nodes to export (i.e. all descendants) and also propagate
%non-active branches down the tree
children = false(1,tr.numLabels);
children(hp) = true;
active = [tr.activeNodes;tr.activeBranches];
for ind = hp-tr.numLeaves:-1:1
    if children(ind+tr.numLeaves)
        children(tr.tree(ind,:))=true;
    end
    active(ind+tr.numLeaves) = all(active(tr.tree(ind,:)) & active(ind+tr.numLeaves));
end
if AllNodes
    active(:) = true;
end

if sum(active(children))>1 % enough leaves to export ?
    publishdlg(h,AllNodes,callerIsContextMenu,find(children),active)
else
    warndlg('Not enough leaves.','Warning','modal')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function publishdlg(~,AllNodes,callerIsContextMenu,permuta,active) 
% dialog window to select options for publishing
vfig = gcbf;
c = get(0,'ScreenSize')*[1 0;0 1;.5 0;0 .5];
fig = figure('WindowStyle','modal','Color',[0.831373 0.815686 0.784314],...
             'Position',[c-[175 125] 350 250],'Resize','off','NumberTitle','off',...
             'Name','Print Phylogenetic Tree to Figure','IntegerHandle','off',...
             'Tag','PhytreetoolPrintDialogBox');
h1=uibuttongroup;h2=uibuttongroup;
set(h1,'Position',[.07 .25 .40 .70],'Title','Rendering Type','backgroundcolor',[0.831373 0.815686 0.784314])
set(h2,'Position',[.53 .35 .40 .45],'Title','Display Labels','backgroundcolor',[0.831373 0.815686 0.784314])
uiFS = get(h1, 'FontSize');
if ~ispc
    uiFS = uiFS -1;
    set(h1, 'FontSize', uiFS);
    set(h2, 'FontSize', uiFS');
end

ui1=uicontrol(h1,'style','radiobutton',...
                 'Position',[5 130 90 20],...
                 'string','Square',...
                 'Tag','SquareRadioButton',...
                 'value',1,...
                 'backgroundcolor',[0.831373 0.815686 0.784314],...
                 'FontSize', uiFS);
ui2=uicontrol(h1,'style','radiobutton',...
                 'Position',[5 100 90 20],...
                 'string','Angular',...
                 'Tag','AngularRadioButton',...
                 'backgroundcolor',[0.831373 0.815686 0.784314],...
                 'FontSize', uiFS);
ui3=uicontrol(h1,'style','radiobutton',...
                 'Position',[5 70 90 20],...
                 'string','Radial',...
                 'Tag','RadialRadioButton',...
                 'backgroundcolor',[0.831373 0.815686 0.784314],...
                 'FontSize', uiFS);
ui7=uicontrol(h1,'style','radiobutton',...
                  'Position',[5 40 90 20],...
                  'string','Equal-Angle',...
                  'Tag','EqualAngleRadioButton',...
                  'backgroundcolor',[0.831373 0.815686 0.784314],...
                 'FontSize', uiFS);
ui8=uicontrol(h1,'style','radiobutton',...
                 'Position',[5 10 109 20],...
                 'string','Equal-Daylight',...
                 'Tag','EqualDaylightRadioButton',...
                 'backgroundcolor',[0.831373 0.815686 0.784314],...
                 'FontSize', uiFS);
ui4=uicontrol(h2,'style','checkbox',...
                 'Position',[5 70 109 20],...
                 'string','Branch Nodes',...
                 'Tag','BranchNodesRadioButton',...
                 'backgroundcolor',[0.831373 0.815686 0.784314],...
                 'FontSize', uiFS);
ui5=uicontrol(h2,'style','checkbox',...
                 'Position',[5 40 109 20],...
                 'string','Leaf Nodes',...
                 'Tag','LeafNodesRadioButton',...
                 'backgroundcolor',[0.831373 0.815686 0.784314],...
                 'FontSize', uiFS);
ui6=uicontrol(h2,'style','checkbox',...
                 'Position',[5 10 109 20],...
                 'string','Terminal Nodes',...
                 'Tag','TerminalNodesRadioButton',...
                 'value',1,'backgroundcolor',[0.831373 0.815686 0.784314],...
                 'FontSize', uiFS);
uicontrol(fig,'style','pushbutton',...
              'Position',[90 20 60 30],...
              'string','Print',...
              'Tag','PrintPushButton',...
              'FontSize', uiFS',...
              'Callback',{@doPublishFigure,AllNodes,callerIsContextMenu,permuta,active});
uicontrol(fig,'style','pushbutton',...
              'Position',[195 20 60 30],...
              'string','Cancel',...
              'Tag','CancelPushButton',...
              'FontSize', uiFS',...
              'Callback',{@doPublishFigure,AllNodes,callerIsContextMenu,permuta,active});
set(fig,'Userdata',[ui1 ui2 ui3 ui4 ui5 ui6 ui7 ui8 vfig]);
set(h1,'SelectionChangeFcn',{@toggleCheckBoxs,[ui3,ui7,ui8],ui6})

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copy the tree to a figure
function doPublishFigure(h,varargin,~,~,permuta,active) 

pfig = get(h,'Parent');
if strcmp(get(h,'string'),'Cancel') 
    close(pfig); 
    return; 
end

% get options from the publishdlg window and close it
hp = get(pfig,'Userdata'); 
vfig = hp(9);
va = get(hp(1:8),'Value');
va = [va{:}];
switch find(va([1:3 7 8]),1)
    case 1; args = {'type','square'};
    case 2; args = {'type','angular'};
    case 3; args = {'type','radial'};
    case 4; args = {'type','equalangle'};
    case 5; args = {'type','equaldaylight'};
end
args = [args,{'bra',va(4),'lea',va(5),'ter',va(6)}];
close(pfig);

tr = get(vfig,'userdata');
newtr = phytree;
ipermuta(permuta) = 1:length(permuta);
numLeaves = (ipermuta(end) + 1)/2;
newtr.tree = ipermuta(tr.tree(permuta(numLeaves+1:end)-tr.numLeaves,:));
newtr.dist = tr.dist(permuta);
newtr.names = tr.names(permuta);
plot(newtr,active(permuta),args{:})

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback function for autofit 
function autoFit(~,varargin) 
tr = get(gcbf,'userdata');
set(tr.ha,'Ylim',[min(tr.y(tr.activeNodes))-1,max(tr.y(tr.activeNodes))+1]);
set(tr.ha,'Xlim',[0  max(tr.x)] + max(tr.x) * [-.1 .05]); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback function to Reset view
function myResetView(~,varargin) 
expandAll(gcbf)
tr = get(gcbf,'userdata');
set(tr.ha,'Ylim',tr.yLim);
set(tr.ha,'Xlim',tr.xLim);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common export function
function exportSubtree(~,varargin,AllNodes,ToWS,callerIsContextMenu) 
tr = get(gcbf,'userdata');
% if called from the uimenu and nothing is selected pick the root
if ~callerIsContextMenu & ~tr.selected %#ok
    tr.selected(end) = true; 
end

selected = find(tr.selected);
commonpath = true(tr.numLabels,1);
for ind = 1:numel(selected)
    commonpath = commonpath & path2root(tr,selected(ind));
end
branchtoexp = find(commonpath,1);
tr.selected(:) = false;
tr.selected(branchtoexp) = true;

set(gcbf,'userdata',tr)
updateTree(gcbf,[]);
doExport(branchtoexp,AllNodes,ToWS)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common export function after having selected a point
function doExport(hp,AllNodes,ToWS) 

tr = get(gcbf,'userdata');

%find all nodes to export (i.e. all descendants)
children = false(1,tr.numLabels);
children(hp) = true;
if AllNodes
    for ind = hp-tr.numLeaves:-1:1
        if children(ind+tr.numLeaves)
            children(tr.tree(ind,:))=true;
        end
    end
    permuta = find(children);
else
    for ind = hp-tr.numLeaves:-1:1
        if children(ind+tr.numLeaves)
            children(tr.tree(ind,:))=tr.activeNodes(tr.tree(ind,:));
        end
    end
    braToLea = find(~tr.activeNodes(tr.tree(:,1)))+tr.numLeaves;
    expBran = find(children(tr.numLeaves+1:end)) + tr.numLeaves;
    permuta = [find(children(1:tr.numLeaves)) ...
        intersect(expBran,braToLea)];
    [dump,hs] = sort(tr.y(permuta)); %#ok
    permuta = [permuta(hs) setdiff(expBran,braToLea)];
end
if sum(children)>1 % enough leaves to export ?
    newtr=phytree;
    ipermuta(permuta) = 1:length(permuta);
    numLeaves = (ipermuta(end) + 1)/2;
    newtr.tree = ipermuta(tr.tree(permuta(numLeaves+1:end)-tr.numLeaves,:));
    newtr.dist = tr.dist(permuta);
    newtr.names = tr.names(permuta);
    if ToWS % export to workspace ?
        s = inputdlg('Workspace variable name ?','Export to Workspace',1);
        while ~(isempty(s) || isvarname(s{1}) || isempty(s{1}))
            s = inputdlg('Not a valid variable name, type a MATLAB variable name ?','Export to Workspace',1);
        end
        if ~(isempty(s) || isempty(s{1}))
            assignin('base',s{1},newtr)
        end
    else % no, then export to other viewer
        view(newtr);
    end
else
    warndlg('Not enough leaves.','Warning','modal')
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateTree(h,anchor)
% Redraws the tree depending on the active Branches
% rather than erase and redraw, we only change specific fields in hlines
% and hdots.

tr = get(h,'userdata');
activeBranches=find(tr.activeBranches)';
oldPos = tr.y(anchor); 
    
% propagate last leaf
lastleaf = 1:tr.numLabels;
for ind = tr.numBranches:-1:1
    if ~tr.activeNodes(tr.tree(ind,1))
        lastleaf(tr.tree(ind,:))=lastleaf(ind+tr.numLeaves);
    end
end

% find x coordinates of branches
tr.x = tr.dist; 
for ind = tr.numBranches:-1:1
    tr.x(tr.tree(ind,:)) = tr.x(tr.tree(ind,:)) + tr.x(ind+tr.numLeaves);
end

% find y coordinates of branches
dummy = lastleaf([true,diff(lastleaf(1:tr.numLeaves))~=0]);
tr.y=zeros(tr.numLabels,1);
tr.y(dummy)=1:length(dummy);
for ind = activeBranches
    tr.y(ind+tr.numLeaves) = mean(tr.y(tr.tree(ind,:)));
end

% update right labels
todis = tr.names(dummy);
set(tr.ha,'ytick',1:length(dummy),'yticklabel',todis)
tr.maxLabelExtent = max(tr.labelExtents(dummy));

% show only active branches
set(tr.hlines,'Visible','off')
set(tr.hlines(tr.activeNodes),'Visible','on')

% update coordinates in lines
for ind = 1:tr.numLabels-1
     set(tr.hlines(ind),'Ydata',tr.y([ind,ind,tr.par(ind)]))
     set(tr.hlines(ind),'Xdata',tr.x([ind,tr.par([ind ind])]))
end
set(tr.hlines(tr.numLabels),'Ydata',tr.y(tr.numLabels)*[1 1 1])
set(tr.hlines(tr.numLabels),'Xdata',tr.x(tr.numLabels)*[1 1 1])
          
% update dots
mask = false(tr.numLabels,1); mask(1:tr.numLeaves) = true;
set(tr.hdots(1),'Ydata',tr.y(tr.activeNodes&~mask),'Xdata',tr.x(tr.activeNodes&~mask))
set(tr.hdots(2),'Ydata',tr.y(tr.activeNodes&mask),'Xdata',tr.x(tr.activeNodes&mask))            

% update red dots
set(tr.hseldots(1),'Ydata',tr.y(tr.activeNodes&~mask&tr.selected),'Xdata',tr.x(tr.activeNodes&~mask&tr.selected))
set(tr.hseldots(2),'Ydata',tr.y(tr.activeNodes&mask&tr.selected),'Xdata',tr.x(tr.activeNodes&mask&tr.selected))            

% set the axis holders          
set(tr.axhold,'Ydata',[0.5,max(tr.y(tr.activeNodes))+0.5])
if numel(oldPos)
    set(tr.ha,'ylim',get(tr.ha,'ylim')+tr.y(anchor)-oldPos);
else
    set(tr.ha,'ylim',get(tr.ha,'ylim')) % just touch 'YLim' such that the listener is triggered
end

% turn on indicative modes
tr.indicativeMode = false;
set([tr.datatip tr.hpathline],'visible','off')
set(tr.hlines,'color','black')
set(tr.hldots(1),'Xdata',[],'Ydata',[])
set(tr.hldots(2),'Xdata',[],'Ydata',[])

% save figure data
set(h,'Userdata',tr)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function path2r = path2root(tr,from)
% helper function, finds path to root
path2r = false(tr.numLabels,1);
if (numel(from)~=1 && sum(from)~=1) 
    return; 
end
path2r(from) = true;
temp = find(path2r);
if numel(temp)
    while temp~=tr.numLabels;
        temp = tr.par(temp);
        path2r(temp) = true;
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pathA,pathB] = path2sel(tr,from)
% helper function, finds path to selected node
path2rt = path2root(tr,from);
commonPath = tr.sel2root & path2rt;
commonPath(find(commonPath,1)) = false;
pathB = tr.sel2root & ~commonPath;
pathA = path2rt & ~commonPath;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tr = doBasicCalculations(tr)
% helper function to compute and find some features of the tree
tr = struct(tr);
tr.numBranches = size(tr.tree,1);
tr.numLeaves = tr.numBranches + 1;
tr.numLabels = tr.numBranches + tr.numLeaves; 

% obtain parents for every node
tr.par(tr.tree(:)) = tr.numLeaves + [1:tr.numBranches 1:tr.numBranches];

% calculate the distance to the closest leaf for every node
% needed for fast threshold cut
tr.dist2Leaf = zeros(tr.numLabels,1);
for ind = 1:tr.numBranches
    tr.dist2Leaf(ind+tr.numLeaves) = ...
       min(tr.dist2Leaf(tr.tree(ind,:))+tr.dist(tr.tree(ind,:)));
end

% calculate drawing coordinates for the tree: x coordinated will never
% change, but y coordinates may change depending on the active branches and
% nodes. 
tr.x = tr.dist; tr.y=[1:tr.numLeaves zeros(1,tr.numBranches)]';
for ind = tr.numBranches:-1:1
    tr.x(tr.tree(ind,:)) = tr.x(tr.tree(ind,:)) + tr.x(ind+tr.numLeaves);
end
for ind =1:tr.numBranches
    tr.y(ind+tr.numLeaves) = mean(tr.y(tr.tree(ind,:)));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function correctFigureSize(fig,recommendedHeight)
% helper function to increase initial figure size depending on the screen &
% tree sizes
screenSize = diff(reshape(get(0,'ScreenSize'),2,2),[],2)-100;
            % 100 gives extra space for the figure header and win toolbar
position = get(fig,'Position');
if recommendedHeight > position(4)
    if recommendedHeight < sum(position([2 4]))
        position(2) = sum(position([2 4])) - recommendedHeight;
        position(4) = recommendedHeight;
    elseif recommendedHeight < screenSize(2)
        position(2) = 30; 
        position(4) = recommendedHeight;
    else 
        position(2) = 30; 
        position(4) = screenSize(2);
    end
    set(fig,'Position',position)
end    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hToggleToolbar = makePhyTreeViewerToolbar(fig)
% helper function to set the toolbar
%
% hToggleToolbar contains handles to easy change the state on/off when
%                changing modes

oldSH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

set(fig,'toolbar','figure')  % needs to update because uicontrols turn it off

% Fix toolbar options, we keep: ZoomIn,ZoomOut,Pan
hw = findall(fig,'type','uitoolbar');
hf = get(hw,'Children');
h1 = findall(hf,'Tag','Exploration.Pan');
h2 = findall(hf,'Tag','Exploration.ZoomOut');
h3 = findall(hf,'Tag','Exploration.ZoomIn');
delete(setxor(hf,[h1,h2,h3]))
set([h1 h2 h3],'Separator','off','clickedCallback',@changeEditMode);

% load icons
icons = load(fullfile(toolboxdir('bioinfo'),'bioinfo','@phytree','phytreeicons'),'icons');

h4 = uitoggletool('ToolTip','Inspect Tool Mode','separator','on',...
                  'Tag','Inspect',        'CData',icons.icons(:,:,1:3)); 
h5 = uitoggletool('ToolTip','Collapse/Expand Branch Mode',...
                  'Tag','Collapse/Expand','CData',icons.icons(:,:,4:6)); 
h6 = uitoggletool('ToolTip','Rotate Branch Mode',...
                  'Tag','Rotate Branch',  'CData',icons.icons(:,:,7:9)); 
h7 = uitoggletool('ToolTip','Rename Leaf/Branch Mode',...
                  'Tag','Rename',         'CData',icons.icons(:,:,10:12)); 
h8 = uitoggletool('ToolTip','Prune (delete) Leaf/Branch Mode',...
                  'Tag','Prune',          'CData',icons.icons(:,:,13:15)); 
set([h4 h5 h6 h7 h8],'clickedCallback',@changeEditMode,'state','off',...
                     'Serializable','off','HandleVisibility','off');   
hToggleToolbar = [h4 h5 h6 h7 h8 h1 h2 h3];
set(0,'ShowHiddenHandles',oldSH)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hToggleUIMenu = makePhyTreeViewerUIMenus(fig)
% helper function to set UI menus
% 
% hToggleUIMenu contains handles to easy check on/off modes in the UI menu

oldSH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

% delete figure menus not used
h1 = findall(fig,'Type','uimenu', 'Tag','figMenuEdit');
h2 = findall(fig,'Type','uimenu', 'Tag','figMenuView');
h3 = findall(fig,'Type','uimenu', 'Tag','figMenuInsert');
h4 = findall(fig,'Type','uimenu', 'Tag','figMenuDesktop');
delete([h1,h2,h3,h4])

% Repair "File" menu
hw = findall(fig,'Type','uimenu', 'Tag','figMenuFile');
hf = get(hw,'children');
h1 = findall(hw,'Tag','figMenuFileExportSetup');
h3 = findall(hw,'Tag','figMenuFilePrintPreview');
h4 = findall(hw,'Tag','printMenu');
delete(setxor(hf,[h1,h3,h4]))

uimenu(hw,'Tag','NewToolUIM',                  'Label','New Viewer...', 'Position',1,'Callback','phytreeviewer')
uimenu(hw,'Tag','OpenUIM',                     'Label','Open...',     'Position',2,'Callback',@loadNewick)
uimenu(hw,'Tag','ImportFromWorkspaceUIM',      'Label','Import from Workspace...','Position',3,'Callback',@importfromwsdlg)
uimenu(hw,'Tag','OpenOriginalNewToolUIM',      'Label','Open Original in New Viewer','Position',4,'Callback',@restoreTree)
uimenu(hw,'Tag','SaveAsUIM',                   'Label','Save As...',  'Position',5,'Callback',@saveNewick,'Separator','on')
item0 = uimenu(hw ,'Tag','PrintToFigureUIM',   'Label','Print to Figure','Position',6);
  uimenu(item0,'Tag','WithHiddenNodesPTFUIM',  'Label','With Hidden Nodes...','Callback',{@preparePublishFigure,1,0});
  uimenu(item0,'Tag','OnlyDisplayedPTFUIM',    'Label','Only Displayed...',   'Callback',{@preparePublishFigure,0,0});
item1 = uimenu(hw,'Tag','ExportToNewToolUIM',  'Label','Export to New Viewer','Position',7);
    uimenu(item1,'Tag','WithHiddenNodesENTUIM','Label','With Hidden Nodes...','Callback',{@exportSubtree,1,0,0});
    uimenu(item1,'Tag','OnlyDisplayedENTUIM',  'Label','Only Displayed...',   'Callback',{@exportSubtree,0,0,0});
item2 = uimenu(hw,'Tag','ExportToWorkspaceUIM', 'Label','Export to Workspace','Position',8);
    uimenu(item2,'Tag','WithHiddenNodesEWUIM', 'Label','With Hidden Nodes...','Callback',{@exportSubtree,1,1,0});
    uimenu(item2,'Tag','OnlyDisplayedEWUIM',   'Label','Only Displayed...',   'Callback',{@exportSubtree,0,1,0});
uimenu(hw,'Tag','ExitUIM',                     'Label','Exit','Separator','on','Position',12,'Callback','close(gcbf)')
set(h1,'Separator','on')
    
% Repair "Tools" menu
%hw = findall(fig,'Type','uimenu','Label','&Tools');
hw = findall(fig,'Type','uimenu','Tag','figMenuTools');
hf = get(hw,'children');
h1 = findall(hw,'Tag','figMenuZoomIn');    set(h1,'Callback',{@changeEditMode,gcbo});
h2 = findall(hw,'Tag','figMenuZoomOut');   set(h2,'Callback',{@changeEditMode,gcbo});
h3 = findall(hw,'Tag','figMenuPan');       set(h3,'Callback',{@changeEditMode,gcbo});
h4 = findall(hw,'Tag','figMenuResetView'); set(h4,'Callback',{@myResetView,gcbo});
h5 = findall(hw,'Tag','figMenuOptions');
set([h1,h4],'separator','off')
delete(setxor(hf,[h1,h2,h3,h4,h5]))
delete(findall(h5,'Tag','figMenuOptionsDatatip'))
delete(findall(h5,'Tag','figMenuOptionsDataBar'))
h6 = uimenu(hw,'Label','Inspect','Position',1,'Callback',@changeEditMode);
h7 = uimenu(hw,'Label','Collapse/Expand','Position',2,...
                                              'Callback',@changeEditMode);
h8 = uimenu(hw,'Label','Rotate Branch','Position',3,...
                                              'Callback',@changeEditMode);
h9 = uimenu(hw,'Label','Rename','Position',4, 'Callback',@changeEditMode);
h10 = uimenu(hw,'Label','Prune','Position',5, 'Callback',@changeEditMode);
item4 = uimenu(hw,'Label','Select','Position',9,'Separator','on');
    uimenu(item4,'Label','Select By Distance','Callback',@selectTool);
    uimenu(item4,'Label','Select Common Ancestor','Callback',@selectCommonAncestor);   
    uimenu(item4,'Label','Select Leaves','Callback',@selectLeaves);   
    uimenu(item4,'Label','Propagate Selection','Callback',@propagateSelection);   
    uimenu(item4,'Label','Swap Selection','Callback',@swapSelection);   
uimenu(hw,'Label','Find Leaf/Branch...','Position',10,'Callback',@findNode); 
uimenu(hw,'Label','Collapse Selected','Position',11,'Callback',@collapseSelected);
uimenu(hw,'Label','Expand Selected','Position',12,'Callback',@expandSelected);
uimenu(hw,'Label','Expand All','Position',13,'Callback',@expandAll);
uimenu(hw,'Label','Fit to Window','Position',14,'Separator','on',...
                                                 'Callback',@autoFit);
set(h1,'Separator','on')

% Repair "Help" menu
hw = findall(fig,'Type','uimenu','Tag','figMenuHelp');
delete(get(hw,'children'));
uimenu(hw,'Label','Bioinformatics Toolbox Help','Position',1,'Callback',...
       'helpview(fullfile(docroot,''toolbox'',''bioinfo'',''bioinfo.map''),''bioinfo_product_page'')')
uimenu(hw,'Label','Phylogenetic Tree Viewer Help','Position',2,'Callback',...
       ['helpview(fullfile(docroot,''toolbox'',''bioinfo'', ''bioinfo.map'')' ...
        ',''phytreetool_reference'')' ] )
uimenu(hw,'Label','Examples','Position',3,'Separator','on',...
       'Callback','demo(''toolbox'',''bioinfo'')')   
tlbx = ver('bioinfo');
mailstr = ['web(''mailto:bioinfo-feedback@mathworks.com?subject=',...
           'Feedback%20for%20Phytreeviewer%20in%20Bioinformatics',...
           '%20Toolbox%20',tlbx(1).Version,''')'];
uimenu(hw,'Label','Send Feedback','Position',4,'Separator','on',...
       'Callback',mailstr);
set(0,'ShowHiddenHandles',oldSH)
hToggleUIMenu = [h6 h7 h8 h9 h10 h1 h2 h3];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hToggleContextMenu,hAxisContextMenu,hDotsContextMenu] = ...
         makePhyTreeViewerContextMenus(~) 
% helper function to set context menus
% 
% hToggleContextMenu contains handles to easy check on/off modes in the
%                    Context Menu (discontinued)
% hAxisContextMenu   contains handle to the Axis Context Menu to used to
%                    reactivate it
% hDotsContextMenu   contains handle to the Dots Context Menu to used to
%                    reactivate it

% set up context menu for the axes (when mouse is not over a node)
hcm1 = uicontextmenu('Callback',@hideActiveIndicators);

% For R2009b we are disallowing changing 'modes' from the contextmenu of
% the axes to avoid user confusion, this aligns better to other tools and
% the user still is able to change modes either from the toolbar of the
% figuremenu:

% h1 = uimenu(hcm1,'Label','Inspect',        'Callback',@changeEditMode);
% h2 = uimenu(hcm1,'Label','Collapse/Expand','Callback',@changeEditMode);
% h3 = uimenu(hcm1,'Label','Rotate Branch',  'Callback',@changeEditMode);
% h4 = uimenu(hcm1,'Label','Rename',         'Callback',@changeEditMode);
% h5 = uimenu(hcm1,'Label','Prune',          'Callback',@changeEditMode);

% For R2009b we are discontinuing the threshold collapse tool:

% item4 = uimenu(hcm1 , 'Label', 'Threshold Collapse','Separator','on');
%    uimenu(item4,'Label','Distance to Leaves',    'Callback',@thresholdCut);
%    uimenu(item4,'Label','Distance to Root',      'Callback',@thresholdCut);

item5 = uimenu(hcm1,'Tag','SelectCM',           'Label','Select');
    uimenu(item5,'Tag','SelectByDistanceCM',    'Label','Select By Distance',    'Callback',@selectTool);
    uimenu(item5,'Tag','SelectCommonAncestorCM','Label','Select Common Ancestor','Callback',@selectCommonAncestor);   
    uimenu(item5,'Tag','SelectLeavesCM',        'Label','Select Leaves',         'Callback',@selectLeaves);   
    uimenu(item5,'Tag','PropagateSelectionCM',  'Label','Propagate Selection',   'Callback',@propagateSelection);   
    uimenu(item5,'Tag','SwapSelectionCM',       'Label','Swap Selection',        'Callback',@swapSelection); 
uimenu(hcm1,'Tag','FindLeafBranchCM',           'Label','Find Leaf/Branch...',   'Callback',@findNode); 
uimenu(hcm1,'Tag','CollapseSelectedCM',         'Label','Collapse Selected',     'Callback',@collapseSelected);
uimenu(hcm1,'Tag','ExpandSelectedCM',           'Label','Expand Selected',       'Callback',@expandSelected);
uimenu(hcm1,'Tag','ExpandAllCM',                'Label','Expand All',            'Callback',@expandAll);
uimenu(hcm1,'Tag','FitToWindowCM',              'Label','Fit to Window',         'Callback',@autoFit,'Separator','on');
uimenu(hcm1,'Tag','ResetToOriginalViewCM',      'Label','Reset to Original View','Callback',@myResetView);

% context menu for dots (when mouse over a dot and right mouse button pressed)
hcm2 = uicontextmenu('Callback',@enableOptionsContextMenu);
   uimenu(hcm2,'Tag','CollapseExpandCM',     'Label','Collapse/Expand',          'Callback',@collapseExpand);
   uimenu(hcm2,'Tag','RotateBranchCM',       'Label','Rotate Branch',            'Callback',@rotateBranch);
   uimenu(hcm2,'Tag','RenameCM',             'Label','Rename',                   'Callback',@renameNode);
   uimenu(hcm2,'Tag','PruneCM',              'Label','Prune',                    'Callback',@pruneTree);
item0 = uimenu(hcm2,'Tag','PrintToFigureCM', 'Label','Print to Figure','Separator','on');   
   uimenu(item0,'Tag','WithHiddenNodesPFCM', 'Label','With Hidden Nodes...',     'Callback',{@preparePublishFigure,1,1});
   uimenu(item0,'Tag','OnlyDisplayedPFCM',   'Label','Only Displayed...',        'Callback',{@preparePublishFigure,0,1});
item1 = uimenu(hcm2,'Tag','ExportToNewToolCM','Label','Export to New Viewer');
   uimenu(item1,'Tag','WithHiddenNodesENTCM','Label','With Hidden Nodes...',     'Callback',{@exportSubtree,1,0,1});
   uimenu(item1,'Tag','OnlyDisplayedENTCM',  'Label','Only Displayed...',        'Callback',{@exportSubtree,0,0,1});
item2 = uimenu(hcm2,'Tag','ExportToWorkspaceCM','Label','Export to Workspace');
   uimenu(item2,'Tag','WithHiddenNodesEWCM', 'Label','With Hidden Nodes...',     'Callback',{@exportSubtree,1,1,1});
   uimenu(item2,'Tag','OnlyDisplayedEWCM',   'Label','Only Displayed...',        'Callback',{@exportSubtree,0,1,1});
   
% save context menus in my data structure to later restore them if desired
%hToggleContextMenu = [h1 h2 h3 h4 h5];
hToggleContextMenu = [];
hAxisContextMenu = hcm1;
hDotsContextMenu = hcm2;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enableOptionsContextMenu(h,varargin)
hideActiveIndicators
selectNode(h,varargin)
tr = get(gcbf,'Userdata');
hc = get(tr.hDotsContextMenu,'Children');
set(hc,'Enable','on')
if sum(tr.selected)~=1 
    set(hc(5),'Enable','off'); 
end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setupYLabelsListeners(hFig,hAxes)
% helper function to setsup listeners for the ylables, so we can detect if
% we would need to change the fontsize

% listens when the Ylim of axes has changed
YLimListener = addlistener(gca,'YLim',...
               'PostSet',@(hSrc,event)ylabelsListener(hSrc,event,hFig,hAxes));
% listens when Position of Figure has changed
PositionListener = addlistener(gcf,'Position',...
                   'PostSet',@(hSrc,event)ylabelsListener(hSrc,event,hFig,hAxes));
% store the listeners
setappdata(gcf,'PhyTreeListeners',[YLimListener, PositionListener]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggleCheckBoxs(~,~,h378,h6)
if any(cell2mat(get(h378,'value')))
    set(h6,'enable','off')
else
    set(h6,'enable','on')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function importfromwsdlg(~,varargin) 
% dialog window to select variable from workspace
mvars = evalin('base','whos');
mvars = mvars(strmatch('phytree',{mvars(:).class}));
vfig = gcbf;
c = get(0,'ScreenSize')*[1 0;0 1;.5 0;0 .5];
fig = figure('WindowStyle','modal','Color',[0.831373 0.815686 0.784314],...
             'Position',[c-[80 100] 160 220],'Resize','off','NumberTitle','off',...
             'Name','Get Phytree Object','IntegerHandle','off' );
ui1=uicontrol(fig,'style','list','Position',[22 70 120 120],'string',char(mvars(:).name),'backgroundcolor','w','Callback',@loadNewick);         
ui2=uicontrol(fig,'style','pushbutton','Position',[15 20 60 30],'string','Import','Callback',@loadNewick);
ui3=uicontrol(fig,'style','pushbutton','Position',[90 20 60 30],'string','Cancel','Callback',@loadNewick);         
uicontrol(fig,'style','text','Position',[20 190 140 20],'string','Select phytree object:','Horizontal','left','backgroundcolor',[0.831373 0.815686 0.784314])
set(fig,'Userdata',[ui1 ui2 ui3 vfig]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = getphytreetoolnumber()
% Computes the index number for this particular tool

% first, finds the used numbers so far
allFigs = findall(0,'tag','PhyTreeTool');
usedNumbers = zeros(1,numel(allFigs)+1);
baseName = 'Phylogenetic Tree ';
baseLen = length(baseName);
for i = 1:numel(allFigs)
    str = get(allFigs(i),'Name');
    usedNumbers(i) = str2double(str(baseLen:end));
end

% This is how we find the next index.  The rule is that we find the lowest
% integer value (non-zero and positive) not yet prescribed to a phytree
% tool, This is the same way MATLAB figures behave.
n = num2str(min(setdiff(1:(max(usedNumbers)+1),usedNumbers)));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = getacceptedwarningfromothertools(~) 
% Finds out if the pruning warning has been accepted
% first, finds the used numbers so far
allFigs = findall(0,'tag','PhyTreeTool');
n = 'NotDone';
for i = 1:numel(allFigs)
    propsForFigure = getappdata(allFigs(i),'propsForFigure');
    if isequal(propsForFigure.PruneWarning,'Done')
        n = 'Done';
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setacceptedwarningtoothertools()
% Set the pruning warning to all the other open tools as 'Done'
allFigs = findall(0,'tag','PhyTreeTool');
for i = 1:numel(allFigs)
    propsForFigure = getappdata(allFigs(i),'propsForFigure');
    propsForFigure.PruneWarning = 'Done';
    setappdata(allFigs(i),'propsForFigure',propsForFigure);
end
    
