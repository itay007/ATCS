function phytreeviewer(tr)
%PHYTREEVIEWER visualize and edit phylogenetic trees.
%
%   PHYTREEVIEWER is an interactive tool that allows to edit phylogenetic
%   trees. The GUI allows to do branch pruning, reorder, renaming, distance
%   exploration and read/write NEWICK formatted files.
%   
%   PHYTREEVIEWER(TREE) loads a phylogenetic tree object into the GUI. 
%
%   PHYTREEVIEWER(FILENAME) loads a NEWICK file into the GUI. 
% 
%   Example:
%      
%       phytreeviewer('pf00002.tree')
%       
%   See also BIRDFLUDEMO, HIVDEMO, PHYTREE, PHYTREE/PLOT, PHYTREE/VIEW,
%   PHYTREEREAD, PHYTREEWRITE.

% Copyright 2003-2012 The MathWorks, Inc.
 
% This function is the gateway entry to the PHYTREE/VIEW method.

if nargin == 0
    mvars = evalin('base','whos');
    mvars = mvars(strncmpi('phytree',{mvars(:).class}, 7));
    c = get(0,'ScreenSize')*[1 0;0 1;.5 0;0 .5];
    fig = figure('WindowStyle','modal','Color',[0.831373 0.815686 0.784314],...
            'Position',[c-[115 165] 230 330],'Resize','off','NumberTitle','off',...
            'Name','Open A Phylogenetic Tree','IntegerHandle','off',...
            'Tag','PhytreetoolOpenPhylogeneticTreeDialogBox');
    h1 = uibuttongroup('BorderWidth',1,'Position',[.03 .14 .94 .83],'Title','Choose tree source:','backgroundcolor',[0.831373 0.815686 0.784314]); %#ok
    ui1 = uicontrol(fig,'style','radiobutton',...
                        'Position',[15 275 190 20],...
                        'string','Import from workspace',...
                        'tag','ImportFromWorkspaceRadioButton',...
                        'value',1,...
                        'callback',@topRadioButtonPressed,...
                        'backgroundcolor',[0.831373 0.815686 0.784314]);
    ui2 = uicontrol(fig,'style','radiobutton',...
                        'Position',[15 110 190 20],...
                        'string','Open phylogenetic tree file',...
                        'tag','OpenPhylogeneticTreeFileRadioButton',...
                        'callback',@bottomRadioButtonPressed,...
                        'backgroundcolor',[0.831373 0.815686 0.784314]);
    ui3 = uicontrol(fig,'style','text',...
                        'Position',[30 250 170 20],...
                        'string','Select phytree object:',...
                        'Horizontal','left','backgroundcolor',[0.831373 0.815686 0.784314]);
    ui4 = uicontrol(fig,'style','list',...
                        'Position',[35 140 120 110],...
                        'BackgroundColor','w',...
                        'string',char(mvars(:).name),...
                        'tag','listUIC',...
                        'Callback',{@doOK,{mvars.name}});
    ui5 = uicontrol(fig,'style','text',...
                        'Position',[30 83 80 20],...
                        'string','File name:',...
                        'Horizontal','left',...
                        'Enable','off',...
                        'backgroundcolor',[0.831373 0.815686 0.784314]);
    ui6 = uicontrol(fig,'style','edit',...
                        'Position',[30 60 125 20],...
                        'BackgroundColor','w',...
                        'horizontalalignment','left',...
                        'Tag','EditBox',...
                        'Enable','off',...
                        'Callback',@doOK);
    ui7 = uicontrol(fig,'style','pushbutton',...
                        'Position',[160 55 60 30],...
                        'string','Browse...',...
                        'tag','BrowsePushButton',...
                        'Enable','off',...
                        'Callback',@browseFile);
    ui8 = uicontrol(fig,'style','pushbutton',...
                        'Position',[40 10 60 25],...
                        'string','OK',...
                        'tag','OKPushButton',...
                        'Callback',{@doOK,{mvars.name}});
    ui9 = uicontrol(fig,'style','pushbutton',...
                        'Position',[130 10 60 25],...
                        'string','Cancel',...
                        'tag','CancelPushButton',...
                        'Callback','delete(gcbf)');        
    setappdata(fig,'h',[ui1 ui2 ui3 ui4 ui5 ui6 ui7 ui8 ui9])
    if isempty(mvars)
        set(ui1,'Value',0,'enable','off');
        set(ui2,'Value',1);
        set([ui3 ui4],'enable','off');
        set([ui5 ui6 ui7],'enable','on');
    end
    
    uiwait(fig) % wait until action is commanded
    return
end
    
if isa(tr,'phytree') 
    view(tr); 
    return; 
end

% at this point tr is a valid filename or some corrupted input, trying to
% read it with phytreeread

if (~ischar(tr) || size(tr,1)~=1)
    error(message('bioinfo:phytreeviewer:InvalidInput'))
end

tr = phytreeread(tr);
view(tr);


function topRadioButtonPressed(fh,~)
 uiresume(get(fh, 'parent'))
 hui = getappdata(gcbf,'h');
 set(hui(5:7),'enable','off');
 set(hui(3:4),'enable','on');
 set(hui(1),'value',1);
 set(hui(2),'value',0);

function bottomRadioButtonPressed(fh,~)
 uiresume(get(fh, 'parent'))
 hui = getappdata(gcbf,'h');
 set(hui(5:7),'enable','on');
 set(hui(3:4),'enable','off');
 set(hui(1),'value',0);
 set(hui(2),'value',1);
 
function browseFile(~,~)
 [filename, pathname] = uigetfile({'*.tree';'*.dnd'},'Select Phylogenetic Tree File');
 if ~filename
     filename=[];
 else
     filename = [pathname, filename];
 end
 hui = getappdata(gcbf,'h');
 set(hui(6),'string',filename)
 
function doOK(h,~,varargin) 
 caller = get(h,'style');
 hui = getappdata(gcbf,'h');
 if (strcmp(caller,'listbox') && strcmp(get(gcbf,'SelectionType'),'open'))...
    || (strcmp(caller,'pushbutton') && get(hui(1),'value')==1)
     mvars = varargin{1};
     view(evalin('base',mvars{get(hui(4),'value')}));
     delete(gcbf)
 elseif get(hui(2),'value')==1
     % it got here because we will try a file
     filename = get(hui(6),'string');
     if ~isempty(filename)
         if ~exist(filename,'file')
             errordlg('File not found.','Open Phylogenetic Tree','modal') 
         else
             try
                 tr = phytreeread(filename);
                 view(tr);
                 delete(gcbf);
             catch allExceptions %#ok<NASGU>
                 errordlg('Unable to determine the file format.',...
                          'Open Phylogenetic Tree','modal') 
             end
         end
     end
 end
    
