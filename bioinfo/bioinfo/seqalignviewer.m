function varargout = seqalignviewer(varargin)
%SEQALIGNVIEWER visualize and edit a sequence alignment.
%
%   SEQALIGNVIEWER is an interactive tool for displaying and editing a
%   multiple sequence alignment.
%
%   SEQALIGNVIEWER(ALIGNMENT) loads a group of previously aligned sequences
%   into the viewer. ALIGNMENT can be a structure with a field Sequence, a
%   character array, or a filename.
%
%   SEQALIGNVIEWER(...,'ALPHABET',ALPHA) specifies the aligned sequences
%   are amino acids ('AA') or nucleotides ('NT'). The default is AA. If
%   ALPHABET is not specified, SEQALIGNVIEWER will guess the alphabet
%   type.
%
%   SEQALIGNVIEWER(...,'SEQHEADERS',NAMES) passes a list of names to label
%   the sequences in the interactive tool. NAMES can be a vector of
%   structures with the fields 'Header' or 'Name', or a cell array of
%   strings. In both cases the number of elements provided must comply with
%   the number of sequences in ALIGNMENT.
%
%   SEQALIGNVIEWER(...,'R2012b',TRUE) uses a previous version of this
%   interactive tool. This option will be removed in the future.
%
%   Example:
%
%       seqalignviewer('aagag.aln')
%
%   See also BIRDFLUDEMO, FASTAREAD, GETHMMALIGNMENT, MULTIALIGN,
%   MULTIALIGNREAD, SEQVIEWER.

%    Copyright 2003-2012 The MathWorks, Inc.

% SEQALIGNVIEWER also provides services to the legacy java viewer:
%
%   seqalignviewer(...,'FROMVIEWER',true) a flag to indicate that the
%   viewer is calling this function.
%
%   seqalignviewer(...,'VARNAME',name) pass in the import variable name.

% Undocumented option to close the legacy java viewer
if nargin == 1 && ischar(varargin{1}) && strncmpi(varargin{1}, 'close', 5)
    seqviewer.MSAUtil.msaMatlab('manage_msaviewer', 'close');
    if nargout > 0
        varargout{1} = [];
        varargout{2} = [];
    end
    return;
end

% Distinguish among different signatures for the input arguments
if nargin == 0
    alignment = [];
elseif ischar(varargin{1}) && isrow(varargin{1}) && strncmpi('r2012b',varargin{1},numel(varargin{1}))
    % Detect if PVP 'R2012b' is used as the first input argument
    alignment = [];
else
    alignment = varargin{1};
    varargin(1) = [];
end

% Parse input PV pairs.
[headers,alphabet,varname,from_viewer,r2012b] =  parse_inputs(varargin{:});

if r2012b % Legacy java viewer needs to check if java is available
    if ~jvm_available
        error(message('bioinfo:seqalignviewer:JavaMissing'));
    end
end

% Shortcut for the case data is not given and we try to open an empty viewer
if isempty(alignment)
    if r2012b % if legacy java viewer and there is a viewer present do NOT open a blank page
        if ~(seqviewer.MSAUtil.msaMatlab('manage_msaviewer', 'state')) % viewer not present
            seqviewer.MSAUtil.msaMatlab('create_ui'); % then open a blank page
        end
        if nargout > 0
            varargout{1} = [];
            varargout{2} = [];
        end
    else % for the HG viewer
        % if there is a viewer present just bring the last viewer to focus
        allFigs = findall(0,'tag','seqAlignViewer.figure');
        if isempty(allFigs)
            sav; % call HG viewer with no data
        elseif numel(allFigs)==1
            figure(allFigs)
        else
            [~,wf] = sort(str2num(char(regexprep(get(allFigs,'Name'),'\D',''))),'descend');
            figure(allFigs(wf(1)))
        end
    end
    return;
end

[sequences,headers,alphabet,varname,ME] = validateAlignment(alignment,headers,alphabet,varname);

if r2012b % Legacy java viewer
    if isempty(ME)
        % open the alignment in viewer
        seqviewer.MSAUtil.msaMatlab('create_ui', varname, headers, char(sequences), alphabet);
        msg = [];
        msgid = [];
    else % there was an error
        msg = ME.message;
        msgid = ME.identifier;
        if ~from_viewer
            if(~seqviewer.MSAUtil.msaMatlab('manage_msaviewer', 'state')) % no viewer present
                rethrow(ME);
            else % show error in viewer
                seqviewer.MSAUtil.msaMatlab('create_ui', msg);
            end
        end
    end
    if nargout > 0
        varargout{1} = msg;
    end
    if nargout == 2
        varargout{2} = msgid;
    end
else % HG viewer
    if isempty(ME) % no error when validating the alignment
        sav(sequences,headers,alphabet); % call HG viewer
    else
        rethrow(ME);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = jvm_available

persistent sRlt;

if isempty(sRlt)
    sRlt=usejava('jvm') & usejava('awt') & usejava('swing');
end
result = sRlt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sequences,headers,alphabet,varname,ME] = validateAlignment(alignment,headers,alphabet,varname)

% All validating code is in a try-catch block because if seqalignviewer was
% called from the GUI itself, then the error message is just handed over as
% an output argument instead of throwing an error and terminating.
ME = [];
sequences =[];
try
    % get the variable name of the first input argument (needed in the
    % legacy java viewer for the window title)
    if isempty(varname)
        varname = inputname(1);
    end
    
    % == handle the case when alignment can be a filename
    if(ischar(alignment) && size(alignment,1)==1 && size(fullfile(pwd,alignment), 2) && (exist(alignment,'file') || exist(fullfile(pwd,alignment),'file')))
        fileName = alignment;
        [~, name, ext] = fileparts(fileName);
        varname = [name ext]; % (needed in the legacy java viewer for the window title)
        try
            alignment = multialignread(fileName);
        catch ME1
            try
                alignment = fastaread(fileName);
            catch ME2 %#ok<NASGU>
                try
                    error(message('bioinfo:seqalignviewer:IncorrectDataFormat'))
                catch ME3
                    rethrow(addCause(ME3,ME1)) %The excepttion in the second try-catch is ignored
                end
            end
        end
    end
    
    % == handle the case can be a char array or a numeric array
    if((ischar(alignment) || (isnumeric(alignment) && isreal(alignment) && isinteger(alignment)))&& size(alignment,1)>1) % char array alignment
        if size(alignment,1)==3 && all(alignment(2,:)==' ' | alignment(2,:)=='|' | alignment(2,:)==':')
            % it is a pairwise alignment with a matching row
            alignment = alignment([1 3],:);
        end
        numseq = size(alignment,1);
        s = [];
        for i=numseq:-1:1
            s(i).Header = sprintf('Sequence %d', i);
            s(i).Sequence = alignment(i,:);
        end
        alignment = s;
    end
    
    % == by now alignment must be a structure and must have the field 'Sequence'
    if ~isstruct(alignment) || ~isfield(alignment,'Sequence')
        error(message('bioinfo:seqalignviewer:UnknownInputType'));
    end
    
    % == try to get headers when they were not specified by the user
    if isempty(headers)
        if isfield(alignment,'Header')
            headers = {alignment(:).Header};
        elseif isfield(alignment,'Name')
            headers = {alignment(:).Name};
        elseif isfield(alignment,'LocusName')
            headers = {alignment.LocusName};
        else % create default headers
            headers = cell(numel(alignment),1);
            for i=1:numel(alignment)
                headers{i} = sprintf('Sequence %d', i);
            end
        end
    end
    
    sequences = {alignment(:).Sequence};
    
    % == handle numeric sequences
    if all(cellfun(@isnumeric,sequences))
        if strcmpi(alphabet,'nt')
            for i = 1:numel(sequences)
                sequences{i}  = int2nt(sequences{i});
            end
        else % for 'aa' or '', if alpha is not specified and input is
            % numeric we do not guess, we assume 'aa'.
            for i = 1:numel(sequences)
                sequences{i}  =  int2aa(sequences{i});
            end
        end
    elseif any(cellfun(@isnumeric,sequences))
        error(message('bioinfo:seqalignviewer:InputNotValid'));
    end
    
    % == polish the sequences - same as in multialign.m
    sequences = strtrim(sequences(:));     % trim sequences
    sequences = strrep(sequences,' ','-'); % alignment gap char can only be '-'
    sequences = strrep(sequences,'.','-'); % alignment gap char can only be '-'
    
    % == check if all the sequences are empty
    if sum(cellfun('length',sequences))==0
        error(message('bioinfo:seqalignviewer:InputNotValid'));
    end
    
    % == check if the input is aligned
    if numel(unique(cellfun('length',sequences)))>1
        error(message('bioinfo:seqalignviewer:SequencesNotAligned'));
    end
    
    % == guess the type of alphabet when alphabet was not specified or has
    % not been guessed (in case of numeric data)
    if isempty(alphabet)
        if bioinfoprivate.isnt(char(sequences))
            alphabet = 'nt';
        elseif bioinfoprivate.isaa(char(sequences))
            alphabet = 'aa';
        else
            error(message('bioinfo:seqalignviewer:InvalidSymbolsInSequences'));
        end
    elseif strcmp(alphabet,'nt') && ~bioinfoprivate.isnt(char(sequences)) ||...
            strcmp(alphabet,'aa') && ~bioinfoprivate.isaa(char(sequences))
        error(message('bioinfo:seqalignviewer:UnmatchingSymbolsInSequences'));
    end
    
    if numel(sequences) < 2
        error(message('bioinfo:seqalignviewer:TooFewSequences'));
    end
    
    if numel(sequences) ~= numel(headers)
        error(message('bioinfo:seqalignviewer:InccorrectSizeSEQHEADERS'));
    end
    
catch ME
    % just continue
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [headers,alphabet,varname,from_viewer,r2012b] =  parse_inputs(varargin)
% Parse input PV pairs.

% default values
headers = []; % by default no sequence names are given
alphabet   = ''; % can be empty because if not given we'll try to guess it from the input
from_viewer = false;
varname = '';
r2012b = false;

if rem(nargin, 2) ~= 0
    error(message('bioinfo:seqalignviewer:IncorrectNumberOfArguments', mfilename));
end

okargs = {'alphabet', 'fromviewer', 'varname', 'seqheaders','r2012b'};
for j=1:2:nargin-1
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1 % ALPHABET
            if  bioinfoprivate.optAlphabet(pval,okargs{k}, mfilename)
                alphabet = 'aa';
            else
                alphabet = 'nt';
            end
        case 2  % FROMVIEWER (undocummented, used in the legacy java viewer)
            from_viewer = bioinfoprivate.opttf(pval,okargs{k}, mfilename);
        case 3  % VARNAME (undocummented, used in the legacy java viewer)
            varname = pval;
        case 4  % SEQHEADERS
            if iscell(pval) || isfield(pval,'Header') || isfield(pval,'Name') || isfield(pval,'LocusName')
                if isfield(pval,'Header') % if struct put them in a cell
                    headers = {pval(:).Header};
                elseif isfield(pval,'Name') % if struct put them in a cell
                    headers = {pval(:).Name};
                elseif isfield(pval,'LocusName') % if struct put them in a cell
                    headers = {pval(:).LocusName};
                else
                    headers = pval;
                end
                headers = headers(:);
            else
                error(message('bioinfo:seqalignviewer:InvalidSEQHEADERS'))
            end
        case 5  % R2012b
            r2012b = bioinfoprivate.opttf(pval,okargs{k}, mfilename);
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HG based viewer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sav(seq,hea,alpha)
% seq is always a cell string with the aligned sequences
% hea is always a cell string with the sequence headers

if nargin
    data = repmat(struct('Header',[],'Sequence',[]),numel(hea),1);
    for j = 1:numel(hea)
        data(j).Sequence = seq{j};
        data(j).Header = hea{j};
    end
else
    data = [];
    alpha = 'aa';
end

% Hauptfigure öffnen
ssz = get(0,'ScreenSize');
h.hFig = figure( ...
    'Tag','seqAlignViewer.figure', ...
    'units','pixel', ...
    'position',[(ssz(3:4)-[1000 600])/2 1000 600], ...
    'HandleVisibility','callback', ...
    'menubar','none', ...
    'NumberTitle','off', ...
    'IntegerHandle','off',...
    'Name',bioinfoprivate.indexedFigureName('seqAlignViewer.figure','Biological Sequence Alignment -'), ...
    'visible','on');

% initialize empty uisav-object
obj = bioinfoprivate.uisav('Parent',h.hFig);

%load icons
iconFile = fullfile(toolboxdir('bioinfo'),'bioinfo','+bioinfoprivate','seqAlignViewerIcons.mat');
icon = load(iconFile);

% ToolBar initialisieren
h.UI.ht = uitoolbar(h.hFig,'Tag','seqAlignViewer.Toolbar');
uipushtool(h.UI.ht,'Tag','seqAlignViewer.Toolbar.ZoomOut',...
    'CData',icon.FontOut,...
    'TooltipString','Zoom out',...
    'clickedcallback',@ZoomOut);
uipushtool(h.UI.ht,'Tag','seqAlignViewer.Toolbar.ZoomDefault',...
    'CData',icon.FontDefault,...
    'TooltipString','Reset to Default FontSize',...
    'clickedcallback',@ZoomDefault);
uipushtool(h.UI.ht,'Tag','seqAlignViewer.Toolbar.ZoomIn',...
    'CData',icon.FontIn,...
    'TooltipString','Zoom in',...
    'clickedcallback',@ZoomIn);
h.UI.htb = uitoggletool(h.UI.ht,'Tag','seqAlignViewer.Toolbar.BackGroundToggle',...
    'CData',icon.FontBackground,...
    'Separator','on',...
    'TooltipString','Turn Background On/Off',...
    'State','on',...
    'clickedcallback',@BackGroundToggleToolbar);
uipushtool(h.UI.ht,'Tag','seqAlignViewer.Toolbar.MoveUp',...
    'CData',icon.SeqMoveUp,...
    'Separator','on',...
    'TooltipString','Move Row(s) up',...
    'clickedcallback',@(~,~) obj.MoveSelectedBlock(-1,2));
uipushtool(h.UI.ht,'Tag','seqAlignViewer.Toolbar.MoveDown',...
    'CData',icon.SeqMoveDown,...
    'TooltipString','Move Row(s) down',...
    'clickedcallback',@(~,~) obj.MoveSelectedBlock(1,2));
uipushtool(h.UI.ht,'Tag','seqAlignViewer.Toolbar.MoveTop',...
    'CData',icon.SeqMoveTop,...
    'TooltipString','Move Row(s) to Top',...
    'clickedcallback',@(~,~) obj.MoveSelectedBlock(-inf,2));
uipushtool(h.UI.ht,'Tag','seqAlignViewer.Toolbar.MoveBottom',...
    'CData',icon.SeqMoveBottom,...
    'TooltipString','Move Row(s) to Bottom',...
    'clickedcallback',@(~,~) obj.MoveSelectedBlock(inf,2));
% Behaivor TBD
% uipushtool(h.UI.ht,'Tag','seqAlignViewer.Toolbar.MoveLeft',...
%     'CData',icon.SeqMoveLeft,...
%     'Separator','on',...
%     'TooltipString','Move Selection Left',...
%     'clickedcallback',@(~,~) obj.MoveSelectedBlock(-1,1));
% uipushtool(h.UI.ht,'Tag','seqAlignViewer.Toolbar.MoveRight',...
%     'CData',icon.SeqMoveRight,...
%     'TooltipString','Move Selection Right',...
%     'clickedcallback',@(~,~) obj.MoveSelectedBlock(1,1));

% UIMENU initialisieren
h.UI.hf = uimenu('parent',h.hFig,'tag','seqAlignViewer.uimenu.File',...
    'Label','File',...
    'Callback',@UpdateFileMenu);
uimenu('parent',h.UI.hf,'tag','seqAlignViewer.uimenu.File.Open',...
    'Label','Open...',...
    'Accelerator','O',...
    'Callback',@Open );
uimenu('parent',h.UI.hf,'tag','seqAlignViewer.uimenu.File.ImportFromWorkspace',...
    'Label','Import from MATLAB Workspace...',...
    'Callback',@ImportFromWorkspace );
h.UI.hfe1 = uimenu('parent',h.UI.hf,'tag','seqAlignViewer.uimenu.File.Export2Fasta',...
    'Label','Export to FASTA File...', ...
    'Accelerator','S', ...
    'Separator','on',...
    'Callback',@(~,~) Export2Fasta('selection') );
h.UI.hfe2 = uimenu('parent',h.UI.hf,'tag','seqAlignViewer.uimenu.File.ExportConsensus2Fasta',...
    'Label','Export Consensus to FASTA File...',...
    'Callback',@(~,~) Export2Fasta('consensus'));
h.UI.hfe3 = uimenu('parent',h.UI.hf,'tag','seqAlignViewer.uimenu.File.Export2Workspace',...
    'Label','Export to MATLAB Workspace...',...
    'Callback',@(~,~) Export2Workspace('selection'));
h.UI.hfe4 = uimenu('parent',h.UI.hf,'tag','seqAlignViewer.uimenu.File.ExportConsensus2Workspace',...
    'Label','Export Consensus to MATLAB Workspace...',...
    'Callback',@(~,~) Export2Workspace('consensus'));
h.UI.hfp = uimenu('parent',h.UI.hf,'tag','seqAlignViewer.uimenu.File.Print',...
    'Label','Print', ...
    'Accelerator','P', ...
    'Separator','on',...
    'Callback','printdlg(gcbf)');
uimenu('parent',h.UI.hf,'tag','seqAlignViewer.uimenu.File.Close',...
    'Label','Close Viewer', ...
    'Accelerator','W', ...
    'Separator','on',...
    'Callback',@(~,~) delete(h.hFig));

h.UI.he = uimenu('parent',h.hFig,'tag','seqAlignViewer.uimenu.Edit',...
    'Label','Edit',...
    'Callback',@UpdateEditMenu);
h.UI.hec = uimenu('parent',h.UI.he,'tag','seqAlignViewer.uimenu.Edit.Copy',...
    'Label','Copy',...
    'Accelerator','C',...
    'Callback',@Export2Clipboard);
h.UI.hed = uimenu('parent',h.UI.he,'tag','seqAlignViewer.uimenu.Edit.DeleteSequences',...
    'Label','Delete Sequences',...
    'Callback',@(~,~) obj.DeleteSequences);
h.UI.hesa = uimenu('parent',h.UI.he,'tag','seqAlignViewer.uimenu.Edit.SelectAll',...
    'Label','Select All',...
    'Separator','on',...
    'Accelerator','A',...
    'Callback',@(~,~) obj.Select(1));
h.UI.heda = uimenu('parent',h.UI.he,'tag','seqAlignViewer.uimenu.Edit.DeselectAll',...
    'Label','Deselect All',...
    'Callback',@(~,~) obj.Select(0));
h.UI.hemu = uimenu('parent',h.UI.he,'tag','seqAlignViewer.uimenu.Edit.MoveRowsUp',...
    'Label','Move Rows(s) up',...
    'Separator','on',...
    'Callback',@(~,~) obj.MoveSelectedBlock(-1,2));
h.UI.hemd = uimenu('parent',h.UI.he,'tag','seqAlignViewer.uimenu.Edit.MoveRowsDown',...
    'Label','Move Rows(s) down',...
    'Callback',@(~,~) obj.MoveSelectedBlock(1,2));
h.UI.hemt = uimenu('parent',h.UI.he,'tag','seqAlignViewer.uimenu.Edit.MoveRowsTop',...
    'Label','Move Rows(s) to Top',...
    'Callback',@(~,~) obj.MoveSelectedBlock(-inf,2));
h.UI.hemb = uimenu('parent',h.UI.he,'tag','seqAlignViewer.uimenu.Edit.MoveRowsBottom',...
    'Label','Move Rows(s) to Bottom',...
    'Callback',@(~,~) obj.MoveSelectedBlock(inf,2));
h.UI.her = uimenu('parent',h.UI.he,'tag','seqAlignViewer.uimenu.Edit.RemoveEmptyColumns',...
    'Label','Remove Empty Columns',...
    'Separator','on',...
    'Callback',@(~,~) obj.DeleteColumns);

h.UI.hd = uimenu('parent',h.hFig,'tag','seqAlignViewer.uimenu.Display',...
    'Label','Display',...
    'Callback',@UpdateDisplayMenu);
h.UI.hdb = uimenu('parent',h.UI.hd,'tag','seqAlignViewer.uimenu.Display.Background',...
    'Label','Background',...
    'Checked','on',...
    'Callback',@BackGroundToggleMenu);
h.UI.hdzi = uimenu('parent',h.UI.hd,'tag','seqAlignViewer.uimenu.Display.ZoomIn',...
    'Label','Zoom In',...
    'Separator','on',...
    'Callback',@ZoomIn);
h.UI.hdzo = uimenu('parent',h.UI.hd,'tag','seqAlignViewer.uimenu.Display.ZoomOut',...
    'Label','Zoom Out',...
    'Callback',@ZoomOut);
h.UI.hdr = uimenu('parent',h.UI.hd,'tag','seqAlignViewer.uimenu.Display.ResetFontSize',...
    'Label','Reset to Default Font Size',...
    'Callback',@ZoomDefault);

h.UI.hdc = uimenu('parent',h.UI.hd,'tag','seqAlignViewer.uimenu.Display.ColorSchemes',...
    'Label','Color Schemes',...
    'Separator','on',...
    'Callback',@UpdateDisplayColorMenu);
h.Colors = getSVColorScheme;
for i = 1:numel(h.Colors)
    h.UI.hdca(i) = uimenu('parent',h.UI.hdc,'visible','off',...
        'Label',h.Colors(i).Name,...
        'UserData',h.Colors(i).Type,...
        'Callback',@SetColorScheme);
end

h.UI.hdt = uimenu('parent',h.UI.hd,'tag','seqAlignViewer.uimenu.Display.ViewTree',...
    'Label','View Tree',...
    'Separator','on',...
    'Callback',@UpdateDisplayTreeMenu);
h.UI.hdta = uimenu('parent',h.UI.hdt,'tag','seqAlignViewer.uimenu.Display.ViewTree.All',...
    'Label','All...',...
    'Callback',@(~,~) DisplayTree(ExportSelectedArea(obj,[],[])));
h.UI.hdts = uimenu('parent',h.UI.hdt,'tag','seqAlignViewer.uimenu.Display.ViewTree.Selected',...
    'Label','Selected...',...
    'Callback',@(~,~) DisplayTree(ExportSelectedArea(obj)));

h.UI.hh = uimenu('parent',h.hFig,'tag','seqAlignViewer.uimenu.Help','Label','Help');
uimenu('parent',h.UI.hh,'tag','seqAlignViewer.uimenu.Help.BioinformaticsToolboxHelp',...
    'Label','Bioinformatics Toolbox Help',...
    'Callback',@(~,~) helpview(fullfile(docroot,'bioinfo','bioinfo.map'),'bioinfo_product_page'));
uimenu('parent',h.UI.hh,'tag','seqAlignViewer.uimenu.Help.GettingStartedTutorial',...
    'Label','Getting Started Tutorial',...
    'Callback',@(~,~) helpview(fullfile(docroot,'bioinfo','bioinfo.map'),'multialignviewer_tutorial'));
uimenu('parent',h.UI.hh,'tag','seqAlignViewer.uimenu.Help.Examples',...
    'Label','Examples','Separator','on',...
    'Callback',@(~,~) demo('toolbox','bioinfo'));
tlbx = ver('bioinfo');
mailstr = sprintf('mailto:bioinfo-feedback@mathworks.com?subject=%s%s%s',...
    'Feedback%20for%20Biological%20Sequence%20Alignment%20Viewer',...
    '%20in%20Bioinformatics%20Toolbox%20',tlbx(1).Version);
uimenu('parent',h.UI.hh,'tag','seqAlignViewer.uimenu.Help.SendFeedback',...
    'Label','Send Feedback...','Separator','on',...
    'Callback',@(~,~) web(mailstr));
aboutstr =  sprintf('%s\nVersion %s %s\nCopyright 1993-%s, MathWorks, Inc.', ...
    tlbx.Name, tlbx.Version, tlbx.Release, datestr(tlbx.Date, 10));
uimenu('parent',h.UI.hh,'tag','seqAlignViewer.uimenu.Help.About',...
    'Label','About','Separator','on',...
    'Callback',@(~,~) msgbox(aboutstr,['About ' tlbx.Name],'modal'));

% uicontextmenu
h.UI.hcm = uicontextmenu('parent',h.hFig,'Tag','seqAlignViewer.uicontextmenu',...
    'Callback',@UpdateContextMenu);
uimenu(h.UI.hcm,'Tag','seqAlignViewer.uicontextmenu.Copy',...
    'Label','Copy',...
    'Callback',@Export2Clipboard);
uimenu(h.UI.hcm,'Tag','seqAlignViewer.uicontextmenu.Export2Fasta',...
    'Label','Export to FASTA File...',...
    'Callback',@(~,~) Export2Fasta('selection'));
uimenu(h.UI.hcm,'Tag','seqAlignViewer.uicontextmenu.Export2Workspace',...
    'Label','Export to MATLAB Workspace...',...
    'Callback',@(~,~) Export2Workspace('selection'));
uimenu(h.UI.hcm,'Tag','seqAlignViewer.uicontextmenu.DeleteSequences',...
    'Separator','on',...
    'Label','Delete Sequences',...
    'Callback',@(~,~) obj.DeleteSequences);
uimenu(h.UI.hcm,'Tag','seqAlignViewer.uicontextmenu.SelectAll',...
    'Separator','on',...
    'Label','Select All',...
    'Callback',@(~,~) obj.Select(1));
uimenu(h.UI.hcm,'Tag','seqAlignViewer.uicontextmenu.DeselectAll',...
    'Label','Deselect All',...
    'Callback',@(~,~) obj.Select(0));
uimenu(h.UI.hcm,'Tag','seqAlignViewer.uicontextmenu.DeleteColumns',...
    'Separator','on',...
    'Label','Remove Empty Columns',...
    'Callback',@(~,~) obj.DeleteColumns);

h.Alphabet = alpha;
if strcmp(alpha,'nt')
    h.Scheme = 'NT';
else
    h.Scheme = 'Function';
end

% UISAV initialisieren
obj.Scheme = h.Colors(strcmp({h.Colors.Name},h.Scheme));
obj.Data = data;
obj.ContextMenu = h.UI.hcm;

% when instantiated UISAV defaulted to the right architecture dependent
% font size, we make the client aware of it because is the client that one
% that handles the font size

h.DefaultFontSize = obj.getFontSize; 
h.FontSize = obj.getFontSize;
ZoomDefault

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callbacks for HG based viewer (nested in sav function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Open(~,~)
        fname = uigetfile( ...
            {'*.fasta;*.fa;*.txt', 'FASTA File Formats (*.fasta, *.fa, *.txt)'; ...
            '*.aln;*.msf;*.txt', 'Line Interleaved Formats (*.aln, *.msf, *.txt)';...
            '*.*','All Files (*.*)'},'Open');
        if fname
            [sequences,headers,alphabet,~,ME] = validateAlignment(fname,[],'','');
            if ~isempty(ME)
                errordlg(ME.message,'Open Alignment File Error','modal')
            elseif obj.IsDataEmpty % put data in current viewer
                obj.Data = struct('Sequence',sequences(:),'Header',headers(:)); % put data in current viewer
            else % send the data to a new viewer
                seqalignviewer(struct('Sequence',sequences(:),'Header',headers(:)),'alphabet',alphabet);
            end
        end
    end

    function ImportFromWorkspace(~,~)
        vars = evalin('base','whos');
        vars = vars(strcmp('struct',{vars.class}));
        for iVar = numel(vars):-1:1
            if evalin('base',sprintf('isfield(%s,''Sequence'')',vars(iVar).name))
                vars(iVar).info = sprintf('%-10s %dx%d',vars(iVar).name,vars(iVar).size(1),vars(iVar).size(2));
            else
                vars(iVar)=[];
            end
        end
        if isempty(vars)
            errordlg('Supported variables were not found in workspace',...
                'Error Importing Sequence','modal')
        else
            answer = listdlg('Name','Import from MATLAB Workspace',...
                'PromptString','Select a workspace variable:',...
                'SelectionMode','single',...
                'ListString',{vars.info});
            if ~isempty(answer)
                [sequences,headers,alphabet,~,ME] = validateAlignment(evalin('base',vars(answer).name),[],'','');
                if ~isempty(ME)
                    errordlg(ME.message,'Error Importing Sequence','modal')
                elseif obj.IsDataEmpty % put data in current viewer
                    obj.Data = struct('Sequence',sequences(:),'Header',headers(:)); % put data in current viewer
                else % send the data to a new viewer
                    seqalignviewer(struct('Sequence',sequences(:),'Header',headers(:)),'alphabet',alphabet);
                end
            end
        end
    end

    function Export2Clipboard(~,~)
        output = obj.ExportSelectedArea;
        clipboard('copy',sprintf('%s\n',output.Sequence));
    end

    function Export2Fasta(what)
        switch what
            case 'selection'
                output = obj.ExportSelectedArea;
            case 'consensus'
                output = obj.ExportConsensus;
            otherwise % all
                output = obj.ExportSelectedArea([],[]);
        end
        [fname, pname] = uiputfile( ...
            {'*.fasta;*.fa;*.txt', 'FASTA File Formats (*.fasta, *.fa, *.txt)'; ...
            '*.*','All Files (*.*)'},'Save as FASTA File');
        if fname
            fastawrite(fullfile(pname,fname),output);
        end
    end

    function Export2Workspace(what)
        switch what
            case 'selection'
                output = obj.ExportSelectedArea;
            case 'consensus'
                output = obj.ExportConsensus;
            otherwise % all
                output = obj.ExportSelectedArea([],[]);
        end
        while 1
            varname = char(inputdlg('Input a Variable Name','Export',1));
            if varname
                if isvarname(varname)
                    if ~isempty(evalin('base',sprintf('who(''%s'')',varname)))
                        answer = questdlg(...
                            sprintf('The variable "%s" already exists.\nDo you want to replace the existing variable?',varname), ...
                            'Replace Existing Variable', ...
                            'Yes', 'No', 'Cancel', 'No');
                        switch answer,
                            case 'Yes'
                                assignin('base',varname,output)
                                break;
                            case 'Cancel',
                                break;
                        end  % switch
                    else
                        assignin('base',varname,output)
                        break;
                    end
                else
                    answer = questdlg(...
                        sprintf('"%s" is an invalid variable name!\n"%s" is more appropriate.\n\nwould you like to use use the suggested name?',varname,genvarname(varname)), ...
                        'Error Exporting Sequence', ...
                        'Yes', 'No', 'Cancel', 'Yes');
                    switch answer,
                        case 'Yes',
                            assignin('base',genvarname(varname),output);
                            break;
                        case 'Cancel',
                            break;
                    end % switch
                end
            else
                break
            end
        end %while
    end

    function ZoomOut(~,~)
       fsr = obj.getFontSizeReel;
       ns = max(fsr(fsr<h.FontSize));
       if isempty(ns)
           ns = min(obj.getFontSizeReel);
       end
       obj.FontSize = ns;
       h.FontSize = obj.getFontSize;
    end

    function ZoomDefault(~,~)
        h.FontSize = h.DefaultFontSize;
        obj.FontSize = h.FontSize;
    end

    function ZoomIn(~,~)
        obj.FontSize = h.FontSize + 1;
        h.FontSize = obj.getFontSize;
    end

    function BackGroundToggleToolbar(hObj,~)
        obj.EnableBackground = get(hObj,'state');
    end

    function BackGroundToggleMenu(this,~)
        set(h.UI.htb,'State',char(setdiff({'on','off'},get(this,'Checked'))))
        obj.EnableBackground = get(h.UI.htb,'State');
    end

    function SetColorScheme(this,~)
        h.Scheme = get(this,'Label');
        clr = getSVColorScheme;
        obj.Scheme = clr(strcmp({clr.Name},get(this,'Label')));
    end

    function DisplayTree(seqs)
        try
            dist = seqpdist(seqs,'Alphabet',h.Alphabet);
            tree = seqneighjoin(dist,'equivar',seqs);
            view(tree);
        catch ME
            errordlg(ME.message,'Error Openning Tree Window','modal')
        end
    end

    function UpdateFileMenu(~,~)
        if obj.IsDataEmpty
            ee = 'off';
        else
            ee = 'on';
        end
        set([h.UI.hfe1,h.UI.hfe2,h.UI.hfe3,h.UI.hfe4,h.UI.hfp],'Enable',ee)
    end

    function UpdateEditMenu(~,~)
        if obj.IsAreaSelected
            ee = 'on';
        else
            ee = 'off';
        end
        set([h.UI.hed h.UI.hemu h.UI.hemd h.UI.hemt h.UI.hemb] ,'Enable',ee)
        if obj.IsDataEmpty
            ee = 'off';
        else
            ee = 'on';
        end
        set([h.UI.hesa h.UI.heda h.UI.hec h.UI.her] ,'Enable',ee)
    end

    function UpdateContextMenu(this,~)
        if obj.IsAreaSelected
            ee = 'on';
        else
            ee = 'off';
        end        
        set(findall(this,'Tag','seqAlignViewer.uicontextmenu.DeleteSequences'),'Enable',ee)
    end

    function UpdateDisplayMenu(~,~)
        set(h.UI.hdb,'Checked',get(h.UI.htb,'state'))
        if obj.IsDataEmpty
            ee = 'off';
        else
            ee = 'on';
        end
        set([h.UI.hdb h.UI.hdzi h.UI.hdzo h.UI.hdr h.UI.hdc h.UI.hdt],'Enable',ee)
    end

    function UpdateDisplayTreeMenu(~,~)
        if obj.IsAreaSelected
            set(h.UI.hdts,'Enable','on')
        else
            set(h.UI.hdts,'Enable','off')
        end
    end

    function UpdateDisplayColorMenu(~,~)
        set(h.UI.hdca,'Visible','off','Checked','off')
        for idx = 1:numel(h.UI.hdca)
            if strcmp(h.Scheme,get(h.UI.hdca(idx),'Label'))
                set(h.UI.hdca(idx),'Checked','on')
            end
            if strcmp(h.Alphabet,get(h.UI.hdca(idx),'UserData'))
                set(h.UI.hdca(idx),'Visible','on')
            end
        end
    end

end

function color = getSVColorScheme
xDoc = xmlread('msacolorschemes.xml'); %toolbox\bioinfo\bioinfo\@seqviewer\msacolorschemes.xml
color = struct('Name',[],'Type',[],'Symbols',[],'Colors',{});
schemes = xDoc.getElementsByTagName('Scheme');
for n = 0:schemes.getLength-1
    name = schemes.item(n).getElementsByTagName('Name');
    color(n+1).Name = char(name.item(0).getFirstChild.getData);
    type = schemes.item(n).getElementsByTagName('Type');
    color(n+1).Type = lower(char(type.item(0).getFirstChild.getData));
    clr = schemes.item(n).getElementsByTagName('option');
    for nClr = 0:clr.getLength-1
        color(n+1).Symbols(end+1) = char(clr.item(nClr).getAttribute('name'));
        color(n+1).Colors(end+1,:) = sscanf(char(clr.item(nClr).getAttribute('value')),'%f %f %f')'/255;
    end
end
end

