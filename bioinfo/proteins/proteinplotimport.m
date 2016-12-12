function varargout = proteinplotimport(varargin)
% PROTEINPLOTIMPORT Imports sequences for PROTEINPLOT

%   Copyright 2003-2009 The MathWorks, Inc.


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @proteinplotimport_OpeningFcn, ...
    'gui_OutputFcn',  @proteinplotimport_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% ----------------------------------------------------------------
function proteinplotimport_OpeningFcn(hObject, eventdata, handles, varargin) %#ok


% Choose default command line output for proteinplotimport
handles.output = hObject;

% field of handles for proteinplot
handles.proteinplot = varargin{1};

% initialize field to hold sequence data
handles.sequence = '';

% center the dialog
centerdlg(handles.proteinplot.proteinplotfig,handles.ppimport);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes proteinplotimport wait for user response (see UIRESUME)
% uiwait(handles.ppimport);


% ----------------------------------------------------------------
function varargout = proteinplotimport_OutputFcn(hObject, eventdata, handles) %#ok
varargout{1} = handles.ppimport;

% ----------------------------------------------------------------
function importfrom_text_CreateFcn(hObject, eventdata, handles) %#ok

set_os_text(hObject)

oldunits = get(hObject,'Units');
set(hObject,'Units','characters');
p = get(hObject,'Position');

set(hObject,'Position',[p(1) p(2) 13 p(4)])
set(hObject,'Units',oldunits)

% ----------------------------------------------------------------
function importfrom_popup_CreateFcn(hObject, eventdata, handles) %#ok
 
set_os_text(hObject)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% ----------------------------------------------------------------
function importfrom_popup_Callback(hObject, eventdata, handles) %#ok
v = get(hObject,'Value');
oldunits = get(handles.misc_static,'Units');
set(handles.misc_static,'Units','characters');
p = get(handles.misc_static,'Position');
if v == 1,
    set([handles.misc_text, handles.misc_button],'Visible','off')
    set(handles.misc_list,'Visible','on')
    set(handles.misc_static,'String','Variable','Position',[p(1) p(2) 10 p(4)])
    set(handles.example_txt,'Visible','off')
    populate_list(handles.misc_list);
elseif v == 2 || v == 3 || v == 4,
    set(handles.misc_list,'Visible','off')
    set(handles.misc_text,'Visible','on')
    set(handles.misc_button,'Visible','on','String','Browse...')
    set(handles.misc_static,'String','File','Position',[p(1) p(2) 5 p(4)])
    set(handles.example_txt,'Visible','off')   
else
    set(handles.misc_list,'Visible','off')
    set(handles.misc_text,'Visible','on')
    set(handles.misc_button,'Visible','off')
    set(handles.misc_static,'String','Access ID','Position',[p(1) p(2) 11 p(4)])
    set(handles.example_txt,'Visible','on')
end

set(handles.misc_static,'Units',oldunits)


% ----------------------------------------------------------------
function ok_button_Callback(hObject, eventdata, handles) %#ok
importdata(hObject,eventdata,handles);

% ----------------------------------------------------------------
function ok_button_KeyPressFcn(hObject, eventdata, handles) %#ok
curkey = get(handles.ppimport,'CurrentKey');
if ~any(strmatch(curkey,{'return';'space'}))
    return
end

importdata(hObject,eventdata,handles);

% ----------------------------------------------------------------
function importdata(hObject,eventdata,handles) %#ok

v = get(handles.importfrom_popup,'Value');

err = '';

if v == 1, % workspace
    val = get(handles.misc_list,'Value');
    names = get(handles.misc_list,'String');
    if ~isempty(names)
        handles.sequence = evalin('base',char(names(val)));
    else
        handles.sequence = '';
    end
elseif v == 2, % text file
    fullname = get(handles.misc_text,'String');
    try
        theSeq = bioinfoprivate.seqread(fullname);
        ftext = theSeq.Sequence;
        if all(isletter(ftext)) % see if it's a plain sequence of characters
            handles.sequence = ftext;
        else
            err = sprintf('%s contains invalid characters for a plain text sequence: %s',fullname, unique(ftext(~isletter(ftext)))');
        end
    catch theException
        err = char(strread(theException.message,'%*s%s','delimiter','\n'));
    end
elseif v == 3, % FASTA file
    fullname = get(handles.misc_text,'String');
    try
        [~,seq] = fastaread(fullname);
        if iscell(seq)
            handles.sequence = seq{1};
            if numel(seq) > 1
                waitfor(warndlg(sprintf('Multiple sequences found in %s. Using first sequence.',fullname),'Multiple Sequences Found','modal'));
            end
        else 
            handles.sequence = seq;
        end        
    catch allExceptions %#ok<NASGU>
        err = ['Cannot use FASTAREAD to import ' fullname];
    end
elseif v == 4, % GenPept file
    fullname = get(handles.misc_text,'String');
    try
        gpstruct = genpeptread(fullname);
        handles.sequence = gpstruct.Sequence;    
    catch allExceptions %#ok<NASGU>
        err = ['Cannot use GENPEPTREAD to import ' fullname];
    end
else % GenPept database
    % search for the given access id
    try
        handles.sequence = getgenpept(char(get(handles.misc_text,'String')),'sequenceonly',true);
    catch theException
        err = theException.message;
    end
end    

if ~isempty(err)
    errordlg(err,'Error');
    return
end

handles.sequence = cleansequence(handles.sequence);

set(handles.proteinplot.sequence,'String',handles.sequence)
close(handles.ppimport);

% ----------------------------------------------------------- -----
function cancel_button_Callback(hObject, eventdata, handles) %#ok
close(handles.ppimport);

% ----------------------------------------------------------------
function cancel_button_KeyPressFcn(hObject, eventdata, handles) %#ok
curkey = get(handles.ppimport,'CurrentKey');
if ~any(strmatch(curkey,{'return';'space'}))
    return
end
close(handles.ppimport);


% ----------------------------------------------------------------
function misc_static_CreateFcn(hObject, eventdata, handles) %#ok

set_os_text(hObject)

oldunits = get(hObject,'Units');
set(hObject,'Units','characters');
p = get(hObject,'Position');
set(hObject,'Position',[p(1) p(2) 10 p(4)])
set(hObject,'Units',oldunits)


% ----------------------------------------------------------------
function misc_list_CreateFcn(hObject, eventdata, handles) %#ok

set_os_text(hObject)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
populate_list(hObject)


% ----------------------------------------------------------------
function misc_list_Callback(hObject, eventdata, handles) %#ok


% ----------------------------------------------------------------
function misc_button_Callback(hObject, eventdata, handles) %#ok

v = get(handles.importfrom_popup,'Value');

if v == 2 || v == 3 || v == 4
    [filename,pathname] = uigetfile({'*.*','All files (*.*)'},'Select file');
    if ~filename
        return
    end    
    set(handles.misc_text,'String',fullfile(pathname,filename));
end

% -------------------------------------------------------------- --
function misc_text_CreateFcn(hObject, eventdata, handles) %#ok

set_os_text(hObject)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% ----------------------------------------------------------------
function ok_button_CreateFcn(hObject, eventdata, handles) %#ok
set_os_text(hObject)

% ----------------------------------------------------------------
function cancel_button_CreateFcn(hObject, eventdata, handles) %#ok
set_os_text(hObject)

% ----------------------------------------------------------------
function example_txt_DeleteFcn(hObject, eventdata, handles) %#ok
set_os_text(hObject)

% ----------------------------------------------------------------
function misc_button_CreateFcn(hObject, eventdata, handles) %#ok
set_os_text(hObject)

% ----------------------------------------------------------------
function example_txt_CreateFcn(hObject, eventdata, handles) %#ok
set_os_text(hObject)


% ----------------------------------------------------------------
function misc_text_Callback(hObject, eventdata, handles) %#ok


% ----------------------------------------------------------------
function populate_list(hObject)
vars = evalin('base','whos');
varsclass = {vars.class}';
ischarvar = strcmp(varsclass,'char');
varnames = {vars(ischarvar).name}';
isstructvar = strcmp(varsclass,'struct');
for n = 1:length(isstructvar)
    if isstructvar(n) && evalin('base',['isfield(' vars(n).name ',''Sequence'')' ])
        varnames = [varnames;{[vars(n).name '.Sequence']}]; %#ok<AGROW>
    end
end
set(hObject,'String',varnames);

% ----------------------------------------------------------------
function set_os_text(hObject)
if ispc
    set(hObject,'FontName','MS Sans Serif','FontSize',8)
else
    set(hObject,'FontName','Helvetica','FontSize',10)
end


% ----------------------------------------------------------------
function ppimport_CloseRequestFcn(hObject, eventdata, handles) %#ok
delete(hObject)
