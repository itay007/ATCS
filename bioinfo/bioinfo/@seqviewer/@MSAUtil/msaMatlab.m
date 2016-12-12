function [result, varargout] = msaMatlab(methodName, varargin) %#ok
%MSAMATLAB Bridge between SEQALIGNVIEWER Java components and MATLAB functions.
%
%   RES = MSAMATLAB(METHODNAME) can be called by SEQALIGNVIEWER Java
%   components by pass in MATLAB function name and input variables for the
%   function. The results returned from called MATLAB functions are
%   returned to the Java component.

%   Copyright 2005-2012 The MathWorks, Inc.


result = [];

try
    switch(methodName)
        case{'get_msaviewer',...
                'jvm_available'}
            result = feval(methodName);
        case{'bioinfo_version'}
            feval(methodName);
        otherwise
            result = feval(methodName, varargin{:});
    end
catch theException
    error(message('bioinfo:MSAViewer:msaMatlabERROR', methodName, theException.message));
end

return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = jvm_available

persistent sRlt;

if isempty(sRlt)
    sRlt=usejava('jvm') & usejava('awt') & usejava('swing');
end
result = sRlt;
return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function viewer = get_msaviewer %#ok

persistent sJavaErrorMessage

viewer = [];
if ~jvm_available
    error(message('bioinfo:MSAViewer:JavaComponentMissing'));
end

if isempty(sJavaErrorMessage)
    try
        viewer = com.mathworks.toolbox.bioinfo.sequence.msaviewer.MSAViewerApi.viewerHandle();
    catch theException
        viewer = [];
        sJavaErrorMessage = sprintf('bioinfo:MSAViewer:FailToOpenError %s', theException.message);
        warning('bioinfo:MSAViewer:FailToOpenError',theException.message);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = manage_msaviewer(varargin)

persistent sJavaErrorOccurred;
result = [];
if ~jvm_available
    return;
end

hViewer = com.mathworks.toolbox.bioinfo.sequence.msaviewer.MSAViewerApi.viewerHandle();
result = hViewer;

if isempty(sJavaErrorOccurred)
    try
        if nargin > 0
            Action = varargin{1};
            switch(Action)
                case 'status'
                    if ~isempty(hViewer)
                        result = 1;
                    else
                        result = 0;
                    end
                case 'state'
                    if ~isempty(hViewer)
                        result = ~hViewer.isEmpty;
                    else
                        result = 0;
                    end

                case 'close'
                    if ~isempty(hViewer)
                        awtinvoke(hViewer, 'viewerTerminate()');
                        hViewer = []; %#ok
                        clear hViewer;
                    end
            end
        end
    catch theException
        sJavaErrorOccurred = sprintf('bioinfo:MSAViewer:JavaErrorOccurred %s', theException.message);
        warning('bioinfo:MSAViewer:JavaErrorOccurred',theException.message);
    end

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = create_ui(varargin) %#ok

if ~seqviewer.MSAUtil.msaMatlab('jvm_available');
    error(message('bioinfo:MSAViewer:JVMMissing'));
end

hViewer =  seqviewer.MSAUtil.msaMatlab('manage_msaviewer');
if(isempty(hViewer))
    error(message('bioinfo:MSAViewer:CanNotInitialize'));
end

if nargin == 0 % open blank
    hViewer.blankOpen();
    result = hViewer;
    return;
end

if nargin == 1 % if errored and viewer is open
    errmsg = varargin{1};
    hViewer.showErrorDialog(errmsg);
    result = hViewer;
    return;
end

% no error, open viewer
title = varargin{1};
headers = varargin{2};
sequences = varargin{3};
alphabet = varargin{4};

[consensus, scores]=getconservationscore(sequences, alphabet);

opened = hViewer.alignmentOpen(alphabet,...
    title,...
    headers,...
    sequences,...
    consensus,...
    scores);

if ~opened
    error(message('bioinfo:MSAViewer:FailToOpenViewerPane'));
end


result = hViewer;
return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = get_consensus(varargin) %#ok
type = varargin{1};
seqs = varargin{2}; % Java string array

try
    [c,s]=getconservationscore(seqs, type);
    result =  struct('Consensus', c, 'Scores', s);
catch theException
    result = {theException.message};
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = save_fasta(varargin) %#ok
file = varargin{1};
seqStruct =getAlignmentStruct(varargin{2});
try
    fastawrite(file,seqStruct);
catch theException
    if(nargout == 1)
        varargout{1} = theException.message;
    end

    return;
end

if(nargout == 1)
    varargout{1} = '';
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = open_file(varargin)%#ok
path = varargin{1};
try
    msg = seqalignviewer( path, 'FROMVIEWER', true, 'R2012b', true);
catch theException
    msg = theException.message;
    return;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = export_alignment(varargin)%#ok
% Assigns variable_data to variable_name in the base workspace

% Grab all the workspace variables and store their names in a cell
% array.
msg = '';

variable_name = varargin{1};
variable_data = getAlignmentStruct(varargin{2})';
check_variable = varargin{3};

ws_vars = evalin('base','whos');
[ws_var_names{1:length(ws_vars)}] = deal(ws_vars.name);

valid_var_name = genvarname(variable_name);

var_already_exists = any(strcmp(valid_var_name, ws_var_names),2);
user_spec_var_name_changed = ~strcmpi(variable_name, valid_var_name);

if user_spec_var_name_changed
    msg = [valid_var_name '_valid'];
    return;
end

if ~isempty(var_already_exists) && var_already_exists && check_variable
    msg = '_exist';
    return;
end

assignin('base',valid_var_name, variable_data);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function msg = import_alignment(var_name)%#ok
% Import variable_name from the base workspace

msg = '';
alignment = [];
if ischar(var_name)
    try
        alignment = evalin('base',var_name);
    catch theException
        msg = theException.message;
        return;
    end
end

if ~isempty(alignment)
    try
        msg = seqalignviewer( alignment, 'varname', var_name, 'fromviewer', true, 'R2012b', true);
    catch theException
        msg = theException.message;
        return;
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = view_tree(varargin)%#ok
msg = '';
type = varargin{1};
seqs =getAlignmentStruct(varargin{2});

try
    dist = seqpdist( seqs, 'Alphabet', type);
    % tree = seqlinkage( dist, 'average', seqs ); %UPGMA method
    tree = seqneighjoin( dist, 'equivar', seqs );
    view(tree);
catch theException
    msg = theException.message;
    return;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alignStruct = getAlignmentStruct(alignCell)
% Convert cell array of alignment from java to a structure

alignStruct = cell2struct(vertcat(alignCell{:})', {'Header','Sequence'});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute conservation score from seqconsensus.
% conservationScore = max(significantSMScore - score,0);
% significantSMScore  = -13.34*(N^-.76) + 11.50;
% conservationScore = 0 if consensus is a gap.

function [c, s] = getconservationscore(seqs, type)
[c, s] = seqconsensus( seqs, 'Alphabet', type, 'gaps', 'all');
s = max(11.5 - s,0);
gapindex = strfind(c, '-');
s(gapindex) = 0;
end
%--------------------------------------------
% Return toolbox version
function bioinfo_version %#ok
tlbx = ver('bioinfo');

mailstr = ['mailto:bioinfo-feedback@mathworks.com?subject=',...
    'Feedback%20for%20Biological%20Sequence%20Alignment%20Viewer%20in%20Bioinformatics',...
    '%20Toolbox%20',tlbx(1).Version];
web(mailstr)
end
