function [result, varargout] = svMatlab(methodName, varargin) %#ok
%SVMATLAB Bridge between SEQVIEWER Java components and MATLAB functions.
%
%   RES = SVMATLAB(METHODNAME) can be called by SEQVIEWER Java components by
%   pass in MATLAB function name and input variables for the function. The
%   results returned from called MATLAB functions are returned to the Java
%   component.

%   Copyright 2005-2012 The MathWorks, Inc.


result = [];

try
    switch(methodName)
        case{'get_seqviewer',...
                'jvm_available'}
            result = feval(methodName);
        case{'bioinfo_version'}
            feval(methodName);
        otherwise
            result = feval(methodName, varargin{:});
    end
catch theException
    error(message('bioinfo:SequenceTool:svMatlabERROR', methodName, theException.message));
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
function viewer = get_seqviewer

persistent sJavaErrorMessage

viewer = [];
if ~jvm_available
    error(message('bioinfo:SequenceTool:JavaComponentMissing'));
end

if isempty(sJavaErrorMessage)
    try
        viewer = com.mathworks.toolbox.bioinfo.sequence.viewer.SequenceViewerApi.viewerHandle();
    catch theException
        viewer = [];
        sJavaErrorMessage = sprintf('bioinfo:SequenceTool:FailToOpenError %s', theException.message);

        warning('bioinfo:SequenceTool:FailToOpenError',theException.message);
    end
end

return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = manage_bioinfodesktop(varargin)

persistent sJavaErrorOccurred;
result = [];
if ~jvm_available
    return;
end

hViewer = com.mathworks.toolbox.bioinfo.sequence.viewer.SequenceViewerApi.viewerHandle();
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
                        result = ~awtinvoke(hViewer, 'isEmpty()');
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
        sJavaErrorOccurred = sprintf('bioinfo:SequenceTool:JavaErrorOccurred %s', theException.message);

        warning('bioinfo:SequenceTool:JavaErrorOccurred',theException.message);
    end

end

return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = create_ui(varargin) %#ok

seqStr = [];
seqFeatures = [];
seqBaseCounts = [];
seqErr = [];

if ~seqviewer.SVUtil.svMatlab('jvm_available');
    error(message('bioinfo:SequenceTool:JavaComponentMissing'));
end

hViewer =  seqviewer.SVUtil.svMatlab('manage_bioinfodesktop');
if(isempty(hViewer))
    error(message('bioinfo:SequenceTool:CanNotInitialize'));
end

if nargin == 1
    sequence = varargin{1};
    if isstruct(sequence)
        seqType = sequence.Alpha;
        seqId = sequence.LocusName;
        seqAccession = sequence.Accession;
        seqStr = sequence.Sequence;
        seqDef = sequence.Definition;
        seqFeatures = sequence.Features;
        seqComment = sequence.PreambleText;
        [seqErr, seqBaseCounts] = getStats(seqStr, seqType);
        if strcmp(sequence.FileType(end-3:end), 'file') == 1
            seqId = sequence.OriginalData;
        end

        if size(seqDef, 1) > 1
            tmpDef = seqDef(1, :);
            for i = 2:size(seqDef)
                tmpDef = [tmpDef, ' ', seqDef(i,:)]; %#ok
            end
            seqDef = tmpDef;
        end
    elseif ischar(sequence)
        seqErr = sequence;
    end
end

%If the viewer is already open and SEQVIEWER pass in an error message
%the viewer shoul display the error dialog.
if ~isempty(seqErr) && isempty(seqBaseCounts)
    hViewer.showErrorDialog(seqErr);
    seqErr = []; %#ok
    result = hViewer;
    return;
elseif ~isempty(seqErr) && ~isempty(seqBaseCounts)
    hViewer.showWarningDialog(seqErr);
    seqErr = []; %#ok
end

if hViewer.sequenceIsOpen(seqId)
    % Already open. Bring it to front.
    awtinvoke(hViewer, 'sequenceToFront(Ljava/lang/String;Z)', seqId,true)
else
    if isempty(seqFeatures)
        opened = hViewer.sequenceOpen(seqType,...
            seqId, ...
            seqAccession,...
            seqStr,...
            seqDef,...
            seqComment,...
            seqBaseCounts);
    else
        if strcmpi(seqType, 'NT')
            cds = featuresparse(seqFeatures, 'feature', 'CDS');
            opened = hViewer.sequenceOpen(seqId, ...
                seqAccession,...
                seqStr,...
                seqDef,...
                cds,...
                seqFeatures,...
                seqComment,...
                seqBaseCounts);
        else
            opened = hViewer.sequenceOpen(seqType,...
                seqId,...
                seqAccession,...
                seqStr,...
                seqDef,...
                seqFeatures,...
                seqComment,...
                seqBaseCounts);
        end
    end

    if ~opened
        error(message('bioinfo:SequenceTool:FailToOpenViewerPane'));
    end
end

result = hViewer;
return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = close_ui(varargin) %#ok
result = [];
seqId = '';
if nargin == 1
    seqId = varargin{:};
end
hviewer = seqviewer.SVUtil.svMatlab('manage_bioinfodesktop');

openflag = awtinvoke(hViewer, 'sequenceIsOpen(Ljava/lang/String;)', seqId);

if ~isempty(hviewer) && openflag
    awtinvoke(hViewer, 'sequenceClose(Ljava/lang/String;)', seqId);
end

hviewer = []; %#ok
clear hviewer;

return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = help_demos %#ok
result =[];
try
    demo('toolbox', 'bioinfomatics');
catch allExceptions %#ok<NASGU>
    return;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function titleString = create_title_string(seqId) %#ok
% the title string

titleString = ['Sequence: ' seqId];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [msg, result] = getStats(seqStr, type)
warnState = warning;
warning('off', 'bioinfo:basecount:UnknownSymbols')
warning('off', 'bioinfo:aacount:UnknownSymbols')
lastwarn('')
try
    if strcmpi(type, 'NT') %NT
        result = basecount(seqStr,'AMBIGUOUS', 'Individual');
    elseif strcmpi(type, 'AMINO') % AA
        result = aacount(seqStr, 'AMBIGUOUS','Bundle');
    end
    msg = lastwarn;
    if ~isempty(msg)
        tmpMsg = strtok(msg, '.');
        msg = sprintf('%s. %s', tmpMsg, 'These will not be included in the counts.');
    end 
    warning(warnState);
catch ME
    msg = ME.message;
    result = [];
end
end % getStats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = open_proteinplot(sequence) %#ok
try
    result = proteinplot(sequence);
catch theException
    result = theException.message;
end
end % end of open_proteinplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = get_dnastrings(varargin) %#ok
funcname = varargin{1};
seqstr = varargin{2};
try
    result =  feval(funcname, seqstr);
catch theException
    result = theException.message;
end
end % end of get_dnastrings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = get_fulltranslation(varargin) %#ok
seqstr = varargin{1};
code = varargin{2};
try
    result =  nt2aa( seqstr, 'frame', 'all', 'geneticcode', code, 'ACGTOnly', false);
catch theException
    result = {theException.message};
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = get_orf(varargin) %#ok
seqstr = varargin{1};
code = varargin{2};
try
    result =  seqshoworfs( seqstr, 'frame', 'all','geneticcode', code,...
        'NoDisplay', true, 'translate', true);
catch theException
    result = {theException.message};
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = get_word(varargin) %#ok
seqstr = varargin{1};
word = varargin{2};
exact = varargin{3};
type = varargin{4};
try
    result =  seqshowwords( seqstr, word, 'Exact', exact, 'Alphabet', type, 'NoDisplay', true);
catch theException
    result = {theException.message};
end
end

function result = get_basecount(varargin) %#ok
seqstr = varargin{1};
warnState = warning;
warning('off', 'bioinfo:basecount:UnknownSymbols')
lastwarn('')
try
    result =  basecount( seqstr, 'Ambiguous', 'Individual');
catch theException
    result = {theException.message};
end
warning(warnState);
end

% Return toolbox version
function bioinfo_version 
tlbx = ver('bioinfo');

mailstr = ['mailto:bioinfo-feedback@mathworks.com?subject=',...
    'Feedback%20for%20Biological%20Sequence%20Viewer%20in%20Bioinformatics',...
    '%20Toolbox%20',tlbx(1).Version];
web(mailstr)
end

