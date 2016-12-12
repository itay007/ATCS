function evalrasmolscript(hviewer, cmd)
%EVALRASMOLSCRIPT sends Rasmol scripts to a molecule viewer.
%
%   EVALRASMOLSCRIPT(H, S) sends Rasmol script S to a molecule viewer with
%   handle H. S can be a string, a file name of a Rasmol script file, or a
%   cell array of strings to be concatenated.
%
%   Examples:
%
%       h = molviewer('aspirin.mol'); 
%       %Wait until the model finished loading 
%       evalrasmolscript(h, 'spin on')
%
%   See also MOLVIEWER, MOLVIEWERDEMO.

% Copyright 2006-2008 The MathWorks, Inc.


bioinfochecknargin(nargin,2,mfilename);

if isempty(hviewer) || ~ishandle(hviewer)
    error(message('bioinfo:evalrasmolscript:InvalidViewerHandle'));
end

if ~isappdata(hviewer, 'BioMolViewer')
    error(message('bioinfo:evalrasmolscript:NotAViewerHandle'));
end

isFileFlag = false;
script = '';
if (ischar(cmd) && size(cmd,1) == 1)
    if (size(fullfile(pwd,cmd), 2) < 1000) && (exist(cmd,'file') || exist(fullfile(pwd,cmd),'file'))
        isFileFlag = true;  
        script = handlefilename(cmd);
    else
        script = cmd;
    end  
elseif iscellstr(cmd)
    for i = 1:length(cmd)
        script = strcat(script, cmd{i}, ';');
    end
else
    error(message('bioinfo:evalrasmolscript:InvalidScriptType'))
end

appdata = getappdata(hviewer, 'BioMolViewer');

if(appdata.isblankflag)
    hviewer = molviewer('resetFigureTool', hviewer);
    appdata = getappdata(hviewer, 'BioMolViewer');
end

if isempty(appdata.viewer)
    error(message('bioinfo:evalrasmolscript:MolViewerNotExist'))
end

if isFileFlag
   awtinvoke(appdata.viewer, 'evalMVFile(Ljava/lang/String;)', script);
else
   awtinvoke(appdata.viewer, 'evalMVScriptT(Ljava/lang/String;)', script);
end

%-------------------------------------------
function fullname = handlefilename(filename)
    wfile = which(filename);
    if ~isempty(wfile)
        filename = wfile;
    end

    [thePath, theName, theExt] = fileparts(filename);
    if isempty(thePath)
        thePath = pwd;
    end

    fullname = fullfile(thePath, [theName theExt]);
