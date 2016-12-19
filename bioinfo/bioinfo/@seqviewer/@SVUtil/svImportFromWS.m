function msg = svImportFromWS(var_name)
%SVIMPORTFROMWS Call by SEQVIEWER Java components to import variables from workspace.
%
%   RES = SVIMPORTFROMWS(VARNAME) can be called by SEQVIEWER Java components
%   to import variables, VARNAME, from the base workspace. The results
%   returned from called MATLAB functions are returned to the Java
%   component.

%   Copyright 2005-2012 The MathWorks, Inc.


msg = '';
sequence = [];
if ischar(var_name)
    try
        sequence = evalin('base',var_name);
    catch theException
        msg = theException.message;
        return;
    end
end

if ~isempty(sequence)
    try
        feval('seqviewer', sequence, var_name);
    catch theException
        msg = theException.message;
        return;
    end
end

 
