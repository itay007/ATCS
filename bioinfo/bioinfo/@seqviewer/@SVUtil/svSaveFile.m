function varargout = svSaveFile(varargin)
%SVSAVEFILE Call by SEQVIEWER Java components to save sequences to file.
%
%   RES = SVSAVEFILE(VARARGIN) can be called by SEQVIEWER Java components to
%   save sequence to file. The results returned from called MATLAB
%   functions are returned to the Java component.

%   Copyright 2005-2012 The MathWorks, Inc.


file = varargin{1};
header = varargin{2};
sequence = varargin{3};
 
try
    if isempty(header)
       feval('fastawrite', file, sequence);
    else
       data = struct('Sequence', sequence, 'Header', header); 
       feval('fastawrite', file, data);
    end
catch theException
    if(nargout == 1)
        varargout{1} = theException.message;
    end
    return;
end

if(nargout == 1)
    varargout{1} = '';
end
