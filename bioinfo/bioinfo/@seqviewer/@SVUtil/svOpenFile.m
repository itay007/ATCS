function varargout = svOpenFile(varargin)
%SVOPENFILE Call by SEQVIEWER Java components to open sequence files.
%
%   RES = SVOPENFILE(VARARGIN) can be called by SEQVIEWER Java components to
%   open sequence files. The results returned from called MATLAB functions
%   are returned to the Java component.

%   Copyright 2005-2012 The MathWorks, Inc.


file = varargin{1};

try
    sequence = feval('genbankread', file);
catch theException
    if(nargout == 1)
        varargout{1} = theException.message;
    end
    return;
end

seqviewer.SVUtil.svMatlab('create_ui', sequence);

if(nargout == 1)
    varargout{1} = '';
end
