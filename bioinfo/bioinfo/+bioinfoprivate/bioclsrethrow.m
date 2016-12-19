function bioclsrethrow(clsname,methodname,ME)
%BIOCLSRETHROW reissues a method exception with a new identifier
% 
%   BIOCLSRETHROW(CLSNAME,METHODNAME,ME) creates a new MException using the
%   same message as ME but with the identifier
%   bioinfo:CLSNAME:METHODNAME:identifier. Then issues the exception as if
%   from the calling method. 
%
%   BIOCLSRETHROW is intended to be used when an error occurs within a
%   try-catch block in a Bioinformatics Toolbox method.

% Copyright 2008-2012 The MathWorks, Inc.

x = MException(sprintf('bioinfo:%s:%s:%s',clsname,methodname,regexp(ME.identifier,'[^:]*$','match','once')),'%s',ME.message);  
x.throwAsCaller;
end 
