function bioerrorrethrow(matlabfile,ME)
%BIOERRORRETHROW reissues a function exception with a new identifier. 
% 
%   BIOERRORRETHROW(MATLABFILE,ME) creates a new MException using the same
%   message as ME but with the identifier bioinfo:MATLABFILE:identifier.
%   Then issues the exception as if from the calling function.
%
%   BIOERRORRETHROW is intended to be used when an error occurs within a
%   try-catch block in a Bioinformatics Toolbox function.

% Copyright 2008-2012 The MathWorks, Inc.

x = MException(sprintf('bioinfo:%s:%s',matlabfile,regexp(ME.identifier,'[^:]*$','match','once')),'%s',ME.message);  
x.throwAsCaller;
end 
