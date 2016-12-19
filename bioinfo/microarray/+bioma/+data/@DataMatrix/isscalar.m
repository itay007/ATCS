function t = isscalar(a)
%ISSCALAR True if DataMatrix is a scalar.
% 
%   TF = ISSCALAR(A) returns true (1) if the DataMatrix object A is a
%   1-by-1 matrix, and false (0) otherwise.
%
%   See also DATAMATRIX/ISEMPTY, DATAMATRIX/ISVECTOR, DATAMATRIX/SIZE.

%   Copyright 2008 The MathWorks, Inc. 


t = isscalar(a.Matrix);
end %DataMatrix/isscalar
