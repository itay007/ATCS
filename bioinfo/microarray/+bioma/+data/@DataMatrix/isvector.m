function t = isvector(a)
%ISVECTOR True if DataMatrix object is a vector.
% 
%   TF = ISVECTOR(A) returns true (1) if the DataMatrix object A is a 1-by-N
%   or N-by-1 vector, where N >= 0, and false (0) otherwise.
%
%   See also DATAMATRIX/ISEMPTY, DATAMATRIX/ISSCALAR, DATAMATRIX/SIZE.

%   Copyright 2008 The MathWorks, Inc. 


t = isvector(a.Matrix);
end %DataMatrix/isvector
