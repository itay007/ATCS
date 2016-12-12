function n = ndims(a)
%NDIMS Number of dimensions of a DataMatrix object.
% 
%   N = NDIMS(A) returns the number of dimensions in the DataMatrix A.
%   The number of dimensions in a DataMatrix object is always 2.  
%
%   See also DATAMATRIX/SIZE.

%   Copyright 2008 The MathWorks, Inc.


n = a.NDims;
end % ndims method
