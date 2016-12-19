function len = length(obj)
%LENGTH Length of a DataMatrix.
% 
%   N = LENGTH(A) returns the number of rows in the DataMatrix A.
%   LENGTH is equivalent to SIZE(A,1).
%  
%   See also DATMATRIX/SIZE.

%   Copyright 2008 The MathWorks, Inc. 


len = length(obj.Matrix);
