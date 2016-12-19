function m = abs(a)
%ABS Overload absolute value for DataMatrix object.
% 
%   B = ABS(A) returns the absolute values of the elements of DataMatrix
%   object A. Output B is a DataMatrix object. The size, row names and
%   column names of B are the same as A.  
% 
%   See also DATAMATRIX/SIGN.

%   Copyright 2008 The MathWorks, Inc. 


m = a;
m.Matrix = abs(a.Matrix);
end %DataMatrix/abs
