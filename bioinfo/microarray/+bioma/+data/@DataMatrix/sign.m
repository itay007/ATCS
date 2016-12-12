function b = sign(a)
%SIGN  Overload signum function for DataMatrix object.
% 
%   B = SIGN(A) returns the signum of the elements of DataMatrix object A.
%   Output B is a DataMatrix object. The size, row names and column names
%   of B are the same as A.
% 
%   See also SIGN, DATAMATRIX/ABS.

%   Copyright 2008 The MathWorks, Inc. 

b= a;
b.Matrix = sign(a.Matrix);
end %DataMatrix/sign
