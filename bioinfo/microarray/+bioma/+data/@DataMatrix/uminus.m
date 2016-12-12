function b = uminus(a)
%UMINUS Unary minus for a DataMatrix object.
% 
%   B = UMINUS(A) negatives the elements of DataMatrix A. B is a
%   DataMatrix object with the same row names and columns names as A.
% 
%   B = UMINUS(A) is called for the syntax '-A' when A is a DataMatrix
%   object.
% 
%   See also DATAMATRIX/UPLUS.

%   Copyright 2008 The MathWorks, Inc. 


b = a;
b.Matrix = uminus(a.Matrix);
end %DataMatrix/uminus
