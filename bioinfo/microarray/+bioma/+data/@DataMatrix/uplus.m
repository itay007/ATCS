function b = uplus(a)
%UPLUS Unary plus for a DataMatrix object.
% 
%   B = UPLUS(A) of DataMatrix A is A. B is a DataMatrix object with the
%   same row names and columns names as A.  
% 
%   B = UPLUS(A) is called for the syntax '+A' when A is an object.
% 
%   See also DATAMATRIX/UMINUS.

%   Copyright 2008 The MathWorks, Inc. 


b = a;
b.Matrix = uplus(a.Matrix);
end %DataMatrix/uplus
