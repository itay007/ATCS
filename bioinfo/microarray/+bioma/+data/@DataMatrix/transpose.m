function b = transpose(a)
%TRANSPOSE Overload transpose for DataMatrix object.
% 
%   B = TRANSPOSE(A) returns the transpose of the 2-dimensional DataMatrix
%   object A.  Note that CTRANSPOSE is identical to TRANSPOSE for
%   DataMatrix objects.
%
%   TRANSPOSE is called for the syntax A.'.
%
%   See also DATAMATRIX/CTRANSPOSE.

%   Copyright 2008 The MathWorks, Inc. 


b = bioma.data.DataMatrix(a.Matrix.', a.ColNames, a.RowNames);
