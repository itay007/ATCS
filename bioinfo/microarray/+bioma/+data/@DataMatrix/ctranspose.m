function b = ctranspose(a)
%CTRANSPOSE Overload transpose for DataMatrix object.
% 
%   B = CTRANSPOSE(A) returns the transpose of the 2-dimensional DataMatrix
%   object A.  Note that CTRANSPOSE is identical to TRANSPOSE for
%   DataMatrix objects.
%
%   CTRANSPOSE is called for the syntax A'.
%
%   See also DATAMATRIX/TRANSPOSE.

%   Copyright 2008 The MathWorks, Inc. 


b = bioma.data.DataMatrix(a.Matrix', a.ColNames, a.RowNames);
