function e = end(a,k,n)
%END Last index in an indexing expression for a DataMatrix.
% 
%   END(A,K,N) is called for indexing expressions involving the DataMatrix
%   object A when END is part of the K-th index out of N indices.  For
%   example, the expression A(end-1,:) calls A's END method with
%   END(A,1,2).
%
%   See also DataMatrix/SIZE.

%   Copyright 2008 The MathWorks, Inc. 


switch k
case 1
    e = a.NRows;
case 2
    e = a.NCols;
otherwise
    e = 1;
end
