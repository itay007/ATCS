function tf = isempty(a)
%ISEMPTY True for empty DataMatrix.
% 
%   TF = ISEMPTY(A) returns true if A is an empty DataMatrix and false
%   otherwise. An empty matrix has no elements, that is PROD(SIZE(A))==0.
%  
%   See also DATASET/SIZE.

%   Copyright 2006 The MathWorks, Inc. 


tf = (a.NRows == 0) || (a.NCols == 0);
