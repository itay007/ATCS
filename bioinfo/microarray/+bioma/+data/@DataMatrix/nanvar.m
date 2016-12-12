function out = nanvar(obj,w,dim) 
%NANVAR  Overload variance, ignoring NaNs, of DataMatrix object.
% 
%   V = NANVAR(A) returns the variance of the values in the DataMatrix
%   object A, treating NaNs as missing values. V is a row vector containing
%   the variance of each column in A. NANVAR normalizes V by N-1 if N>1,
%   where N is the sample size.  For N=1, V is normalized by N.
%
%   V = NANVAR(A,FLAG) specifies how to normalize the data. If FLAG = 0,
%   normalizes V by N-1 if N>1. For N=1, V is normalized by N. If FLAG =1,
%   normalizes V by N, where N is the sample size. NANVAR(A,0) is the same
%   as NANVAR(A).
%
%   V = NANVAR(A,W) computes the variance using the weight vector W. The
%   length of W must equal the length of the dimension over which NANVAR
%   operates, and its elements must be nonnegative. NANVAR normalizes W to
%   sum to one.
%
%   V = NANVAR(A,W,DIM) takes the variance along the dimension DIM of A.
%   Default dimension is 1. If DIM = 2, S is a column vector containing the
%   variance from each row. Pass in 0 for W to use the default
%   normalization by N-1, or 1 to use N. 
% 
%   See also DATAMATRIX/VAR, DATAMATRIX/NANSTD, DATAMATRIX/NANMEAN,
%   DATAMATRIX/NANMEDIAN, DATAMATRIX/NANMIN, DATAMATRIX/NANMAX,
%   DATAMATRIX/NANSUM.

%   Copyright 2008-2012 The MathWorks, Inc. 


bioinfochecknargin(nargin,1,'DataMatrix:nanvar')

if nargin < 2 || isempty(w)
    w = 0;
end

if nargin < 3 || isempty(dim)
    dim = 1;
end

try
    validateattributes(dim, {'numeric'}, {'integer', 'scalar', 'positive'})
catch ME 
    error(message('bioinfo:DataMatrix:nanvar:dimensionMustBePositiveInteger'));
end

try 
    out = nanvar(obj.Matrix, w, dim);
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','nanvar', ME)
end

end %DataMatrix/var
