function out = var(obj,w,dim,nanflag) 
%VAR Return the variances of a DataMatrix object.
% 
%   V = VAR(A) returns the variance of the values in the DataMatrix object
%   A,  treating NaNs as missing values. V is a row vector containing the
%   variance of each column in A. VAR normalizes V by N-1 if N>1, where N
%   is the sample size.  For N=1, V is normalized by N.
% 
%   V = VAR(A,FLAG) specifies how to normalize the data. If FLAG = 0,
%   normalizes V by N-1 if N>1. For N=1, V is normalized by N. If FLAG =1,
%   normalizes V by N, where N is the sample size. VAR(A,0) is the same as
%   VAR(A).
%
%   V = VAR(A,W) computes the variance using the weight vector W. The
%   length of W must equal the length of the dimension over which VAR
%   operates, and its elements must be nonnegative. VAR normalizes W to
%   sum to one.
%
%   V = VAR(A,W,DIM) takes the variance along the dimension DIM of A.
%   Default dimension is 1. If DIM = 2, V is a column vector containing the
%   variance for each row. Pass in 0 for W to use the default normalization
%   by N-1, or 1 to use N. 
% 
%   V = VAR(A,W,DIM,IGNORENAN) indicates whether NaNs should be ignored.
%   IGNORENAN can be true (default) or false. 
% 
%   The variance is the square of the standard deviation (STD).
% 
%   See also DATAMATRIX/IQR, DATAMATRIX/MEAN, DATAMATRIX/MEDIAN,
%            DATAMATRIX/STD. 

%   Copyright 2008-2012 The MathWorks, Inc. 


bioinfochecknargin(nargin,1,'DataMatrix:var')

if nargin < 2 || isempty(w)
    w = 0;
end

if nargin < 3 || isempty(dim)
    dim = 1;
end

if nargin < 4 || isempty(nanflag)
    nanflag = true;
end

try
    validateattributes(dim, {'numeric'}, {'integer', 'scalar', 'positive'})
catch ME 
    error(message('bioinfo:DataMatrix:var:dimensionMustBePositiveInteger'));
end
nanflag =  bioinfoprivate.opttf(nanflag, 'IGNOREFLAG', 'DataMatrix:var');

try
    if nanflag
        out = nanvar(obj.Matrix, w, dim);
    else
        out = var(obj.Matrix, w, dim);
    end
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','var', ME)
end

end %DataMatrix/var
