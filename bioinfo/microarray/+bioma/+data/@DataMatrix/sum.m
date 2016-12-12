function out = sum(obj,dim,nanflag) 
%SUM  Sum of elements of a DataMatrix object.
% 
%   S = SUM(A) returns the sum of the elements in the DataMatrix object A,
%   treating NaNs as missing values. S is a row vector with the sum for
%   each column in A. If the values of A are double or single, S has the
%   same class, otherwise S is a double.
% 
%   S = SUM(A,DIM) sums along the dimension DIM. Default dimension is 1. If
%   DIM = 2, S is a column vector with the sum for each row in A.
%
%   S = SUM(A,DIM,IGNORENAN) indicates whether NaNs should be ignored.
%   IGNORENAN can be true (default) or false. 
% 
%   See also DATAMATRIX/MAX, DATAMATRIX/MIN

%   Copyright 2008-2012 The MathWorks, Inc. 


if nargin < 2 || isempty(dim)
    dim = 1;
end

if nargin < 3 || isempty(nanflag)
    nanflag = true;
end

try
    validateattributes(dim, {'numeric'}, {'integer', 'scalar', 'positive'})
catch ME 
    error(message('bioinfo:DataMatrix:sum:dimensionMustBePositiveInteger'));
end
nanflag =  bioinfoprivate.opttf(nanflag, 'IGNOREFLAG', 'DataMatrix:sum');

try
    if nanflag
        out = nansum(obj.Matrix, dim);
    else
        out = sum(obj.Matrix, dim);
    end
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','sum', ME)
end
end %DataMatrix/sum
