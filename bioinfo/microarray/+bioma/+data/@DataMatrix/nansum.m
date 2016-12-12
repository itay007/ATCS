function out = nansum(obj,dim) 
%NANSUM  Overload sum, ignoring NaNs, for DataMatrix object.
% 
%   S = NANSUM(A) returns the sum of the elements of the DataMatrix A,
%   treating NaNs as missing values. S is a row vector with the sum over
%   each column. If the values of A is double or single, S has the same
%   class, otherwise S is double.
%
%   S = NANSUM(A,DIM) sums along the dimension DIM. Default dimension is 1.
%   If DIM = 2, S is a column vector with the sum over each row.
%
%   See also DATAMATRIX/SUM, DATAMATRIX/NANMEAN, DATAMATRIX/NANVAR,
%   DATAMATRIX/NANSTD, DATAMATRIX/NANMIN, DATAMATRIX/NANMAX,
%   DATAMATRIX/NANMEDIAN.  

%   Copyright 2008-2012 The MathWorks, Inc. 


if nargin < 2 || isempty(dim)
    dim = 1;
end

try
    validateattributes(dim, {'numeric'}, {'integer', 'scalar', 'positive'})
catch ME 
    error(message('bioinfo:DataMatrix:nansum:dimensionMustBePositiveInteger'));
end
try
    out = nansum(obj.Matrix, dim);
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','nansum', ME)
end
end %DataMatrix/nansum
