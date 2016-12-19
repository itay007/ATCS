function m = nanmedian(x,dim)
%NANMEDIAN Overload median value, ignoring NaNs, for DataMatrix object.
% 
%   M = NANMEDIAN(A) returns the median value of the elements in DataMatrix
%   object A, treating NaNs as missing values. M is a row vector containing
%   the median values from each column in A.
%
%   NANMEDIAN(A,DIM) returns the median value along the dimension DIM of A,
%   treating NaNs as missing values. Default dimension is 1. If DIM = 2, M
%   is a column vector containing the median value from each row.
%
%   See also DataMatrix/MEDIAN, DataMatrix/NANMEAN, DataMatrix/NANSTD,
%   DataMatrix/NANVAR, DataMatrix/NANMIN, DataMatrix/NANMAX,
%   DataMatrix/NANSUM.  

%   Copyright 2008-2012 The MathWorks, Inc.


% Find NaNs and set them to zero
try
    if nargin < 2 || isempty(dim)
        dim = 1;
    end
    
    try
        validateattributes(dim, {'numeric'}, {'integer', 'scalar', 'positive'})
    catch ME 
        error(message('bioinfo:DataMatrix:nanmedian:dimensionMustBePositiveInteger'));
    end

    m = nanmedian(x.Matrix, dim);
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','nanedian', ME); 
end
end %DataMatrix/nanmedian
