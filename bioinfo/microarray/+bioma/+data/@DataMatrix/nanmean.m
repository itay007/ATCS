function m = nanmean(x,dim)
%NANMEAN  Overloaded mean value, ignoring NaNs, for DataMatrix.
% 
%   M = NANMEAN(A) returns the mean value of the elements in DataMatrix
%   object A, treating NaNs as missing values. M is a row vector containing
%   the mean values of each column in A. 
%
%   NANMEAN(A,DIM) returns the mean value along the dimension DIM of A,
%   treating NaNs as missing values. Default dimension is 1. If DIM = 2, M
%   is a column vector containing the mean value of each row.
%
%   See also DataMatrix/MEAN, DataMatrix/NANMEDIAN, DataMatrix/NANSTD,
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
        error(message('bioinfo:DataMatrix:nanmean:dimensionMustBePositiveInteger'));
    end

    m = nanmean(x.Matrix, dim);
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','nanmean', ME); 
end
end %DataMatrix/nanmean
