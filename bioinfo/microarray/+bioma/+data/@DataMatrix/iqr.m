function out = iqr(obj,dim) 
%IQR Return the interquartile range of the values of a DataMatrix object.
% 
%   Y = IQR(A) returns the interquartile range of the values in DataMatrix
%   A. Y is a row vector containing the interquartile range of each column
%   of A.
%   
%   Y = IQR(A,DIM) returns the interquartile range along the dimension DIM
%   of A. Default dimension is 1. If DIM = 2, M is a column vector
%   containing the interquartile range of each row.
% 
%   See also DATAMATRIX/PRCTILE, DATAMATRIX/STD, DATAMATRIX/VAR

%   Copyright 2008-2012 The MathWorks, Inc. 


bioinfochecknargin(nargin,1,['DataMatrix:' mfilename])

if nargin < 2 || isempty(dim)
    dim = 1;
end

try
    validateattributes(dim, {'numeric'}, {'integer', 'scalar', 'positive'})
catch ME 
    error(message('bioinfo:DataMatrix:iqr:dimensionMustBePositiveInteger'));
end

if dim > obj.NDims
    error(message('bioinfo:DataMatrix:iqr:ExceedDimension', dim, obj.NDims));
end

out = diff(prctile(obj.Matrix, [25; 75], dim), [], dim);
end %DataMatrix/iqr
