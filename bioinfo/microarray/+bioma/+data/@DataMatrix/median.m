function out = median(obj,dim, nanflag) 
%MEDIAN Return the median values of a DataMatrix object.
% 
%   M = MEDIAN(A) returns the median value of the elements in DataMatrix
%   object A, treating NaNs as missing values. M is a row vector containing
%   the median value from each column in A.
%   
%   M = MEDIAN(A,DIM) returns the median value along the dimension DIM of
%   A. Default dimension is 1. If DIM = 2, M is a column vector containing
%   the median value from each row in A.
% 
%   M = MEDIAN(A,DIM,IGNORENAN) indicates whether NaNs should be ignored.
%   IGNORENAN can be true (default) or false. 
% 
%   See also DATAMATRIX/MAX, DATAMATRIX/MEAN, DATAMATRIX/MIN,
%   DATAMATRIX/SUM

%   Copyright 2008-2012 The MathWorks, Inc. 


bioinfochecknargin(nargin,1,'DataMatrix:median')
%== input checking
if nargin < 2 || isempty(dim)
    dim = 1;
end

if nargin < 3 || isempty(nanflag)
    nanflag = true;
end

try
    validateattributes(dim, {'numeric'}, {'integer', 'scalar', 'positive'})
catch ME 
    error(message('bioinfo:DataMatrix:median:dimensionMustBePositiveInteger'));
end
nanflag =  bioinfoprivate.opttf(nanflag, 'IGNOREFLAG','DataMatrix:median');
%== Perform operation
try
    if nanflag
        out = nanmedian(obj.Matrix, dim);
    else
        out = median(obj.Matrix, dim);
    end
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','median', ME);
end    
end %DataMatrix/median
