function out = mean(obj,dim, nanflag) 
%MEAN Return the mean values of a DataMatrix object.
% 
%   M = MEAN(A) returns the mean value of the elements in DataMatrix object
%   A, treating NaNs as missing values. M is a row vector containing the
%   mean value of each column in A. 
%   
%   M = MEAN(A,DIM) returns the mean value along the dimension DIM of A.
%   Default dimension is 1. If DIM = 2, M is a column vector containing the
%   mean value of each row in A.
% 
%   M = MEAN(A,DIM,IGNORENAN) indicates whether NaNs should be ignored.
%   IGNORENAN can be true (default) or false. 
% 
%   See also DATAMATRIX/MAX, DATAMATRIX/MEDIAN, DATAMATRIX/MIN,
%   DATAMATRIX/SUM

%   Copyright 2008-2012 The MathWorks, Inc. 


bioinfochecknargin(nargin,1,'DataMatrix:mean')
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
    error(message('bioinfo:DataMatrix:mean:dimensionMustBePositiveInteger'));
end
nanflag =  bioinfoprivate.opttf(nanflag, 'IGNOREFLAG', 'DataMatrix:mean');

%== Perform operation
try
    if nanflag
        out = nanmean(obj.Matrix, dim);
    else
        out = mean(obj.Matrix, dim);
    end
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','mean', ME);
end
end %DataMatrix/mean
