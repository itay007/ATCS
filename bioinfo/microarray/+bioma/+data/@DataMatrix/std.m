function out = std(obj,w,dim,nanflag) 
%STD Return the standard deviations of a DataMatrix object.
% 
%   S = STD(A) returns the standard deviation of the DataMatrix object A. S
%   is a row vector containing the standard deviation from each column in
%   A, treating NaNs as missing values. STD normalizes S by (N-1) if N>1,
%   where N is the sample size. For N=1, S is normalized by N.
% 
%   S = STD(A, FLAG) specifies how to normalize the data. If FLAG = 0,
%   normalizes S by N-1 if N>1. For N=1, S is normalized by N. If FLAG =1,
%   normalizes S by N, where N is the sample size. STD(A,0) is the same as
%   STD(A).
%
%   S = STD(A,FLAG,DIM) takes the standard deviation along the dimension
%   DIM of A. Default dimension is 1. If DIM = 2, S is a column vector
%   containing the standard deviation from each row in A.
% 
%   S = STD(A,FLAG,DIM,IGNORENAN) indicates whether NaNs should be ignored.
%   IGNORENAN can be true (default) or false.
% 
%   See also DATAMATRIX/IQR, DATAMATRIX/MEAN, DATAMATRIX/MEDIAN,
%   DATAMATRIX/VAR 

%   Copyright 2008-2012 The MathWorks, Inc. 


bioinfochecknargin(nargin,1,'DataMatrix:std')

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
    error(message('bioinfo:DataMatrix:std:dimensionMustBePositiveInteger'));
end
nanflag =  bioinfoprivate.opttf(nanflag, 'IGNOREFLAG', 'DataMatrix:std');

try 
    if nanflag
        out = nanstd(obj.Matrix, w, dim);
    else
        out = std(obj.Matrix, w, dim);
    end
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','std', ME)
end

end %DataMatrix/std
