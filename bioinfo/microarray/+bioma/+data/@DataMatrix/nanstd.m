function y = nanstd(varargin)
%NANSTD  Overload standard deviation, ignoring NaNs, for DataMatrix object.
% 
%   S = NANSTD(X) returns the standard deviation of the DataMatrix object
%   A. S is a row vector containing the standard deviation from each column
%   in A, treating NaNs as missing values. NANSTD normalizes S by (N-1) if
%   N>1, where N is the sample size. For N=1, S is normalized by N.
% 
%   S = NANSTD(A, FLAG) specifies how to normalize the data. If FLAG = 0,
%   normalizes S by N-1 if N>1. For N=1, S is normalized by N. If FLAG =1,
%   normalizes S by N, where N is the sample size. NANSTD(A,0) is the same
%   as NANSTD(A).
%
%   S = NANSTD(A,FLAG,DIM) takes the standard deviation along the dimension
%   DIM of A, treating NaNs as missing values. Default dimension is 1. If
%   DIM = 2, S is a column vector containing the standard deviation from
%   each row.
%
%   See also DATAMATRIX/STD, DATAMATRIX/NANVAR, DATAMATRIX/NANMEAN,
%   DATAMATRIX/NANMEDIAN, DATAMATRIX/NANMIN, DATAMATRIX/NANMAX,
%   DATAMATRIX/NANSUM.

%   Copyright 2008-2012 The MathWorks, Inc.


% Call nanvar(x,flag,dim) with as many inputs as needed
try 
    y = sqrt(nanvar(varargin{:}));
catch ME
     bioinfoprivate.bioclsrethrow('DataMatrix','nanstd', ME)
end
end %DataMatrix/nanstd
