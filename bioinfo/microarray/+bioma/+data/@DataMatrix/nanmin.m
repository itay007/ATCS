function [varargout] = nanmin(varargin)
%NANMIN Overloaded minimum value, ignoring NaNs, for DataMatrix object.
% 
%   M = NANMIN(A) returns the smallest element in DataMatrix object A with
%   NaNs treated as missing. M is a row vector containing the minimum
%   element from each column.
% 
%   [M,I] = NANMIN(A) returns the indices of the minimum values in vector
%   I. If the values along the first non-singleton dimension contain more
%   than one minimal element, the index of the first one is returned.
%
%   [M,I,MN] = NANMIN(A)  returns the row names of the rows with minimal
%   elements.
% 
%   C = NANMIN(A,B)  returns a MATLAB array C the same size as A or B with
%   smallest elements taken from A or B. A or B must be a DataMatrix
%   object. A and B must have the same size or one can be a scalar.
% 
%   [M,I,MN] = NANMIN(A,[], DIM) returns the minimum values of A along the
%   dimension DIM. Default dimension is 1. If DIM = 2, M is a column vector
%   containing the minimum values from each row, MN contains the column
%   names of columns with minimal elements.
%
%   See also DATAMATRIX/MIN, DATAMATRIX/NANMAX, DATAMATRIX/NANMEAN,
%   DATAMATRIX/NANMEDIAN, DATAMATRIX/NANVAR, DATAMATRIX/NANSTD. 

%   Copyright 2008-2012 The MathWorks, Inc. 


% Call [m,ndx] = min(a,b) with as many inputs and outputs as needed
try
    [varargout{1:nargout}] = min(varargin{:});
catch ME
     bioinfoprivate.bioclsrethrow('DataMatrix','nanmin', ME)
end
end %DataMatrix/nanmin
