function [varargout] = nanmax(varargin)
%NANMAX  Overloaded for maximum value, ignoring NaNs, for DataMatrix object.
% 
%   M = NANMAX(A) returns the largest element in DataMatrix object A with
%   NaNs treated as missing. M is a row vector containing the maximum
%   element values from each column in A.
% 
%   [M,I] = NANMAX(A) returns the indices of the maximum values in vector
%   I. If the values along the first non-singleton dimension contain more
%   than one maximal element, the index of the first one is returned.
%
%   [M,I,MN] = NANMAX(A) returns the row names of the rows with maximal
%   elements.
% 
%   C = NANMAX(A,B) returns a MATLAB array C the same size as A and B with
%   the largest elements taken from A or B. A or B must be a DataMatrix
%   object. A and B must have the same size or one can be a scalar.
%   
%   [M,I,MN] = NANMAX(A,[],DIM) returns the maximum values of A along the
%   dimension DIM. Default dimension is 1. If DIM = 2, M is a column vector
%   containing the maximum values from each row, and MN contains the column
%   names corresponding to the columns with maximal elements.
%
%   See also DATAMATRIX/MAX, DATAMATRIX/NANMEAN, DATAMATRIX/NANMEDIAN,
%   DATAMATRIX/NANMIN, DATAMATRIX/NANVAR, DATAMATRIX/NANSTD. 

%   Copyright 2008-2012 The MathWorks, Inc. 


% Call [m,ndx] = max(a,b) with as many inputs and outputs as needed
try
    [varargout{1:nargout}] = max(varargin{:});
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','nanmax', ME)
end
end %DataMatrix/nanmax
