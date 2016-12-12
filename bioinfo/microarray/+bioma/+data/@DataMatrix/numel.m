function n = numel(a,varargin)
%NUMEL Number of elements in a DataMatrix object.
% 
%   N = NUMEL(A) returns 1. To find the number of elements, N, in the
%   DataMatrix object A, use PROD(SIZE(A)) or NUMEL(A,':',':').
%
%   N = NUMEL(A, VARARGIN) returns the number of subscripted elements, N,
%   in A(index1, index2), where VARARGIN is a cell array whose elements are
%   index1 and index2.
%
%   See also DATAMATRIX/SIZE, DATAMATRIX/LENGTH.

%   Copyright 2008 The MathWorks, Inc. 


switch nargin
    case 1
        %== Return 1 for subsref NARGOUT 
        n = 1;
    otherwise
        try
            n = numel(a.Matrix, varargin{:});
        catch ME
            bioinfoprivate.bioclsrethrow('DataMatrix','numel', ME);
        end
end

end % numel
