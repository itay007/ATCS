function [c, varargout] = max(a, b, dim) 
%MAX Return the maximum values of a DataMatrix object.
% 
%   M = MAX(A) returns the largest element in DataMatrix object A. M is a
%   row vector containing the maximum element value from each column in A.
% 
%   [M,I] = MAX(A) returns the indices of the maximum values to vector I.
%   If the values along the first non-singleton dimension contain more
%   than one maximal element, the index of the first one is returned.
%
%   [M,I,MN] = MAX(A) returns the row names of the rows with maximal
%   elements.
% 
%   C = MAX(A,B) returns a MATLAB array C the same size as A and B with the
%   largest elements taken from A or B. A or B must be a DataMatrix object.
%   A and B must have the same size, or one can be a scalar.
%   
%   [M,I,MN] = MAX(A,[],DIM) returns the maximum values of A along the
%   dimension DIM. Default dimension is 1. If DIM = 2, M is a column vector
%   containing the maximum values from each row, and MN contains the column
%   names corresponding to the columns with maximal elements.
% 
%   See also DATAMATRIX/MEAN, DATAMATRIX/MEDIAN, DATAMATRIX/MIN,
%   DATAMATRIX/SUM

%   Copyright 2008-2012 The MathWorks, Inc. 


bioinfochecknargin(nargin,1,'DataMatrix:max')

if nargin < 3
    dim = 1;
end

if nargin < 2 || isequal(b,[])
    b = [];
elseif isa(a, 'bioma.data.DataMatrix') && isa(b, 'bioma.data.DataMatrix')
    [c, varargout{1:nargout-1}] = doMax(a.Matrix, b.Matrix, nargin);
elseif isa(a, 'bioma.data.DataMatrix') && (isnumeric(b) || islogical(b))
    [c, varargout{1:nargout-1}] = doMax(a.Matrix, b, nargin);
elseif isa(b, 'bioma.data.DataMatrix') && (isnumeric(a) || islogical(a))
    [c, varargout{1:nargout-1}] = doMax(a, b.Matrix, nargin);
end

%== Only one matrix
try
    if isempty(b)
        [c, varargout{1}] = max(a.Matrix, [], dim);
        
        if nargout > 2
            if dim == 1
                varargout{2} = a.RowNames(varargout{1})';
            elseif dim ==2
                varargout{2} = a.ColNames(varargout{1})';
            end
        end
    end
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','max', ME)
end
end %DataMatrix/max

function varargout = doMax(p1, p2, ninarg)
% If ninarg > 2, it should error
try
    if ninarg > 2
        varargout{:} = max(p1, p2, 1);
    else
        varargout{:} = max(p1, p2);
    end
catch ME
   bioinfoprivate.bioclsrethrow('DataMatrix','max', ME)  
end
end % doMax
