function [c, varargout] = min(a, b, dim) 
%MIN Return the minimum values of a DataMatrix object.
% 
%   M = MIN(A) returns the smallest element in DataMatrix object A. M is a
%   row vector containing the minimum element value from each column in A.
% 
%   [M,I] = MIN(A) returns the indices of the minimum values to vector I.
%   If the values along the first non-singleton dimension contain more
%   than one minimal element, the index of the first one is returned.
%
%   [M,I,MN] = MIN(A) returns the row names of the rows with minimal
%   elements.
% 
%   C = MIN(A,B) returns a MATLAB array C the same size as A or B with the
%   smallest elements taken from A or B. A or B must be a DataMatrix
%   object. A and B must have the same size, or one can be a scalar.
%   
%   [M,I,MN] = MIN(A,[], DIM) returns the minimum values of A along the
%   dimension DIM. Default dimension is 1. If DIM = 2, M is a column vector
%   containing the minimum values from each row, and MN contains the column
%   names corresponding to the columns with minimal elements. 
% 
%   See also DATAMATRIX/MAX, DATAMATRIX/MEAN, DATAMATRIX/MEDIAN,
%   DATAMATRIX/SUM

%   Copyright 2008 The MathWorks, Inc. 


bioinfochecknargin(nargin,1,'DataMatrix:min')

if nargin < 3
    dim = 1;
end

if nargin < 2 || isequal(b,[])
    b = [];
elseif isa(a, 'bioma.data.DataMatrix') && isa(b, 'bioma.data.DataMatrix')
    [c, varargout{1:nargout-1}] = doMin(a.Matrix, b.Matrix, nargin);
elseif isa(a, 'bioma.data.DataMatrix') && (isnumeric(b) || islogical(b))
    [c, varargout{1:nargout-1}] = doMin(a.Matrix, b, nargin);
elseif isa(b, 'bioma.data.DataMatrix') && (isnumeric(a) || islogical(a))
    [c, varargout{1:nargout-1}] = doMin(a, b.Matrix, nargin);
end

%== Only one matrix
if isempty(b)
    [c, varargout{1}] = min(a.Matrix, [], dim);

    if nargout > 2
        if dim == 1
            varargout{2} = a.RowNames(varargout{1})';
        elseif dim ==2
            varargout{2} = a.ColNames(varargout{1})';
        end        
    end
end
end %DataMatrix/min

function varargout = doMin(p1, p2, ninarg)
% If ninarg > 2, it should error
try
    if ninarg > 2
        varargout{:} = min(p1, p2, 1);
    else
        varargout{:} = min(p1, p2);
    end
catch ME
   bioinfoprivate.bioclsrethrow('DataMatrix','min', ME)  
end
end % doMin
