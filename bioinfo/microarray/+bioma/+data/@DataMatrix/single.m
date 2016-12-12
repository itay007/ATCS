function m = single(obj,rowindices, colindices)
% SINGLE Convert DataMatrix object to a single precision array.
% 
%   B = SINGLE(A) converts the DataMatrix object A to a single precision
%   array B.
%   
%   B = SINGLE(A, ROWS) converts specified rows of the DataMatrix object A
%   specified by ROWS to a single precision array B. ROWS can be a positive
%   integer, a vector of positive integers, a string specifying a row name,
%   a cell array containing one or more row names, or a logical vector.
%
%   M = SINGLE(A, ROWS, COLS) converts a subset of the DataMatrix object A
%   specified by ROWS and COLS to a single precision array B. COLS can be a
%   positive integer, a vector of positive integers, a string specifying a
%   column name, a cell array containing one or more column names, or a
%   logical vector.  
% 
%   See also DATAMATRIX/DOUBLE.

%   Copyright 2008-2012 The MathWorks, Inc. 


%== Input check
try
    m = single(obj.Matrix);

    if nargin == 2
        rowindices = getDimensionIndices(obj, rowindices, 1);
        m = m(rowindices, :);
    elseif nargin >= 3
        rowindices = getDimensionIndices(obj, rowindices, 1);
        colindices = getDimensionIndices(obj, colindices, 2);
        m = m(rowindices, colindices);
    end
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','single', ME);
end
end % DataMatrix/single
