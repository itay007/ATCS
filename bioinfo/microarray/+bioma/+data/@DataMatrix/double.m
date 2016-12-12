function m = double(obj, rowindices, colindices)
% DOUBLE Convert DataMatrix object to a double precision array.
% 
%   B = DOUBLE(A) converts the DataMatrix object A to a double precision
%   array B.
%
%   B = DOUBLE(A, ROWS) converts specified rows of the DataMatrix object A
%   specified by ROWS to a double precision array B. ROWS can be a positive
%   integer, a vector of positive integers, a string specifying a row name,
%   a cell array containing one or more row names, or a logical vector.
%
%   B = DOUBLE(A, ROWS, COLS) converts a subset of the DataMatrix object A
%   specified by ROWS and COLS to a double precision array B. COLS can be a
%   positive integer, a vector of positive integers, a string specifying a
%   column name, cell array containing one or more column names, or a
%   logical vector.  
% 
%   See also DATAMATRIX/SINGLE

%   Copyright 2008 The MathWorks, Inc. 

try
    m = double(obj.Matrix);

    if nargin == 2
        rowindices = getDimensionIndices(obj, rowindices, 1);
        m = m(rowindices, :);
    elseif nargin >= 3
        rowindices = getDimensionIndices(obj, rowindices, 1);
        colindices = getDimensionIndices(obj, colindices, 2);
        m = m(rowindices, colindices);
    end
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','double', ME);
end
end % DataMatrix/double
