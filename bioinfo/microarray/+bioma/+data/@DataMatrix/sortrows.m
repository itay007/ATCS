function [b,idx] = sortrows(a,col,mode)
%SORTROWS Sort rows of a DataMatrix object in ascending or descending order.
% 
%   B = SORTROWS(A) returns a copy of the DataMatrix A with rows of the
%   matrix sorted in ascending order by all the columns in A. 
%
%   B = SORTROWS(A,COL) sorts the rows in A by the columns specified by
%   COL. COL is a positive integer, a vector of positive integers, a column
%   name, a cell array containing one or more column names, or a logical
%   vector.
% 
%   B = SORTROWS(A,'ROWNAME') sorts the rows in A by the row names.
% 
%   B = SORTROWS(A, COL, MODE) sort in the direction specified by MODE. If
%   MODE is 'ascend' the rows in A will be sorted by the specified columns
%   in ascending order (default); if MODE is 'descend', the rows in A will
%   be sorted by specified columns in descending order.
%
%   [B,I] = SORTROWS(A, ...) also returns an index vector IDX such that B =
%   A(IDX,:).
%
%   See also DATAMATRIX/SORTCOLS.

%   Copyright 2008-2012 The MathWorks, Inc. 


%== Input check
bioinfochecknargin(nargin,1,'DataMatrix:sortrows')

%== Initialize
if nargin < 2 || isempty(col)
    col = 1:a.NCols;
end

if nargin < 3 || isempty(mode)
    descend = false;
elseif strcmpi(mode,'ascend')
    descend = false;
elseif strcmpi(mode,'descend')
    descend = true;
else
    error(message('bioinfo:DataMatrix:sortrows:UnrecognizedMode'));
end

%== Sort on each index columns, last to first
if strcmpi(col, 'rowname') 
    if isempty(a.RowNames)
        error(message('bioinfo:DataMatrix:sortrows:EmptyRowNames'));
    else
        [~, idx] = sort(a.RowNames);
    end
else
    col = getDimensionIndices(a, col, 2);
    if strcmp(':', col)
        col = 1:a.NCols;
    end
    idx = (1:a.NRows)';
    for i = fliplr(col)
        col_i = a.Matrix(idx, i);
        [~, ord] = sort(col_i);
        idx = idx(ord);
    end
end
if descend
    idx = flipud(idx);
end

b = a;
b.Matrix = a.Matrix(idx, :);
if ~isempty(a.RowNames)
    b.RowNames = a.RowNames(idx);
end
end % DataMatrix/sortrow
    
