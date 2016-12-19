function [b,idx] = sortcols(a,row,mode)
%SORTCOLS Sort columns of a DataMatrix object in ascending or descending order.
% 
%   B = SORTCOLS(A) returns a copy of the DataMatrix A with columns of
%   the matrix sorted in ascending order by all the rows in A. 
%
%   B = SORTCOLS(A,ROW) sorts the columns in A by the rows specified by
%   ROW. ROW is a positive integer, a vector of positive integers, a row
%   name, a cell array containing one or more row names, or a logical
%   vector.
% 
%   B = SORTCOLS(A,'COLNAME') sorts the columns in A by the column names.
% 
%   B = SORTCOLS(A, ROW, MODE) sort in the direction specified by MODE. If
%   MODE is 'ascend', the columns in A will be sorted by the specified rows
%   in ascending order (default); if MODE is 'descend', the columns in A
%   will be sorted by specified rows in descending order.
%
%   [B,I] = SORTCOLS(A, ...) also returns an index vector IDX such that
%   B = A(:,IDX).
%
%   See also DATAMATRIX/SORTROWS.

%   Copyright 2008-2012 The MathWorks, Inc. 


%== Input check
bioinfochecknargin(nargin,1,'DataMatrix:sortcols')

%== Initialize
if nargin < 2 || isempty(row)
    row = (1:a.NRows)';
end

if nargin < 3 || isempty(mode)
    descend = false;
elseif strcmpi(mode,'ascend')
    descend = false;
elseif strcmpi(mode,'descend')
    descend = true;
else
    error(message('bioinfo:DataMatrix:sortcols:UnrecognizedMode'));
end

%== Sort on each indexed rows, last to first
if strcmpi(row, 'colname')
    if isempty(a.ColNames)
        error(message('bioinfo:DataMatrix:sortcols:EmptyColNames'));
    else
        [~, idx] = sort(a.ColNames);
    end
else
    row = getDimensionIndices(a, row, 1);
    if strcmp(':', row)
        row = 1:a.NRows;
    end
    idx = (1:a.NCols);
    for i = flipud(row)
        row_i = a.Matrix(i,idx);
        [~, ord] = sort(row_i,2);
        idx = idx(ord);
    end
end
if descend
    idx = fliplr(idx);
end

b = a;
b.Matrix = a.Matrix(:, idx);
if ~isempty(a.ColNames)
    b.ColNames = a.ColNames(idx);
end
end % DataMatrix/sortcols
