function a = subsasgn(a,s,b)
%SUBSASGN Subscripted assignment to a DataMatrix object.
% 
%   A = SUBSASGN(A,S,B) is called for the syntax A(I,J)=B, A.(ROW)(COL)=B
%   when A is a DataMatrix.  S is a structure array with the fields:
%       type -- string containing '()', '{}', or '.' specifying the
%               subscript type.
%       subs -- Cell array or string containing the actual subscripts.
%
%   A(I,J) = B assigns the values of the DataMatrix object or a valid
%   MATLAB numeric or logical array, B into the elements of the rectangular
%   submatrix of DataMatrix object A specified by the subscript vectors I
%   and J. I and J can be a positive integer, a vector of positive
%   integers, a string of name, a cell arrays containing one or more unique
%   names, or a logical vector.  The assignment does not use row names,
%   column names, or any other properties of B to modify properties of A;
%   Elements of B are assigned into A by position, not by matching names.
%   The assignment expanding the number of rows or columns in A is not
%   allowed. No internal reshaping of A is allowed: the size of B must
%   match the size of submatrix of A to be replaced. 
%
%   A.(I)(J) = B assigns the values of B into the elements of the
%   rectangular submatrix of DataMatrix object A specified by the subscript
%   vectors I and J. The subscripts I and J provide subscripting into a
%   DataMatrix matrix. I and J can be a positive integer, a vector of
%   positive integers, a string of name, a cell arrays containing one or
%   more unique names, or a logical vector.. A.(I) are equivalent to
%   A.(I)(':'). ROWNAME is a cell array of column names. B can be a
%   DataMatrix object or a valid MATLAB numeric or logical array.
%
%   Note: Cell indexing A{I,J}= B is not supported.
%
%   A.PROPERTYNAME = P assigns to a DataMatrix property. PROPERTYNAME is
%   'NAME', 'ROWNAMES', or 'COLNAMES'.  To assign into an element of the
%   property, A.PROPERTYNAME may also be followed by further subscripting
%   as supported by the property.
%
%   See also DATAMATRIX, DATAMATRIX/SUBSREF, SET.

%   Copyright 2008-2012 The MathWorks, Inc.


%== Empty DataMatrix
creating = isequal(a,[]);
if creating
    a = bioma.data.DataMatrix;
end

switch s(1).type
    case '()'
        % '()' is assignment into a subset of rows/columns from another DataMatrix
        if numel(s(1).subs) ~= a.ndims
            error(message('bioinfo:DataMatrix:subsasgn:NDSubscript'));
        elseif ~isscalar(s)
            error(message('bioinfo:DataMatrix:subsasgn:InvalidSubscriptExpression'));
        end
        
        % Translate row names into indices (leaves ':' alone)
        checkCellNameUnique(s(1).subs{1}, 1);
        [rowIndices, numRowIndices] = getDimensionIndices(a,s(1).subs{1}, 1, 0);
        if any(rowIndices == 0) && ~islogical(rowIndices)
            rowSubNames = s(1).subs{1};
            error(message('bioinfo:DataMatrix:subsasgn:UnrecognizedRowName', commaSeparatedList( rowSubNames( rowIndices==0 ) )))  
        elseif isempty(rowIndices) && ~isempty(a.RowNames)
            error(message('bioinfo:DataMatrix:subsasgn:UnrecognizedRowName', commaSeparatedList( s( 1 ).subs{ 1 } )))              
        end
        % Translate column names into indices (translates ':')
        checkCellNameUnique(s(1).subs{2}, 2);
        [colIndices, numColIndices] =  getDimensionIndices(a, s(1).subs{2}, 2, 0);
        if any(colIndices == 0) && ~islogical(colIndices)
            colSubNames = s(1).subs{2};
            error(message('bioinfo:DataMatrix:subsasgn:UnrecognizedColumnName', commaSeparatedList( colSubNames( colIndices==0 ) )))  
        elseif isempty(colIndices) && ~isempty(a.ColNames)
            error(message('bioinfo:DataMatrix:subsasgn:UnrecognizedColumnName', commaSeparatedList( s( 1 ).subs{ 2 } )))  
        end
        %== Empty assignment
        %Syntax: a(rows,:) = [];
        %        a(:, cols) = [];
        %        a(rows, cols) = [] is not allowed
        if issqrbrkt(b)
            try
               a.Matrix(rowIndices, colIndices) = b;
            catch ME
               bioinfoprivate.bioclsrethrow('DataMatrix','subsasgn', ME)
            end   
            
            if bioma.util.isColon(s(1).subs{2})
                a.NRows = size(a.Matrix, 1);
                if ~isempty(a.RowNames)
                    a.RowNames(rowIndices) = [];
                end
            elseif bioma.util.isColon(s(1).subs{1})
                a.NCols = size(a.Matrix, 2);
                if ~isempty(a.ColNames)
                    a.ColNames(colIndices) = [];
                end
            else
                error(message('bioinfo:DataMatrix:subsasgn:InvalidEmptyAssignment'));
            end
            
            %== Assignment from a DataMatrix
            %Syntax: a(rows, cols) = b
            %
            %This operation is supposed to replace the values in DataMatrix a
            %by values in DataMatrix b. No internal reshaping is allowed: the
            %size of b must match the size of submatrix of A to be replaced.
        elseif isa(b, 'bioma.data.DataMatrix')
            if ~isscalar(b) %allow scalar expansion
                if b.NRows ~= numRowIndices
                    error(message('bioinfo:DataMatrix:subsasgn:RowDimensionMismatch'));
                end
                
                if b.NCols ~= numColIndices
                    error(message('bioinfo:DataMatrix:subsasgn:ColumnDimensionMismatch'));
                end
            end
            
            s.subs{1} = rowIndices;
            s.subs{2} = colIndices;
            if isempty(a)
                a = bioma.data.DataMatrix(subsasgn(a.Matrix, s, b.Matrix), b.RowNames, b.ColNames);
            else
                a.Matrix = subsasgn(a.Matrix, s, b.Matrix);
            end
            %== Assignment from a numeric or a logical matrix
            %Syntax: a(rows, cols) = b
            %
            %This operation is supposed to replace the values in DataMatrix a
            %by values in matrix b. No internal reshaping is allowed: the
            %size of b must match the size of submatrix of A to be replaced.
        elseif isnumeric(b) || islogical(b)
            if ~isscalar(b) %allow scalar expansion
                [b_nRows, b_nCols] = size(b);
                
                if b_nRows ~= numRowIndices
                    error(message('bioinfo:DataMatrix:subsasgn:RowDimensionMismatch'));
                end
                
                if b_nCols ~= numColIndices
                    error(message('bioinfo:DataMatrix:subsasgn:ColumnDimensionMismatch'));
                end
            end
            
            if isempty(a)
                a = bioma.data.DataMatrix(subsasgn(a.Matrix, s, b));
            else
                s.subs{1} = rowIndices;
                s.subs{2} = colIndices;
                a.Matrix = subsasgn(a.Matrix, s, b);
            end
            
        else
            error(message('bioinfo:DataMatrix:subsasgn:InvalidAssignment'));
        end
    case '{}'
        error(message('bioinfo:DataMatrix:subsasgn:CellSubscript'));
    case '.'
        % Using '.' to reference to elements in a DataMatrix
        rowSubName = s(1).subs;
        %== Get row indices first
        checkCellNameUnique(rowSubName, 1);
        rowIndices =  getDimensionIndices(a, rowSubName, 1, 0);
        
        if ~isempty(rowIndices)
            % Using cell strings, if there are unrecognized names, it should error:
            if any(rowIndices == 0) && ~islogical(rowIndices)
                ukNames = rowSubName(rowIndices == 0);
                error(message('bioinfo:DataMatrix:subsasgn:UnrecognizedRowName', commaSeparatedList( ukNames )))  
            end
            
            %== Number of subscripts must be less or equal to NDims
            if numel(s) > a.NDims
                error(message('bioinfo:DataMatrix:subsasgn:DotNDimSubscript'));
            end
            
            if numel(s) == 1
                colSubName = ':';
            elseif numel(s) == 2
                if strcmp(s(2).subs{:}, ':')
                    colSubName = ':';
                else
                    colSubName = s(2).subs{:};
                end
            end
            checkCellNameUnique(colSubName, 2);
            colIndices =  getDimensionIndices(a, colSubName, 2);
            if any(colIndices == 0) && ~islogical(colIndices)
                colSubNames = s(1).subs{1};
                error(message('bioinfo:DataMatrix:subsasgn:UnrecognizedColumnName', commaSeparatedList( colSubNames( colIndices==0 ) )))
            elseif isempty(colIndices) && ~isempty(a.ColNames)
                error(message('bioinfo:DataMatrix:subsasgn:UnrecognizedColumnName', commaSeparatedList( colSubName )))
            end
            s(1).type = '()';
            s(1).subs = {rowIndices, colIndices};
            a = subsasgn(a, s(1), b);
        else
            % Try to figure out the best error to throw
            name = s(1).subs;
            if ~(ischar(name) && isvector(name) && (size(name,1)==1))
                error(message('bioinfo:DataMatrix:subsasgn:InvalidPropertyName'));
            end
            bioinfoprivate.pvpair(name, [],[fieldnames(obj); {'Matrix'}], 'DataMatrix:subsasgn');
            error(message('bioinfo:DataMatrix:subsasgn:InvalidPropertyAssignment', rowSubName));
        end
end %switch
end %DataMatrix/subsasgn

function tf = issqrbrkt(a)
%ISSQRBRKT Check if an array is the square bracket '[]'

% This checking does not work for the ones(0,0) like cases
tf = isnumeric(a)&&isempty(a)&&~isobject(a)&&~isstruct(a);
end

function list = commaSeparatedList(list)
% Returns a comma separated list as a single string
if iscellstr(list)
   list = [sprintf('%s, ',list{1:end-1}) list{end}];
end
end

function checkCellNameUnique(nameCell, dim)
% Check if the names in the cell array nameCell is unique. Error is if not
% unique.
if iscellstr(nameCell) && (numel(unique(nameCell))~=numel(nameCell))
    dimStr = {'Row', 'Column'};
    switch dimStr{dim}
        case 'Row'
            error(message('bioinfo:DataMatrix:subsasgn:NonUniqueRowNames'));
        case 'Column'
            error(message('bioinfo:DataMatrix:subsasgn:NonUniqueColumnNames'))
    end
end
end
