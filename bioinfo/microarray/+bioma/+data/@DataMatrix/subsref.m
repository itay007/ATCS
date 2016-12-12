function [varargout] = subsref(a,s)
%SUBSREF Subscripted reference for a DataMatrix object.
%
%   B = SUBSREF(A,S) is called for the syntax A(I,J), or A.(Row)(Col)
%   when A is a DataMatrix object.  S is a structure array with the fields:
%       type -- string containing '()', or '.' specifying the
%               subscript type.
%       subs -- Cell array or string containing the actual subscripts.
%
%   B = A(I,J) returns a DataMatrix object that contains a subset of the
%   rows and columns in the DataMatrix object A. I and J can be a positive
%   integer, a vector of positive integers, a string of names, a cell array
%   containing one or more unique names, or a logical vector. B contains
%   the same property values as A, subsetted for rows or columns where
%   appropriate. () subscribe type does not allow linear indexing.
%
%   B = A.(I)(J) returns a subset of matrix elements, B, of a DataMatrix
%   object A. The subscripts I and J provide subscripting into a DataMatrix
%   matrix. I and J can be a positive integer, a vector of positive
%   integers, a string of names, a cell array containing one or more unique
%   names, or a logical vector. A.ROWNAME or A.(I) are equivalent to
%   A.(I)(':'). ROWNAME is a row name literal or a cell array of column
%   names.
%
%   Note: Cell indexing A{I,J} is not supported.
%
%   P = A.PROPERTYNAME returns a DataMatrix property. Further subscripting
%   to properties 'RowNames', 'ColNames' is also followed.
%
%   See also DATAMATRIX, DATAMATRIX/SUBSASGN, DATAMATRIX/GET.

%   Copyright 2008-2012 The MathWorks, Inc.


switch s(1).type
    case '()'  % '()' is a reference to a subset of a DataMatrix object.
        
        % Parenthesis indexing can only return a single thing.
        if nargout > 1
            error(message('bioinfo:DataMatrix:subsref:TooManyOutputs'));
            % No cascaded subscripts are allowed to follow parenthesis indexing.
        elseif ~isscalar(s)
            error(message('bioinfo:DataMatrix:subsref:InvalidSubscriptExpr'));
        elseif numel(s(1).subs) ~= a.NDims
            error(message('bioinfo:DataMatrix:subsref:NDSubscript'));
        end
        
        %== Return if rowIndices is empty -[] case
        if issqrbrkt(s(1).subs{1})
            b = subsref(a.Matrix,s);
            varargout{1} = bioma.data.DataMatrix(b);
            return;
        end
        
        try
            % Translate row names into indices (leaves ':' alone)
            checkCellNameUnique(s(1).subs{1}, 1);           
            rowIndices = getDimensionIndices(a,s(1).subs{1}, 1, 0);
            if any(rowIndices == 0) && ~islogical(rowIndices)
                rowSubNames = s(1).subs{1};
                error(message('bioinfo:DataMatrix:subsref:UnrecognizedRowName', commaSeparatedList( rowSubNames( rowIndices==0 ) )))  
            elseif isempty(rowIndices) && ~isempty(a.RowNames)
                error(message('bioinfo:DataMatrix:subsref:UnrecognizedRowName', commaSeparatedList( s( 1 ).subs{ 1 } )))
            end
            
            %== Return if rowIndices is empty -[] case
            if issqrbrkt(s(1).subs{2})
                b = subsref(a.Matrix,s);
                varargout{1} = bioma.data.DataMatrix(b);
                return;
            end
            
            % Translate column names into indices (translates ':')
            checkCellNameUnique(s(1).subs{2}, 2);
            colIndices =  getDimensionIndices(a, s(1).subs{2}, 2, 0);
            if any(colIndices == 0) && ~islogical(colIndices)
                colSubNames = s(1).subs{2};
                error(message('bioinfo:DataMatrix:subsref:UnrecognizedColumnName', commaSeparatedList( colSubNames( colIndices==0 ) )))  
            elseif isempty(colIndices) && ~isempty(a.ColNames)
                error(message('bioinfo:DataMatrix:subsref:UnrecognizedColumnName', commaSeparatedList( s( 1 ).subs{ 2 } )))
            end
            
            % Create the output DataMatrix and copy non-indexed properties
            % over.
            if ~isempty(a.RowNames)
                rowNames = a.RowNames(rowIndices);
            else
                rowNames = [];
            end
            if ~isempty(a.ColNames)
                colNames = a.ColNames(colIndices);
            else
                colNames = [];
            end
            b = bioma.data.DataMatrix(a.Matrix(rowIndices, colIndices),...
                rowNames, colNames);
            
            varargout{1} = b;
        catch ME
            if strfind(ME.identifier, 'IndexOutOfBounds')
                error(message('bioinfo:DataMatrix:subsref:IndexOutOfBounds', inputname( 1 ), size( a.Matrix, 1 ), size( a.Matrix, 2 )));
            else
                ME.rethrow;
            end
        end
    case '{}'
        error(message('bioinfo:DataMatrix:subsref:CellSubscript'));
    case '.'
        % Using '.' to reference to elements in a DataMatrix      
        rowSubName = s(1).subs;
        %== Get row indices first. If is char, getDimensionIndices returns
        %empty
        % Translate row names into indices (leaves ':' alone)
        checkCellNameUnique(rowSubName, 1);
        rowIndices =  getDimensionIndices(a, rowSubName, 1, 0);       
        if ~isempty(rowIndices) % Indexing to rowindices
            % Using cell strings, if there are unrecognized names, it
            % should error:
            if any(rowIndices == 0) && ~islogical(rowIndices)
                ukNames = rowSubName(rowIndices == 0);
                error(message('bioinfo:DataMatrix:subsref:UnrecognizedRowName', commaSeparatedList( ukNames ))) 
            end
            
            b = a.Matrix(rowIndices, :);
            if ~isscalar(s)
                if size(b, 1)>1 && isscalar(s(2).subs)
                    s(2) = getSubscripts(a, s(2), 2, true);
                else
                    s(2) = getSubscripts(a, s(2), 2, false);
                end
         
                try
                    [varargout{1:nargout}] = subsref(b,s(2:end));
                catch ME
                    ME.rethrow;
                end
            else
                varargout{1} = b;
            end
            return;
        else % Find out if it is indexing into properties
            try
                [varargout{1:nargout}] = builtin('subsref',a,s);
            catch ME
                bioinfoprivate.bioclsrethrow('DataMatrix','subsref', ME)
            end
        end
end % switch
end % DataMatrix/subsref method

function s = getSubscripts(obj, s, dim, rows)
% Get the corresponding indices from obj, and replace the subscripts.
if nargin < 4
    rows = false;
end
if isequal(s.type, '()')
    indices = s.subs{1};
    checkCellNameUnique(indices, dim);
    indices =  getDimensionIndices(obj, indices, dim);
elseif isequal(s.type, '{}')
    checkCellNameUnique(s.subs, dim);
    indices =  getDimensionIndices(obj, s.subs, dim);
else isequal(s(2).type, '.')   
    checkCellNameUnique(s.subs, dim);
    indices =  getDimensionIndices(obj, s.subs, dim);
end
s.type = '()';
if rows
    s.subs = {':', indices};
else
    s.subs = {indices};
end
end % getSubsscripts

function tf = issqrbrkt(a)
%ISSQRBRKT Check if an array is the square bracket '[]'

% This checking does not work for the ones(0,0) like cases
tf = ~ischar(a)&&~isnumeric(a)&&isempty(a)&&~isobject(a)&&~isstruct(a);
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
            error(message('bioinfo:DataMatrix:subsref:NonUniqueRowNames'));
        case 'Column'
            error(message('bioinfo:DataMatrix:subsref:NonUniqueColumnNames'));
    end
end
end
