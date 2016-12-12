function varargout = rownames(obj,idx, names)
%ROWNAMES Retrieve or set the row names of a DataMatrix object.
% 
%   RN = ROWNAMES(A) returns the row names of a DataMatrix object A. 
%
%   RN = ROWNAMES(A,I) returns the names of rows specified by I of a
%   DataMatrix object A. I can be a positive integer, a vector of positive
%   integers, a string specifying a row name, a cell array containing one
%   or more row names, or a logical vector.
%
%   B = ROWNAMES(A,I,NAMES) returns a DataMatrix object B with names of
%   specified rows set to NAMES. The number of the names in NAMES must
%   match the number of rows specified by I. NAMES can be a cell array of
%   strings, a character array, or a numeric vector. NAMES can also be a
%   single string as a prefix for row names; row numbers will be appended
%   to the prefix. NAMES can also be a logical true or false (default); if
%   true, default unique names will be assigned to the rows.
%
%   See also DATAMATRIX/COLNAMES.

%   Copyright 2008-2012 The MathWorks, Inc. 


%== Input check
bioinfochecknargin(nargin,1,'DataMatrix:rownames')

if nargin == 1
    varargout{1} = obj.RowNames;
elseif nargin == 2
    if isempty(idx)
        idx = ':';
    end
    
    %== Check an get the row index
    rowIndices = getDimensionIndices(obj, idx, 1, 0);
    
    if any(rowIndices == 0) && ~islogical(rowIndices)
        error(message('bioinfo:DataMatrix:rownames:UnrecognizedRowName', commaSeparatedList( idx( rowIndices==0 ) )))
    elseif isempty(rowIndices) && ~isempty(obj.RowNames)
        error(message('bioinfo:DataMatrix:rownames:UnrecognizedRowName', commaSeparatedList( idx )))
    end
    
    %== with only two input argument return the row names
    if ~isempty(obj.RowNames)
        varargout{1} = obj.RowNames(rowIndices);
    else
        varargout{1} = [];
    end
elseif  nargin == 3 %== set the row names
    if isempty(idx) || bioma.util.isColon(idx)
        varargout{1:nargout} = set(obj, 'RowNames', names);
    else
        [rowIndices, numIndices] = getDimensionIndices(obj, idx, 1, 0);
        
        if any(rowIndices == 0) && ~islogical(rowIndices)
            error(message('bioinfo:DataMatrix:rownames:UnrecognizedRowName', commaSeparatedList( idx( rowIndices==0 ) )))
        elseif isempty(rowIndices) && ~isempty(obj.RowNames)
            error(message('bioinfo:DataMatrix:rownames:UnrecognizedRowName', commaSeparatedList( idx )))
        end
        
        %== Check the names
        if isempty(names)
            if numIndices == obj.NRows
                obj.RowNames = [];
            else
                obj.RowNames(rowIndices) = cellstr(repmat('', numIndices,1));
            end
            varargout{1:nargout} = obj;
            return;
        end
        
        if (iscellstr(names) || isnumeric(names) || ischar(names)) && isvector(names)
            if ischar(names)
                % returns a cell array of index strings
                if islogical(rowIndices)
                    widx = (1:obj.NRows)';
                    idxStr = num2str(widx(rowIndices));
                else
                    idxStr = num2str(rowIndices(:));
                end
                names = cellstr([ repmat(names, numIndices,1)  idxStr]);
            elseif isnumeric(names)
                names = cellstr(num2str(names(:)));
            end
            
            if numel(names) ~= numIndices
                error(message('bioinfo:DataMatrix:rownames:NumIndicesNotMatchNumNames'))
            end
        elseif islogical(names)
            %default unique names.
            % returns a cell array of index strings
            if islogical(rowIndices)
                widx = (1:obj.NRows)';
                idxStr = num2str(widx(rowIndices));
            else
                idxStr = num2str(rowIndices(:));
            end
            names = cellstr([ repmat('Row', numIndices,1)  num2str(idxStr)]);
        else
            error(message('bioinfo:DataMatrix:rownames:InvalidRowNames'))
        end
        
        obj.RowNames(rowIndices) = names;
        varargout{1:nargout} = obj;
    end
end
end % DataMatrix/rownames

function list = commaSeparatedList(list)
% Returns a comma separated list as a single string
if iscellstr(list)
   list = [sprintf('%s, ',list{1:end-1}) list{end}];
end
end
