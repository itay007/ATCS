function varargout = colnames(obj,idx, names)
%COLNAMES Retrieve or set the column names of a DataMatrix object.
% 
%   RN = COLNAMES(A) returns the column names of a DataMatrix object A. 
%
%   RN = COLNAMES(A,I) returns the names of columns specified by I of a
%   DataMatrix object A. I can be a positive integer, a vector of positive
%   integers, a string specifying a column name, a cell array containing
%   one or more column names, or a logical vector.
%
%   B = COLNAMES(A,I,NAMES) returns a DataMatrix object B with names of
%   specified columns set to NAMES. The number of the names in NAMES must
%   match the number of columns specified by I. NAMES can be a cell array
%   of strings, a character array, or a numeric vector. NAMES can also be a
%   single string as a prefix for column names; column numbers will be
%   appended to the prefix. NAMES can also be a logical true or false
%   (default); if true, default unique names will be assigned to the
%   columns.
%
%   See also DATAMATRIX/ROWNAMES.

%   Copyright 2008-2012 The MathWorks, Inc. 


%== Input check
bioinfochecknargin(nargin,1,['DataMatrix:' mfilename])

if nargin == 1
    varargout{1} = obj.ColNames;
elseif nargin == 2
    if isempty(idx)
        idx = ':';
    end
    
    %== Check an get the column index
    colIndices = getDimensionIndices(obj, idx, 2, 0);
    
    if any(colIndices == 0) && ~islogical(colIndices)
        error(message('bioinfo:DataMatrix:colnames:UnrecognizedColumnName', commaSeparatedList( idx( colIndices==0 ) )))
    elseif isempty(colIndices) && ~isempty(obj.ColNames)
        error(message('bioinfo:DataMatrix:colnames:UnrecognizedColumnName', commaSeparatedList( idx )))
    end
    
    %== with only two input argument return the column names
    if ~isempty(obj.ColNames)
        varargout{1} = obj.ColNames(colIndices);
    else
        varargout{1} = [];
    end
elseif  nargin == 3 %== set the column names
    if isempty(idx) || bioma.util.isColon(idx)
        varargout{1:nargout} = set(obj, 'ColNames', names);
    else
        [colIndices, numIndices] = getDimensionIndices(obj, idx, 2, 0);
        
        if any(colIndices == 0) && ~islogical(colIndices)
            error(message('bioinfo:DataMatrix:colnames:UnrecognizedColumnName', commaSeparatedList( idx( colIndices==0 ) )))
        elseif isempty(colIndices) && ~isempty(obj.ColNames)
            error(message('bioinfo:DataMatrix:colnames:UnrecognizedColumnName', commaSeparatedList( idx )))
        end
        
        %== Check the names
        if isempty(names) 
            if numIndices == obj.NCols
                obj.ColNames = [];
            else
                obj.ColNames(colIndices) = cellstr(repmat('', 1, numIndices));
            end
            varargout{1:nargout} = obj;
            return;
        end
        
        if (iscellstr(names) || isnumeric(names) || ischar(names)) && isvector(names)
            if ischar(names)
                % returns a cell array of index strings
                if islogical(colIndices)
                    widx = (1:obj.NCols)';
                    idxStr = num2str(widx(colIndices));
                else
                    idxStr = num2str(colIndices(:));
                end
                names = cellstr([ repmat(names, numIndices,1)  idxStr]);
            elseif isnumeric(names)
                names = cellstr(num2str(names(:)));
            end
            
            if numel(names) ~= numIndices
                error(message('bioinfo:DataMatrix:colnames:NumIndicesNotMatchNumNames'))
            end
        elseif islogical(names)
            %default unique names.
                % returns a cell array of index strings
                if islogical(colIndices)
                    widx = (1:obj.NCols)';
                    idxStr = num2str(widx(colIndices));
                else
                    idxStr = num2str(colIndices(:));
                end
            names = cellstr([ repmat('Col', numIndices,1)  idxStr]);
        else
            error(message('bioinfo:DataMatrix:colnames:InvalidColNames'))
        end
        
        obj.ColNames(colIndices) = names;
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
