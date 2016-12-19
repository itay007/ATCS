function [data, rowNames, colNames, varargout] = tableFileRead(file, varargin)
%TABLEFILEREAD Read table-oriented text files.
%
%   [DATA, ROWNAMES, COLNAMES] = TABLEFILEREAD(FILE) reads a table-oriented
%   text file, FILE, using IMPORTDATA and returns the values in a numeric
%   array DATA, and row names in ROWNAMES and column names in COLNAMES.
% 
%   TABLEFILEREAD(..., 'DELIMITER', DELIMITER) uses the specified
%   delimiter. DELIMITER can be any symbol. '\t' is the default.
% 
%   TABLEFILEREAD(..., 'HLINE', HLINE) where HLINE is a number indicating
%   which row of the file the column header text (column names) is located,
%   loads data from row HLINE+1 to the end of the file. The default HLINE =
%   1. If the data file contains no column names, set HLINE = 0.
% 
%   TABLEFILEREAD(..., 'ROWS', ROWNAMES, ...) reads only from a subset of
%   the rows, specified by ROWNAMES, in the file. ROWNAMES can be a cell
%   array of strings, a character array, or a numeric or logical vector.
% 
%   TABLEFILEREAD(..., 'COLUMNS', COLNAMES, ...) reads only from a subset
%   of the columns, specified by COLNAMES, in the file. COLNAMES can be a
%   cell array of strings, a character array, or a numeric or logical
%   vector.

%   Copyright 2008-2012 The MathWorks, Inc.


%== Check input arguments
inStruct = parse_inputs(varargin{:});
try 
    a = importdata(file, inStruct.Delimiter, inStruct.HLine);
catch ME 
    bioinfoprivate.bioclsrethrow('DataMatrix', mfilename, ME)
end

if isnumeric(a)
    data = a;
    rowNames = [];
    colNames = [];
elseif iscell(a) %if hline incorrect, a is a cell array of lines read.
     error(message('bioinfo:DataMatrix:tableFileRead:WrongFormat'))
elseif isstruct(a)
    if isempty(a.data)
        error(message('bioinfo:DataMatrix:tableFileRead:EmptyData'))
    end
    
    try
        %== parse the output from importdata
        [nrowD, ncolD] = size(a.data);
        if ~isfield(a, 'textdata')
            data = a.data;
            rowNames = [];
            colNames = [];
            return;
        end
        % Special case for multi-sheet Excel files
        if isstruct(a.data) && isstruct(a.textdata)
            [~,~,ext] = fileparts(file);
            if strcmpi(ext,'.xls')
                dataFieldnames = fieldnames(a.data);
                textFieldnames = fieldnames(a.textdata);
                if isequal(dataFieldnames,textFieldnames)
                   if numel(dataFieldnames) > 1
                       warning(message('bioinfo:DataMatrix:tableFileRead:MultiSheet'))
                   end
                   a.data = a.data.(dataFieldnames{1});
                   a.textdata = a.textdata.(dataFieldnames{1});
                   [nrowD, ncolD] = size(a.data);
                end
            end
        end
                   
                   
        %== parse cases
        [nrowTD, ncolTD] = size(a.textdata);
        
        %===No column names, with row names and with column+row names but
        %with delimiters other than '\t' 
        if nrowTD-nrowD == inStruct.HLine && ncolTD == 1
            data = a.data;
            rowNames = a.textdata(inStruct.HLine+1:end,1);
            if inStruct.HLine == 0
                colNames = [];
            else
                colNames = a.textdata(inStruct.HLine,:);
            end
         %===With column names, not row names
        elseif nrowTD == inStruct.HLine && ncolTD == ncolD
            data = a.data;
            rowNames = [];
            if inStruct.HLine == 0
                colNames = [];
            else
                colNames = a.textdata(inStruct.HLine, :);
            end
        %=== With both column names and row names
        elseif nrowTD-nrowD == inStruct.HLine && ncolTD-ncolD == 1
            data = a.data;
            rowNames = a.textdata(inStruct.HLine+1:end,1);
            if inStruct.HLine == 0
                colNames = [];
            else
                colNames = a.textdata(inStruct.HLine,2:end);
            end
        %===No column names and no row names with header
        elseif nrowD ~= nrowTD && nrowTD == inStruct.HLine && ncolTD == 1
            data =  a.data;
            rowNames = [];
            colNames = [];
        end
      
        if iscell(rowNames) && size(rowNames, 1) == 1 && size(data,1) >1
            rowNames = regexp(rowNames{:}, [inStruct.Delimiter '+'], 'split');
        end
        
        if iscell(colNames) && size(colNames, 2) == 1 && size(data,2) >1
            colNames = regexp(colNames{:}, [inStruct.Delimiter '+'], 'split');
        end
        
        if ischar(rowNames)
            rowNames = cellstr(rowNames);
        end
        if ischar(colNames)
            colNames = cellstr(colNames);
        end
        %==
        if ~isempty(rowNames) && size(data, 1) ~= numel(rowNames) 
            rowNames = [];
        end
        
        %== Case with importdata read the row name column as an empty
        %column name
        if ~isempty(rowNames) && (numel(colNames) - size(data, 2)) == 1
            colNames = colNames(2:end);
        end
        
        if ~isempty(colNames) && size(data, 2) ~= numel(colNames) 
            colNames = [];
        end
        
        %== parse rowNames and colNames for DataMatrix
        if ~isempty(inStruct.Rows)
            specIdx = getDimIndices(inStruct.Rows, rowNames, 1);
            rowNames = rowNames(specIdx);
            data = data(specIdx, :);
        end
        if ~isempty(inStruct.Columns)
            specIdx = getDimIndices(inStruct.Columns, colNames, 2);
            colNames = colNames(specIdx);
            data = data(:, specIdx);
        end
    catch ME
        bioinfoprivate.bioclsrethrow('DataMatrix', mfilename,ME);
    end
end

if nargout > 3
    if isempty(inStruct.Name)
        [~, varargout{1}] = fileparts(file);
    else
        varargout{1} = inStruct.Name;
    end
end
end %DataMatrix/tableFileRead

function inputStruct = parse_inputs(varargin)
% Check for the right number of inputs
if rem(nargin, 2)== 1
    error(message('bioinfo:DataMatrix:DataMatrix:IncorrectNumberOfArguments'))
end

% The allowed inputs
okargs = {'name','rows','columns','delimiter', 'hline'};

% Defaults
inputStruct.Name = [];           % DataMatrix object name
inputStruct.Rows = [];          % rownames
inputStruct.Columns = [];       % columnnames
inputStruct.Delimiter = '\t';   % tab-delimited
inputStruct.HLine = 1;

for i=1:2:nargin
    pname = varargin{i};
    pval = varargin{i+1};
    k = find(strncmpi(pname, okargs,numel(pname)));
    if isempty(k)
        error(message('bioinfo:DataMatrix:DataMatrix:UnknownParameterName', pname));
    elseif length(k)>1
        error(message('bioinfo:DataMatrix:DataMatrix:AmbiguousParameterName', pname));
    else
        switch(k)
            case 1 %name
                if ischar(pval)
                    inputStruct.Name = pval;
                elseif iscellstr(pval)
                    inputStruct.Name  = pval{1};
                else
                    error(message('bioinfo:DataMatrix:tableFileRead:BadNameInputFormat'))
                end
            case 2 % rows
                inputStruct.Rows = parseDimNames(pval, 1);
            case 3 % columns
                inputStruct.Columns = parseDimNames(pval, 2);
            case 4 % delimiter
                inputStruct.Delimiter = pval;
            case 5 % hline
                if isnumeric(pval) && isscalar(pval)
                    inputStruct.HLine = pval;
                else
                    error(message('bioinfo:DataMatrix:tableFileRead:InvalidHLine'))
                end
        end %switch
    end
end %for
end % parse_inputs


function dimVal = parseDimNames(pval, dim)
% Parse the row/column names
if ischar(pval)
    dimVal = {pval};
elseif iscellstr(pval)
    dimVal= pval;
elseif (isnumeric(pval) || islogical(pval) )&& isvector(pval)
    dimVal = pval;
else
    switch dim
        case 1
            error(message('bioinfo:DataMatrix:tableFileRead:BadROWSNames'));
        case 2
            error(message('bioinfo:DataMatrix:tableFileRead:BadCOLUMNSNames'));
    end
end
end % parseDimNames

function dimidx = getDimIndices(specIndices, dimNames, dim)
% Return the specified row/column name indices
if isnumeric(specIndices) || islogical(specIndices)
    dimidx = specIndices;
elseif iscellstr(specIndices)
    dimidx = zeros(1,numel(specIndices));
    for i = 1:numel(dimidx)
        idx = find(strcmp(specIndices{i},dimNames));
        if isempty(idx)
            error(message('bioinfo:DataMatrix:tableFileRead:UnrecognizedName', specIndices{ i }));
        else
            dimidx(i) = idx;
        end
    end  
    
    if all(dimidx==0)
        switch dim
            case 1
                warning(message('bioinfo:DataMatrix:DataMatrix:NoMatchingRowNames'))
            case 2
                warning(message('bioinfo:DataMatrix:DataMatrix:NoMatchingColumnNames'))
        end
        dimidx = true(numel(dimNames), 1);
    end
end
end
