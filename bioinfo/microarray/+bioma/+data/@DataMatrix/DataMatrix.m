classdef DataMatrix
%DATAMATRIX Two dimensional numeric or logical array with row and column names.     
%   
%   DataMatrix is a data structure that encapsulates two dimensional
%   numeric or logical arrays together with row names and column names. It
%   can be considered as a table-format matrix that can be indexed with its
%   row and column names. Most of the basic operations are parallel those
%   of MATLAB numeric arrays.
% 
%   Genomic data from a microarray based experiment can be best represented
%   as a tables (2D-matrix) with each row representing a gene and each
%   column representing a sample. DataMatrix is a MATLAB data container for
%   storing microarray experiment measurement data.
% 
%   DataMatrix properties:
%       Name            - Name of the DataMatrix.
%       RowNames        - Names of the rows.
%       ColNames        - Names of the columns.
%       NRows           - Number of rows in a DATAMATRIX.
%       NCols           - Number of columns in a DATAMATRIX.
%       NDims           - Number of dimensions of a DATAMATRIX.
%       ElementClass    - Class name of the data elements in a DataMatrix.
% 
%   DataMatrix methods:
%       DataMatrix  - Create a DataMatrix object.
%       colnames    - Retrieve or set the column names of a DataMatrix.
%       rownames    - Retrieve or set the row names of a DataMatrix.
%       dmarrayfun  - Apply a function to each element of a DataMatrix.
%       dmbsxfun    - Binary singleton expansion function for DataMatrix.
%       dmwrite     - Writes DataMatrix to a text file.
%       
%   Note: The basic MATLAB array like operations are not listed here. For
%   more details, please see the reference page of DATAMATRIX.
% 
%   Example:
%       d =bioma.data.DataMatrix(rand(3,4),'Rownames',true,'Colnames','Test')
%       % Index the DataMatrix d
%       d(2:3, {'Test1','Test3'})
%       % Sum the DataMatrix d
%       sum(d)
%
%   See also DATASET, EXPRESSIONSET, EXPTDATA, METADATA.

%   Copyright 2008-20012 The MathWorks, Inc.


properties (SetAccess = 'public', GetAccess = 'public')
    %NAME Name of the DataMatrix.
    %    The Name property is a string specifying the name of a DataMatrix.
    %
    %    See also DATAMATRIX.
    Name = '';
    
    %ROWNAMES Names of the rows.
    %   The RowNames property is a cell array of strings containing the
    %   row names.
    %
    %   See also DATAMATRIX.
    RowNames =[];
    
    %COLNAMES Names of the columns.
    %   The ColNames property is a cell array of strings containing the
    %   column names.
    %
    %   See also DATAMATRIX.
    ColNames =[];
end

properties (SetAccess = 'private', GetAccess = 'public')
    %NROWS Number of rows in a DATAMATRIX.
    %    The NRows property is a scalar equals to the number of rows in
    %    a DataMatrix.
    %
    %    See also DATAMATRIX.
    NRows = 0;
    
    %NCOLS Number of columns in a DATAMATRIX.
    %    The NCols property is a scalar equals to the number of columns in
    %    a DataMatrix.
    %
    %    See also DATAMATRIX.
    NCols = 0;
    
    %NDIMS Number of dimensions of a DATAMATRIX.
    %    The NDims property is a scalar of the number of dimensions of a
    %    DataMatrix. It is a constant of 2.
    %
    %    See also DATAMATRIX.
    NDims = 2;
    
    %ELEMENTCLASS Class name of the data elements in a DataMatrix.
    %    The ElementClass property is string specifies the class type of
    %    the elements in a DataMatrix.
    %
    %    See also DATAMATRIX.
    ElementClass = [];
end

properties(SetAccess = 'private', GetAccess = 'protected', Hidden=true)
    Matrix = [];
end

methods
    %==Constructor
    function obj = DataMatrix(varargin)
        % DATAMATRIX Create a DataMatrix object.
        %
        %   DM = DATAMATRIX(DATA) creates an object, DM, of the DataMatrix
        %   class from DATA, a two-dimensional MATLAB numeric or logical
        %   array. DATA can also be a DataMatrix object.
        %
        %   DM = DATAMATRIX(DATA, ROWNAMES, COLUMNNAMES) creates an object,
        %   DM, of the DataMatrix class from DATA, a two-dimensional MATLAB
        %   numeric or logical array with specified row names ROWNAMES and
        %   column names COLUMNNAMES. The number of the names in ROWNAMES
        %   must match the number of rows in DATA. The number of names in
        %   COLUMNNAMES must match the number of columns in DATA. The row
        %   names and column names can be a cell array of strings, a
        %   character array, or a numeric vector. If ROWNAMES or
        %   COLUMNNAMES is empty or not specified, default unique names
        %   will be assigned to the rows or columns.
        %
        %   DM = DATAMATRIX(..., 'ROWNAMES',ROWNAMES) creates an object,
        %   DM, of the DataMatrix class with row names specified by
        %   ROWNAMES. The number of the names in ROWNAMES must match the
        %   number of rows in DATA. The row names can be a cell array of
        %   strings, a character array, or a numeric or logical vector.
        %   ROWNAMES can also be a single string as a prefix for row names;
        %   row numbers will be appended to the prefix. ROWNAMES can also
        %   be a logical true or false (default); if true, default unique
        %   names will be assigned to the rows.
        %
        %   DM = DATAMATRIX(..., 'COLNAMES',COLUMNNAMES) creates an object,
        %   DM, of the DataMatrix class with column names specified by
        %   COLUMNNAMES. The number of names in COLUMNNAMES must match the
        %   number of columns in DATA. The column names can be a cell array
        %   of strings, a character array, or a numeric or logical vector.
        %   COLUMNNAMES can also be a single string as a prefix for column
        %   names; column numbers will be appended to the prefix.
        %   COLUMNNAMES can also be a logical true or false (default); if
        %   true, default unique names will be assigned to the columns.
        %
        %   DM = DATAMATRIX(..., 'NAME', NAME) creates an object, DM, of
        %   the DataMatrix class with specified name NAME. The name must be
        %   a string.
        %
        %   DM = DATAMATRIX('FILE',FILENAME) creates a DataMatrix object,
        %   DM, using the IMPORTDATA function to read table-oriented data
        %   in a text file. The file typically consists of column names in
        %   the first row, row names in the first column, and numeric data
        %   starting in the (2,2) position, and uses '\t' as the delimiting
        %   character.
        %
        %   DM = DATAMATRIX('FILE',FILENAME, 'DELIMITER', DELIMITER, ...)
        %   uses the specified delimiter. DELIMITER can be any symbol. '\t'
        %   is the default.
        %
        %   DM = DATAMATRIX('FILE',FILENAME, 'HLINE', HLINE) where HLINE is
        %   a number indicating which row of the file the column header
        %   text (column names) is located, loads data from row HLINE+1 to
        %   the end of the file. The default HLINE = 1. If the data file
        %   contains no column names, set HLINE = 0.
        %
        %   DM = DATAMATRIX('FILE',FILENAME, 'ROWS', ROWNAMES, ...) creates
        %   an object, DM, of the DataMatrix class from a subset of the
        %   rows, specified by ROWNAMES, in the table-oriented data in a
        %   file. ROWNAMES can be a cell array of strings, a character
        %   array, or a numeric or logical vector.
        %
        %   DM = DATAMATRIX('FILE',FILENAME, 'COLUMNS', COLNAMES, ...)
        %   creates an object, DM, of the DataMatrix class from a subset of
        %   the columns, specified by COLNAMES, in the table-oriented data
        %   in a file. COLNAMES can be a cell array of strings, a character
        %   array, or a numeric or logical vector.
        %
        %   See also DATAMATRIX.

        if nargin == 0
            return
        end
        
        %== Check if the first argument is matrix
        if bioma.util.isString(varargin{1})
            % File name case
            obj = createDataMatrixFromFile(obj, varargin{:});
        else
            obj = createDataMatrix(obj, varargin{:});
        end
    end
end
end % DataMatrix class

%== Helper functions
function obj = createDataMatrix(obj, varargin)
% Create a DataMatrix object with input data matrix
argCount = 1;
nInput = nargin -1;

%== The first argument must be a numeric or logical MATLAB array
arg = varargin{argCount};
if isa(arg, 'bioma.data.DataMatrix')
    obj = arg;
else
    try
        validateattributes(arg, {'numeric', 'logical'},{'real', '2d'});
    catch ME 
        error(message('bioinfo:DataMatrix:DataMatrix:InvalidArgumentMatrix'))
    end
    
    obj.Matrix = arg;
    obj.NRows = size(obj.Matrix, 1);
    obj.NCols = size(obj.Matrix, 2);
    obj.NDims = ndims(obj.Matrix);
    obj.ElementClass = class(obj.Matrix);
end

%== Processing individual input arguments
s_argCount = argCount;
while argCount < nInput
    argCount = argCount + 1;
    arg = varargin{argCount};
    %== Guess if the input is param name/value pairs
    if bioma.util.isString(arg)
        %== Start of param name/value pairs
        try
            if ~(ischar(arg) && isvector(arg) && (size(arg,1)==1))
                error(message('bioinfo:DataMatrix:DataMatrix:InvalidPropertyName'));
            end
            bioinfoprivate.pvpair(arg, [],[fieldnames(obj); {'Matrix'}], 'DataMatrix:DataMatrix');
            s_argCount = argCount;
            argCount = argCount - 1;
            break;
        catch ME %#ok
            %== Add the Names for dimensions
            if argCount == 2
                obj = set(obj, 'RowNames',arg);
            elseif argCount == 3
                obj = set(obj, 'ColNames', arg);
            else
                argCount = argCount - 1;
                s_argCount = argCount;
                break;
            end
        end
    else
        %== Add the Names for dimensions
        if argCount == 2
            obj = set(obj, 'RowNames',arg);
        elseif argCount == 3
            obj = set(obj, 'ColNames', arg);
        end
        s_argCount = argCount;
    end
end % while argCount < nargin processing individual input argument

%==Processing param name/value pair
if argCount < nInput
   obj = parse_inputs(obj, varargin{s_argCount:end});
end
end % createDataMatrix

function obj = parse_inputs(obj, varargin)
% Parse input PV pairs.
if nargin < 2
    return;
end
% Check for the right number of inputs
if rem(nargin-1, 2)== 1
    error(message('bioinfo:DataMatrix:DataMatrix:IncorrectNumberOfArguments'))
end

obj = set(obj, varargin{:});
end % parse_inputs

function obj = createDataMatrixFromFile(obj, varargin)
% Create DataMatrix object by reading data from a file
argCount = 1;
nInput = nargin -1;

%== The first argument must be a string "file"
arg = varargin{argCount};
try
    validatestring(arg, {'File', 'file'});
catch ME 
    error(message('bioinfo:DataMatrix:DataMatrix:InvalidArgumentFile'))
end

% Check for the right number of inputs
if rem(nInput, 2)== 1
    error(message('bioinfo:DataMatrix:DataMatrix:IncorrectNumberOfArguments'))

end
argCount = argCount + 1;
%== Get the file names
arg = varargin{argCount};
[data, rowNames, colNames, name] = tableFileRead(arg, varargin{argCount+1:end});
%== Create the data matrix
obj = createDataMatrix(obj, data, rowNames, colNames, 'Name', name);
end % createDataMatrixFromFile
