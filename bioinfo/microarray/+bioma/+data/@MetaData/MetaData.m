classdef MetaData
    %METADATA Class to contain metadata collection of variables and their description.
    %
    %   The MetaData class consists of two parts: a dataset array of
    %   samples (rows or observations) and the values of variables measured
    %   on those samples; and a dataset array containing a description and
    %   other meta information of each variable measured.
    %
    %   In the dataset array containing the measurement values, each column
    %   corresponds to a variable and each row corresponds to a sample
    %   (row or observation). While in the variable meta information
    %   dataset, each row corresponds to a variable, and each column
    %   corresponds to a meta information labels. At least one column,
    %   named VariableDescription, contains a textual description of each
    %   variable. To access the dataset arrays in a MetaData object , you
    %   can use the methods VARIABLEVALUES and VARIABLEDESC.
    %
    %   MetaData properties:
    %       Name             - Name of the MetaData object.
    %       Description      - Textual description of the MetaData object.
    %       NSamples         - Number of samples (rows or observations).
    %       NVariables       - Number of variables measured.
    %       DimensionLabels  - Dimension labeling for rows and columns.
    %
    %   MetaData methods:
    %       MetaData           - Create a MetaData object.
    %       variableValues     - Access and set the measured values of the variables.
    %       variableDesc       - Access and set the meta information about the variables.
    %       variableNames      - Access and set the variable names.
    %       sampleNames        - Access and set the sample names.
    %       combine            - Combine data from two MetaData objects.
    %       isempty            - True for empty MetaData object.
    %       size               - Size of a MetaData object.
    %       varValuesTable     - Create a graphic table showing variable values.
    %
    %   Examples:
    %       load fisheriris
    %       iris = dataset({nominal(species),'species'},{meas,'SL','SW','PL','PW'});
    %       vardesc ={'Iris species', 'Sepal Length', 'Sepal Width',...
    %                                 'Petal Length', 'Petal Width'}';
    %       variris = dataset(vardesc, 'VarNames', {'VariableDescription'},...
    %                                  'ObsNames', {'species', 'SL', 'SW', 'PL', 'PW'});
    %       a = bioma.data.MetaData(iris, variris);
    %       % From a text file
    %       b = bioma.data.MetaData('File', 'sampleMetaData.txt');
    %
    %   See also DATAMATRIX, DATASET, EXPRESSIONSET, EXPTDATA, IMPORTDATA, MIAME. 
    
    %   Copyright 2009-2012 The MathWorks, Inc.

    
    properties(SetAccess = 'public', GetAccess='public', Hidden=false)
        %NAME Name of the MetaData object.
        %
        %    The Name property is a string specifying the name of a
        %    MetaData object.
        %
        %    See also METADATA.
        Name = '';
        
        %DESCRIPTION Textual description of the MetaData object.
        %
        %    The Description property is a string of a description about a
        %    MetaData object.
        %
        %    See also METADATA.
        Description = '';
        
        %DIMENSIONLABELS Dimension labeling for rows and columns.
        %
        %    The DimensionLabels property is a two-element cell array of
        %    strings containing the names of the dimensions of a MetaData
        %    object. The default labels are Samples for dimension 1 (row)
        %    and Variables for dimension 2 (column).
        %
        %    See also METADATA.
        DimensionLabels = {'Samples', 'Variables'};
    end
    
    properties (SetAccess = 'private', GetAccess = 'public')
        %NSAMPLES Number of samples (rows).
        %
        %    The NSamples property is a scalar equal to the number of
        %    samples (observations or rows) in a MetaData object.
        %
        %    See also METADATA.
        NSamples = 0;
        
        %NVARIABLES Number of variables (columns).
        %
        %    The NVariables property is a scalar equal to the number of
        %    variables (columns) in a MetaData object.
        %
        %    See also METADATA.
        NVariables = 0;
    end
    
    properties(SetAccess = 'private', GetAccess = 'protected', Hidden=true)
        %DATA The metadata in the MetaData.
        %    The Data property is a dataset containing samples (rows) and
        %    the values of the measured variables (columns) on the samples.
        Data = dataset;
        
        %VARMETAINFO Variable meta information in the MetaData.
        %   The VarMetaInfo property is a dataset containing the meta
        %   information about the measured variables in the MetaData. At
        %   least one column, named VariableDescription, containing a
        %   textual description of each variable. The number of rows in
        %   VarMetaInfo dataset is equal to the number of columns in Data
        %   property.
        VarMetaInfo = dataset;
    end
    
    methods
        function obj = MetaData(varargin)
            %METADATA Create a MetaData object.
            %
            %   S = METADATA(DATA) creates an object, S, of MetaData class
            %   from a dataset array DATA. The dataset DATA contains a
            %   collection of samples (observations) and the values of
            %   measured variables with each column corresponding to a
            %   variable and each row corresponding to a sample
            %   (observation).
            %
            %   S = METADATA(DATA, VARMETAINFO) creates a MetaData object S
            %   from dataset arrays DATA and VARMETAINFO. VARMETAINFO
            %   contains meta information about the variables in DATA. In
            %   the VARMETAINFO dataset, the variables are the rows. The
            %   number of rows in the VARMETAINFO dataset must be equal to
            %   the number of columns in DATA. At least one column, named
            %   VariableDescription, contains a textual description of each
            %   variable.
            %
            %   S = METADATA(DATA, VARDESC) creates a MetaData object S
            %   from a dataset array DATA and VARDESC, a cell array of
            %   strings that describe the variables in DATA. The length of
            %   VARDESC must equal the number of variables in DATA.
            %
            %   S = METADATA('FILE',FILENAME) creates a MetaData object S
            %   from data imported from text file, FILENAME, containing a
            %   table with mixed numeric and text columns.
            %
            %   S = METADATA('FILE',FILENAME,..., 'PATH', FILEPATH)
            %   specifies the path to FILENAME.
            %
            %   S = METADATA('FILE',FILENAME,..., 'DELIMITER', DELIMITER)
            %   specifies the column separator. It must be a string.
            %   Default is '\t' for tab.
            %
            %   S = METADATA('FILE',FILENAME,..., 'ROWNAMES', ROWNAMES).
            %   ROWNAMES can be a cell array of strings giving the row
            %   names, or a single number giving the column of the table
            %   which contains the row names, or a character string giving
            %   the name of the table column containing the row names. By
            %   default the first column in the table is used for the row
            %   names. If ROWNAMES is set to empty or zero, the rows will
            %   be numbered.
            %
            %   S = METADATA('FILE',FILENAME,..., 'COLNAMES', COLNAMES)
            %   only reads the data from columns specified by the names in
            %   cell array COLNAMES. The default is to read data from all
            %   columns in the table, and the first row contains headers
            %   for the columns. Set COLNAMES to empty to specify column
            %   names are missing from the table. The column will be
            %   numbered.
            %
            %   S = METADATA('FILE',FILENAME,..., 'VARDESCCHAR', VC)
            %   specifies header lines beginning with character VC as
            %   variable name descriptions used for creating the
            %   VarMetaInfo dataset. By default no variable meta
            %   information will be read.
            %
            %   S = METADATA(..., 'NAME', NAME) creates a MetaData object
            %   S, with specified name NAME. NAME must be a string.
            %
            %   S = METADATA(..., 'DESCRIPTION', DESC) creates a MetaData
            %   object S, with specified description DESC. DESC must be a
            %   string.
            %
            %   S = METADATA(..., 'SAMPLENAMES',SNAMES) creates a MetaData
            %   object S, with specified the sample names. SNAMES can be a
            %   cell array of strings with length equal to the number of
            %   samples in S.
            %
            %   S = METADATA(..., 'VARIABLENAMES' VNAMES) creates a
            %   MetaData object S, with specified the variable names.
            %   VNAMES can be a cell array of strings with length equal to
            %   the number of variables in S.
            %
            %   See also METADATA, DATAMATRIX, DATASET, IMPORTDATA.
            
            %   S = METADATA('FILE',FILENAME,..., 'MAXHEADERLINES', MH)
            %   specifies maximum number of header lines in the file. The
            %   default is 100.
            
            
            %== Process inputs
            if nargin == 0
                return;
            end
            
            %= Check if the first argument is a dataset or a string
            if bioma.util.isString(varargin{1})
                % File name case
                obj = createMetadataSetFromFile(obj, varargin{:});
            else
                obj = createMetadataSet(obj, varargin{:});
            end
        end % MetaData
    end
    
    methods
        function data = variableValues(obj, value)
            %VARIABLEVALUES Access and set the measured values of the variables.
            %
            %   DS = VARIABLEVALUES(M) returns a dataset array, DS, with
            %   samples as rows (observations), variables as columns from a
            %   MetaData object M.
            %
            %   B = VARIABLEVALUES(M, DS) returns a MetaData object B with
            %   the variable measurement values set to dataset array DS.
            %
            %   See also METADATA, SAMPLENAMES, VARIABLEDESC, VARIABLENAMES.
            
            %== Input check
            bioinfochecknargin(nargin, 1, 'MetaData:variableValues')
            if nargin == 1
                data = obj.Data;
            elseif nargin == 2
                if ~isa(value, 'dataset')
                    error(message('bioinfo:MetaData:variableValues:InvalidDataType'))
                end
                data = setData(obj, value);
            end
        end
        
        function data = variableDesc(obj, value)
            %VARIABLEDESC Access and set the meta information about the variables.
            %
            %   DS = VARIABLEDESC(M) returns a dataset array, DS, with
            %   variable names as rows, description labels as columns of
            %   MetaData object M.
            %
            %   B = VARIABLEDESC(M, DS) returns a MetaData object B with
            %   the variable meta information data set to a dataset array
            %   DS. DS can also be a cell array of strings containing
            %   descriptions of the variables for the MetaData.
            %
            %   See also METADATA, SAMPLENAMES, VARIABLEVALUES, VARIABLENAMES. 
            
            %== Input check
            bioinfochecknargin(nargin, 1, 'MetaData:variableDesc')
            if nargin == 1
                data = obj.VarMetaInfo;
            elseif nargin == 2
                try
                    data = setVarMetaInfo(obj, value);
                catch ME
                    bioinfoprivate.bioclsrethrow('MetaData', 'variableDesc', ME)
                end
            end
        end
        
        function out = variableNames(obj, idx, names)
            %VARIABLENAMES Access and set the variable names.
            %
            %   VNAMES = VARIABLENAMES(M) returns a cell array of variable
            %   names from MetaData object M.
            %
            %   VNAMES = VARIABLENAMES(M, I) returns the names of variables
            %   specified by I in MetaData object M. I can be a positive
            %   integer, a vector of positive integers, a string specifying
            %   a variable name, a cell array containing one or more
            %   variable names, or a logical vector.
            %
            %   B = VARIABLENAMES(M, I, NAMES) returns MetaData object B
            %   with the names of specified variables set to NAMES. The
            %   number of names in NAMES must equal the number of variables
            %   specified by I. NAMES can be a cell array of strings, a
            %   character array, or a numeric vector. NAMES can also be a
            %   single string as a prefix for variable names; variable
            %   numbers will be appended to the prefix. NAMES can also be a
            %   logical true or false (default); if true, default unique
            %   names will be assigned to the variables. The variable names
            %   (as rows) in variable meta information dataset array are
            %   updated automatically.
            %
            %   See also METADATA, SAMPLENAMES, VARIABLEVALUES, VARIABLEDESC. 

            %== Input check
            bioinfochecknargin(nargin, 1, 'MetaData:variableNames')
            try
                if nargin < 2
                    out = get(obj.Data, 'VarNames');
                elseif nargin < 3
                    out = bioma.util.nameList(obj, 'variableNames',idx);
                elseif nargin ==3
                    if isempty(idx) || bioma.util.isColon(idx)
                        out = setVariableNames(obj, names);
                    else
                        if ischar(names)
                            names = cellstr(names);
                        end
                        N = numel(idx);
                        if iscellstr(names)
                            if isempty(names) || numel(names) ~= N || (numel(unique(names))~=numel(names))
                                error(message('bioinfo:MetaData:variableNames:InvalidInput', 'VariableNames'));
                            end
                        else
                            error(message('bioinfo:MetaData:variableNames:InvalidInputType', 'VariableNames'));
                        end
                        tNames = variableNames(obj);
                        tNames(idx) = names;
                        out = setVariableNames(obj, tNames);
                    end
                end
            catch ME
                bioinfoprivate.bioclsrethrow('MetaData', 'variableNames', ME)
            end
        end
        
        function out = sampleNames(obj, idx, names)
            %SAMPLENAMES Access and set the sample names.
            %
            %   SNAMES = SAMPLENAMES(M) returns a cell array of sample
            %   (observation or row) names from MetaData object M.
            %
            %   SNAMES = SAMPLENAMES(M, I) returns the names of samples
            %   specified by I in MetaData object M. I can be a positive
            %   integer, a vector of positive integers, a string specifying
            %   a sample name, a cell array containing one or more sample
            %   names, or a logical vector.
            %
            %   B = SAMPLENAMES(M, I, NAMES) returns a MetaData object M
            %   with the names of specified samples set to NAMES. The
            %   number of names in NAMES must equal the number of samples
            %   specified by I. NAMES can be a cell array of strings, a
            %   character array, or a numeric vector. NAMES can also be a
            %   single string as a prefix for sample names; sample numbers
            %   will be appended to the prefix. NAMES can also be a logical
            %   true or false (default); if true, default unique names will
            %   be assigned to the samples.
            %
            %   See also METADATA, VARIABLEVALUES, VARIABLEDESC, VARIABLENAMES. 
            
            %== Input check
            bioinfochecknargin(nargin, 1, 'MetaData:sampleNames')
            try
                if nargin < 2
                    out = get(obj.Data, 'ObsNames');
                elseif nargin < 3
                    out = bioma.util.nameList(obj, 'sampleNames',idx);
                elseif nargin ==3
                    if isempty(idx) || bioma.util.isColon(idx)
                        out = setSampleNames(obj, names);
                    else
                        if ischar(names)
                            names = cellstr(names);
                        end
                        N = numel(idx);
                        if iscellstr(names)
                            if isempty(names) || numel(names) ~= N || (numel(unique(names))~=numel(names))
                                error(message('bioinfo:MetaData:sampleNames:InvalidInput', 'SampleNames'));
                            end
                        else
                            error(message('bioinfo:MetaData:sampleNames:InvalidInputType', 'SampleNames'));
                        end
                        
                        tNames = sampleNames(obj);
                        tNames(idx) = names;
                        out = setSampleNames(obj, tNames);
                    end
                end
            catch ME
                bioinfoprivate.bioclsrethrow('MetaData', 'sampleNames', ME)
            end
        end
        
        function c = combine(a, b)
            %COMBINE Combine data from two MetaData object.
            %
            %   C = COMBINE(A, B) combines data from MetaData object A with
            %   MetaData object B and returns the result in MetaData object
            %   C. Sample names in A and B must be unique. Variable names
            %   present in both A and B occupy a single column in the
            %   resulting MetaData object C. Variable names unique to
            %   either A or B create columns with values to those samples
            %   where the variable is present. Variable meta inforamtion
            %   data in the resulting MetaData object C is updated to
            %   reflect the combine.
            %
            %   See also METADATA.
            
            bioinfochecknargin(nargin, 2, 'MetaData:combine')
            
            if ~isa(a, 'bioma.data.MetaData') || ~isa(b, 'bioma.data.MetaData')
                error(message('bioinfo:MetaData:combine:InvalidInput'))
            end
            %== Find the common variable names in A and B
            try
                data = horzCombineDataset(a.Data, b.Data);
                vmi = vertCombineDataset(a.VarMetaInfo, b.VarMetaInfo);
                c = bioma.data.MetaData(data, vmi);
            catch ME
                bioinfoprivate.bioclsrethrow('MetaData', 'combine', ME)
            end
        end
        
        function varargout = varValuesTable(obj, varargin)
            %VARVALUESTABLE Create a graphic table showing variable values.
            %
            %   H = VARVALUESTABLE(M) creates a two dimensional graphic
            %   table with variable values of a MetaData object M and
            %   returns the uitable handle H.
            %
            %   VARVALUESTABLE(DM, FH) specifies the parent handle of the
            %   uitable object. The parent can be a figure or uipanel
            %   handle.
            %
            %   See also METADATA.
            
            %== Check input arguments
            if nargin > 1
                arg = varargin{1};
                if ishandle(arg)
                    hParent = arg;
                end
            else
                figureName = bioinfoprivate.indexedFigureName('MetaData', 'MetaData');
                hParent = figure('Renderer',' zbuffer',...
                    'Name', figureName,...
                    'NumberTitle','off',...
                    'Tag', 'MetaData',...
                    'IntegerHandle','off',...
                    'Visible', 'on');
            end
            
            hTable = uitable('Parent', hParent);
            set(hParent, 'ResizeFcn', {@bioma.util.parentResizeCB,...
                hTable, 1/50, 1/50});
            
            %== Show data
            set(hTable, 'Data', datasetToCell(obj.Data))
            set(hTable, 'ColumnName', get(obj.Data, 'VarNames'),...
                'RowName', get(obj.Data, 'ObsNames'));
            %== Output
            if nargout >  0
                varargout{1} = hTable;
            end
        end % end of varValuesTable
        
        function t = isempty(obj)
            %ISEMPTY True for empty MetaData object.
            %
            %   TF = ISEMPTY(M) returns true (1) if M is an empty MetaData
            %   object and false (0) otherwise. An empty MetaData object
            %   has empty dataset arrays.
            %
            %   See also METADATA.
            
            t = isempty(obj.Data);
        end
        
        function varargout = size(a,dim)
            %SIZE Size of a MetaData object.
            %
            %   D = SIZE(M) returns a two-elment row vector D containing
            %   the number of samples and variables in MetaData object M.
            %
            %   [M,N] = SIZE(M) returns the number of samples and variables
            %   in MetaData object M as separate output variables.
            %
            %   D = SIZE(M, DIM) returns the length of the dimension
            %   specified by the scalar DIM. For example, SIZE(M,1) returns
            %   the number of samples, and SIZE(M,2) returns the number of
            %   variables. If DIM > NDIMS(M), D will be 1.
            %
            %   See also METADATA.
            
             if nargin == 1
                if nargout < 2
                    varargout{:} = size(a.Data);
                elseif nargout == 2
                    [varargout{1}, varargout{2}] = size(a.Data);
                else
                    [varargout{1}, varargout{2}] = size(a.Data);
                    varargout(3:nargout) = {1};
                end
                
            else
                if nargout > 1
                    error(message('bioinfo:MetaData:size:TooManyOutputs'));
                end
                
                try
                    varargout{1} = size(a.Data, dim);
                catch ME
                    bioinfoprivate.bioclsrethrow('MetaData','size', ME)
                end
             end
        end % size method
    end % end of methods
    
    methods(Hidden = true)
        function [varargout] = subsref(a, s)
            %SUBSREF Subscripted reference for a MetaData object.
            %
            %   B = A(I,J) returns a MetaData that contains a subset of the
            %   rows and columns in the MetaData A. I and J can be a
            %   positive integer, a vector of positive integers, a string
            %   of names, a cell array containing one or more unique names,
            %   or a logical vector. B contains the same property values as
            %   A, subsetted for rows or columns where appropriate. ()
            %   subscribe type does not allow linear indexing.
            %
            %   Note: Cell indexing A{I,J} is not supported.
            %
            %   P = A.PROPERTYNAME returns a MetaData property. 
            
            % Undocumented:
            % D = A.VARIABLENAME returns the data for named columns.
            
            switch s(1).type
                case '()'
                    try
                        % Change subscript for VarMetaInfo
                        s_vm = s;
                        s_vm.subs = {s.subs{2},':'};
                        [varargout{1:nargout}] = bioma.data.MetaData(subsref(a.Data, s),...
                            subsref(a.VarMetaInfo, s_vm));
                    catch ME
                        bioinfoprivate.bioclsrethrow('MetaData','subsref', ME)
                    end
                case '{}'
                    error(message('bioinfo:MetaData:subsref:CellSubscript'));
                case '.'                  
                    queryName = s(1).subs;
                    % Check if it is a valid variable name (undocumented)
                    if any(strcmp(a.variableNames,queryName))
                        [varargout{1:nargout}] = subsref(a.Data,s);
                    % Check if it is a valid property name    
                    elseif any(strcmp(properties(a),queryName))
                        [varargout{1:nargout}] = builtin('subsref',a,s);
                    elseif any(strcmp(methods(a),queryName))
                        [varargout{1:nargout}] = builtin('subsref',a,s);
                    else
                        error(message('bioinfo:MetaData:subsref:NotAProperty'));
                    end                        
            end
        end %subsref
        
        function obj = subsasgn(obj,s,b)
            %SUBSASGN Subscripted assignment to a MetaData object.
            %
            %   A(I,J) = B assigns the contents of the MetaData B to a
            %   subset of the samples and variables in the MetaData A
            %   specified by the subscript vectors I and J. I and J can be
            %   a positive integer, a vector of positive integers, a string
            %   of name, a cell arrays containing one or more unique names,
            %   or a logical vector.  The assignment does not use row
            %   names, column names, or any other properties of B to modify
            %   properties of A; Elements of B are assigned into A by
            %   position, not by matching names. However properties of A
            %   are extended with default values if the assignment expands
            %   the number of observations or variables in A.
            %
            %   Note: Cell indexing A{I,J} = B is not supported.
            %
            %   A.PROPERTYNAME = P assigns to a MetaData property.
            
            %== Empty MetaData
            creating = isequal(obj,[]);
            if creating
                obj = bioma.data.MetaData;
            end
            
            switch s(1).type
                case '()'
                    try
                        % Change subscript for VarMetaInfo
                        tmp = obj;
                        tmp.Data = subsasgn(obj.Data, s, b.Data);
                        tmp.VarMetaInfo = subsasgn(obj.VarMetaInfo, s, b.VarMetaInfo);
                        obj = tmp;
                    catch ME
                        bioinfoprivate.bioclsrethrow('MetaData','subsref',ME)
                    end
                case '{}'
                    error(message('bioinfo:MetaData:subsasgn:CellSubscript'));
                case '.'
                    % Intercept only when it is a valid property name,
                    % otherwise error out.
                    propName = s(1).subs;
                    if any(strcmp(properties(obj),propName))
                        switch propName
                            case 'Name' % Names
                                obj = setName(obj,p);
                            case 'Description' % Description
                                obj = setDescription(obj, p);
                            case 'DimensionLabels' % DimemsionLabels
                                obj = setDimensionLabels(obj, p);
                            otherwise
                                error(message('bioinfo:MetaData:subsasgn:SetProhibited', propName));
                        end                                
                    else
                        error(message('bioinfo:MetaData:subsasgn:NotAProperty'));
                    end
            end %switch
        end %subsasgn
        
        function disp(a)
            %DISP Display MetaData object.
            %
            %   DISP(A) displays the MetaData A, including the first two
            %   and last row names, a table formatted variable labels and
            %   the variable metadata, without printing the MetaData object
            %   name.
            
            colPad = repmat(' ', 1, 4);
            rowNames = a.sampleNames;
            amChars = bioma.util.printOneLineMultiNames(rowNames, colPad,...
                [a.NSamples, a.NVariables], 'MetaData');
            
            %== Display row names
            fprintf('Sample Names:\n')
            fprintf('%s\n', amChars)
            
            fprintf('Variable Names and Meta Information:\n')
            if isempty(a.VarMetaInfo)
                varNames = get(a.Data, 'VarNames');
                
                for i = 1:a.NVariables
                    if isempty(a.VarMetaInfo) || isempty(a.VarMetaInfo(i, 'VariableDescription'))
                        varDesc = 'NA';
                    else
                        varDesc = a.VarMetaInfo{i, 'VariableDescription'};
                    end
                    fprintf('%s%s: %s\n',colPad, varNames{i}, varDesc)
                end
            else
                disp(a.VarMetaInfo)
            end
        end % disp method
    end % hidden methods block
end % MetaData class

%== Helper functions
function obj = createMetadataSet(obj, varargin)
% Create a MetaData object with input data
argCount = 1;
nInput = nargin -1;

%== The first argument must be a dataset or a MetaData
arg = varargin{argCount};
if isa(arg, 'bioma.data.MetaData')
    obj = arg;
elseif isa(arg, 'dataset')
    obj.Data = arg;
    [obj.NSamples, obj.NVariables] = size(obj.Data);
else
    error(message('bioinfo:MetaData:MetaData:InvalidInputDataType'))
end

%== Process individual input arguments
s_argCount = argCount;
while argCount < nInput
    argCount = argCount + 1;
    arg = varargin{argCount};
    %== Guess if the input is param name/value pairs
    if bioma.util.isString(arg)
        %== Start of param name/value pairs
        s_argCount = argCount;
        argCount = argCount - 1;
        break;
    elseif isa(arg, 'dataset') && argCount == 2
        %== Add VarMetaInfo
        % Check if there is a VariableDescription column name
        obj.VarMetaInfo = arg;
        %== Validate row numbers in VarMetaInfo equals col number in Data.
        validateVMIRowNums(obj);
    elseif iscellstr(arg) && argCount == 2
        if length(arg) ~= obj.NVariables
            error(message('bioinfo:MetaData:MetaData:NumOfVarDescriptionsNotMatch'))
        end
        obj.VarMetaInfo = dataset(arg(:), 'VarNames', {'VariableDescription'},...
                                          'ObsNames', get(obj.Data, 'VarNames'));
    else
        if argCount == 2
            error(message('bioinfo:MetaData:MetaData:InvalidVarMetaInfo'))
        end
    end
end % while argCount < nargin processing individual input argument

%==Processing param name/value pair
if argCount < nInput
    obj = parseInputs(obj, varargin{s_argCount:end});
end

if isempty(obj.VarMetaInfo)
    obj.VarMetaInfo = getDefaultVarMetadata(obj, []);
end
end % end of createMetadataSet

function obj = createMetadataSetFromFile(obj, varargin)
% Create MetaData object by reading data from a file
argCount = 1;

%== The first argument must be a string "file"
arg = varargin{argCount};
try
    validatestring(arg, {'File', 'file'});
catch ME 
    error(message('bioinfo:MetaData:MetaData:InvalidArgumentFile'))
end

argCount = argCount + 1;
%== Get the file names
filenameArg = varargin{argCount};
[inStruct, inToPass] = parseFileInputs(varargin{argCount+1:end});
[data, varMetaInfo] = readMixedData(filenameArg, inStruct);

%== Create the data matrix
obj = createMetadataSet(obj, data, varMetaInfo, inToPass{:});

if isempty(obj.VarMetaInfo)
    obj.VarMetaInfo = getDefaultVarMetadata(obj, []);
end
end % end of createMetadataSetFromFile

function [data, midata] = readMixedData(filename, inPV)
%READMIXEDDATA returns dataset objects by reading a file.
% 
%   [DATA, VARMETAINFO]=READMIXEDDATA(FILENAME) creates a dataset object
%   DATA from a text file containing a table with mixed numeric and text
%   columns and a dataset object VARMETAINFO if there are description
%   information about the variables at the beginning of the file.
%
%   [...] = READMIXEDDATA(..., 'PATH', FILEPATH) specifies the path to
%   FILENAME.
%
%   [...] = READMIXEDDATA(..., 'DELIMITER', DELIMITER) specifies the column
%   separator. It must be a string. Default is '\t' for tab.
%
%   [...] = READMIXEDDATA(..., 'ROWNAMES', ROWNAMES). ROWNAMES can be a
%   cell array of strings giving the row names, or a single number giving
%   the column of the table which contains the row names, or a character
%   string giving the name of the table column containing the row names. By
%   default the first column in the table is used for the row names. If
%   ROWNAMES is set to empty, the row will be numbered.
%
%   [...] = READMIXEDDATA(..., 'COLNAMES', COLNAMES) only reads the data
%   from columns with names in cell array COLNAMES. The default is to read
%   data from all columns in the table. Set COLNAMES to empty to specify
%   column names are missing from the table, the column will be numbered.
%
%   [...] = READMIXEDDATA(..., 'VARDESCCHAR', VC) specifies lines beginning
%   with character VC as variable name descriptions used for creating the
%   VarMetaInfo dataset. By default no variable meta information will be
%   read.

%   [...] = READMIXEDDATA(..., 'MAXHEADERLINES', MH) specifies maximum
%   number of header lines in the file. The default is 100.

%== Get full file name including path
full_filename = [inPV.Path, filesep, filename];

%== Try open the file
try
    fopenMessage = '';
    [fid, fopenMessage] = fopen(full_filename,'rt'); %#ok
catch ME 
    error(message('bioinfo:MetaData:MetaData:CanNotOpenFile', full_filename, fopenMessage));
end

%== Handle the headers and variable meta information lines
numHeaderLines = 1;
varDescLines = cell(inPV.MaxHeaderLines, 1);
numVarDescLines = 1;
colNames = [];
allData = [];

try
    %== Figure out header lines
    while numHeaderLines <= inPV.MaxHeaderLines || feof(fid)
        fpos = ftell(fid);
        theLine = fgetl(fid);
        
        % If reach the end-of-file, get out
        if feof(fid)
            break;
        end
        
        if ~isempty(theLine)
            theLine = strtrim(theLine);
            colNames = strread(theLine, '%s', 'delimiter', inPV.Delimiter);
            % if the size is greater than one it could be a column name lines
            if numel(colNames) > 1
                if ~isempty(inPV.VarDescChar)
                    varDescLines(numVarDescLines:end) = [];
                end
                
                % Rewind back one line if there is row for column names
                if isempty(inPV.ColNamesToRead)
                    fseek(fid, fpos, -1);
                end
                break
            else
                numHeaderLines = numHeaderLines + 1;
            end
            %== Get the variable name descriptions if available
            if ~isempty(inPV.VarDescChar) && strncmp(theLine, inPV.VarDescChar, 1)
                varDescLines{numVarDescLines} = theLine;
                numVarDescLines = numVarDescLines+1;
            end
        end
    end % end of while
    
    %== Error if there is table or too many header lines
    if numHeaderLines >= inPV.MaxHeaderLines || isempty(colNames)
        fclose(fid);
        error(message('bioinfo:MetaData:MetaData:NoColumnLikeTable', filename, inPV.MaxHeaderLines));
    end
    
    %== Figure out the number of variables
    numColumns = numel(colNames);
    colToRead = true(1, numColumns);
    isNumericCol = false(1,numColumns);
    
    %== Decide the columns to read
    if ~isempty(inPV.ColNamesToRead) && iscellstr(inPV.ColNamesToRead)
        matches = cell2mat(cellfun(@(x)(strncmpi(x, colNames, length(x))),...
            inPV.ColNamesToRead, 'UniformOutput', false));
        % Check that we found all the columns
        foundCols = any(matches);
        if any(~foundCols)
            notFound = find(~foundCols);
            for count = 1:numel(notFound)
                warning(message('bioinfo:MetaData:MetaData:ColumnNotFound', inPV.ColNamesToRead{ notFound( count ) }));
            end
        end
        colToRead = any(matches,2);
        if isnumeric(inPV.RowNames) && inPV.RowNames ~= 0
            colToRead(inPV.RowNames) = true;
        end
    end
    colNames = colNames(colToRead);
    
    %== Read one line and use this to identify text columns
    fpos = ftell(fid);
    line1 = strread(fgetl(fid),'%s','delimiter', inPV.Delimiter);
    % Rewind back one line
    fseek(fid, fpos, -1);
    % See if the columns contain numeric data
    for count = 1:numel(line1)
        if colToRead(count) && (~isnan(str2double(line1{count})) ||...
                strcmp(line1{count},'NaN'))
            isNumericCol(count) = true;
        end
    end
    
    %==Read data from file using textscan
    allData = bioinfoprivate.bioReadMixedData(fid,isNumericCol,inPV.Delimiter,~colToRead);
    %== Close file
    fclose(fid);
catch ME
    fclose(fid);
    bioinfoprivate.bioclsrethrow('MetaData','MetaData', ME)
end

%== Create dataset object of the data
% Figure out which column contains row names
numColsRead = numel(allData);
numRowsRead = numel(allData{1});
rowNames = {};

if isnumeric(inPV.RowNames)
    rowNamesColIdx = inPV.RowNames;
elseif ischar(inPV.RowNames)
    rowNamesColIdx = strmatch(inPV.RowNames, colNames);
    
    if isempty(rowNamesColIdx)
        rowNamesColIdx = 0;
        warning(message('bioinfo:MetaData:MetaData:RowNameColumnNotFound', inPV.RowNames));
    end
    
elseif iscellstr(inPV.RowNames)
    if numel(inPV.RowNames) == numRowsRead
        rowNames = inPV.RowNames;
        rowNamesColIdx = 0;
    else
        error(message('bioinfo:MetaData:MetaData:NumberRowNamesNotMatchRowNumber', numel( inPV.RowNames ), numRowsRead));
    end
end

dataColIdx = find((1:numColsRead) ~= rowNamesColIdx);
numDataCols = numel(dataColIdx);
if isempty(rowNames) && rowNamesColIdx ~= 0
    rowNames = allData{rowNamesColIdx};
end
%== Convert cell array allData into a dataset
dataColNames = colNames(dataColIdx);
data = dataset({cell(numRowsRead, numDataCols), dataColNames{:}});%#ok
for loop = 1:numDataCols
    data = replacedata(data, allData{dataColIdx(loop)}, dataColNames{loop});
end

if ~isempty(rowNames)
    if (numel(unique(rowNames))~=numel(rowNames))
        error(message('bioinfo:MetaData:MetaData:DuplicatedRowNames'))
    end
    data = set(data, 'ObsNames', rowNames);
end

%== Process variable meta information
if ~all(arrayfun(@(x)isempty(x{:}), varDescLines))
    varLabelDescriptions = cell(numDataCols,1);
    % Find the end index of the variable name label before the descriptions
    endIndices = cellfun(@(x)regexp(varDescLines,...
        ['[' inPV.VarDescChar ' ]\s?' x '[:]?\s?'], 'end'), colNames, 'UniformOutput', false);
    
    for loop=1:numDataCols
        endIdx = find(~cellfun(@(x)isempty(x), endIndices{dataColIdx(loop)}));
        miLineStr = varDescLines{endIdx};
        varLabelDescriptions{loop} = miLineStr(endIndices{dataColIdx(loop)}{endIdx}:end);
    end
    
    midata = dataset(varLabelDescriptions, 'VarNames', 'VariableDescription',...
                                           'ObsNames', (get(data, 'VarNames'))');
else
    midata = dataset;
end

end %end of readMixedData

function [inStruct, inToPass] = parseFileInputs(varargin)
% Parse input PV pairs for read data files.

% Check for the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:MetaData:MetaData:IncorrectNumberOfArguments'))
end

% Allowed inputs
okargs = {'path', 'delimiter', 'rownames', 'colnames',...
    'vardescchar','maxheaderlines',...
    'name', 'description', 'samplenames', 'varnames','vardescription'};

% Defaults
inStruct.Path = '';            % Path to file
inStruct.Delimiter = '\t';      % Delimiter
inStruct.RowNames = 1;          % The column number of the column contains row names
inStruct.ColNamesToRead =1;     % Cell array of colnames to be read
inStruct.VarDescChar = '';      % Variable name description lines starting character
inStruct.MaxHeaderLines = 100;   % The maximum number of header lines
% Defaults to pass to the object
inToPass = varargin(:);
passIdx = true(1,numel(varargin));
for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1 % path
            inStruct.Path = pval;
            passIdx(j:j+1) = false;
        case 2 % delimiter
            if ischar(pval)
                inStruct.Delimiter = pval;
            else
                error(message('bioinfo:MetaData:MetaData:InvalidDelimiter'))
            end
            passIdx(j:j+1) = false;
        case 3 % Row names
            if isempty(pval)
                inStruct.RowNames = 0;
            elseif isnumeric(pval) && isscalar(pval)
                inStruct.RowNames = pval;
            elseif ischar(pval)
                inStruct.RowNames = strtrim(pval);
            elseif iscellstr(pval)
                inStruct.RowNames = pval(:);
                
                if numel(pval) == 1
                    inStruct.RowNames = pval{:};
                end
            else
                error(message('bioinfo:MetaData:MetaData:InvalidRowNamesInput'))
            end
            passIdx(j:j+1) = false;
        case 4 % column names - column names to be read
            if isempty(pval) || iscellstr(pval)
                inStruct.ColNamesToRead = pval;
            elseif ischar(pval)
                inStruct.ColNamesToRead = {pval};
            else
                error(message('bioinfo:MetaData:MetaData:ColNamesMustBeCellString'));
            end
            passIdx(j:j+1) = false;
        case 5 % variable meta information character
            if ischar(pval)
                inStruct.VarDescChar = pval;
            else
                error(message('bioinfo:MetaData:MetaData:InvalidVarDescChar'))
            end
            passIdx(j:j+1) = false;
        case 6 % Max number of header lines
            if isnumeric(pval) && isscalar(pval)
                inStruct.MaxHeaderLines = pval;
            else
                error(message('bioinfo:MetaData:MetaData:InvalidMaxHeaderLines'))
            end
            passIdx(j:j+1) = false;
    end
end
inToPass = inToPass(passIdx);
end %parseFileInputs

%----------
function validateVMIRowNums(obj)
% Validate the number of rows in VarMetaInfo equals the number of columns
% in Data.
if ~isempty(obj.VarMetaInfo) && size(obj.VarMetaInfo,1) ~= size(obj.Data, 2)
    error(message('bioinfo:MetaData:MetaData:NumOfVariablesNotMatch'))
end
end

%------------
function data = horzCombineDataset(a, b)
% Horizontal combine datasets A and B. The variables in DA and DB must be
% unique. The function pads the missing observation names in A or B and
% returns a dataset C with variables in A and B with all the unique
% observation names.

%== Transform dataset A with missing observations from B
data_a = [a; getPaddingRows(a, b)];
%== Transform dataset B with missing observations from A
data_b = [getPaddingRows(b, a); b];
%== Horizontal concatenate padded A and B
data = [data_a, data_b];
end

%------------
function data = vertCombineDataset(a, b)
% Vertical combine datasets A and B. The variables in DA and DB must be
% unique. The function pads the missing observation names in A or B and
% returns a dataset C with variables in A and B with all the unique
% observation names.

a_varNames = get(a, 'VarNames');
b_varNames = get(b, 'VarNames');

[~, idxA, idxB] = intersect(a_varNames, b_varNames);
d_comm = [a(:, idxA); b(:, idxA)];
idxDA = 1:size(a,2) ~= idxA;
idxDB = 1:size(b,2) ~= idxB;

%== Horizontal concatenate A and B
if idxDA == 0 && idxDB == 0
    data = d_comm;
elseif idxDA == 0 && idxDB ~= 0
    d_diff = horzCombineDataset(a(:, idxA), b(:, idxDB));
    data = [d_comm, d_diff(:, 2:end)];
elseif idxDA ~=0 && idxDB == 0
    d_diff = horzCombineDataset(a(:, idxDA), b(:, idxB));
    data = [d_comm, d_diff(:, 2:end-1)];
end
end

%-------------
function data = getPaddingRows(a, b)
%Create a dataset data for padding dataset A with unique observations
% (rows) from dataset B. If the column is numeric or logical NaNs are
% added. If the column is a cellstr 'NA' is added. Columns with other data
% types will be padded with empty cell.

%== Obervation (row) names from A and B.
a_ObsNames = get(a, 'ObsNames');
b_ObsNames = get(b, 'ObsNames');
%== VarNames from A
a_VarNames = get(a, 'VarNames');

%== Find the common observation names
comm_BA = ismember(b_ObsNames, a_ObsNames);
%== Number of unique observations in B not in A
num_UB = sum(~comm_BA);
if num_UB == 0
    data = dataset;
else
    data = dataset({cell(num_UB, size(a, 2)), a_VarNames{:}},...
                    'ObsNames', b_ObsNames(~comm_BA));%#ok
    for loop = 1:size(a,2)
        if isnumeric(a.(a_VarNames{loop})) || islogical(a.(a_VarNames{loop}))
            data = replacedata(data, nan(num_UB, 1), a_VarNames{loop});
        elseif ischar(a.(a_VarNames{loop})) || iscellstr(a.(a_VarNames{loop}))
            data = replacedata(data, cellstr(repmat('NA', num_UB, 1)), a_VarNames{loop});
        end
    end
end
end

function dcell = datasetToCell(ds)
[nRows, nCols] = size(ds);
dcell = cell(nRows, nCols);
for irow = 1:nRows
    for jcol = 1:nCols
        dcell{irow, jcol} = ds{irow, jcol};
    end
end
end

function obj = setData(obj, data)
if isa(data, 'dataset')
    obj.Data = data;   
    [obj.NSamples, obj.NVariables] = size(obj.Data);
    obj.VarMetaInfo = getDefaultVarMetadata(obj);
else
    error(message('bioinfo:MetaData:MetaData:InvalidDataInputFormat'))
end
end

function obj = setVarMetaInfo(obj, vmi)
if isempty(vmi)
     obj.VarMetaInfo = getDefaultVarMetadata(obj, vmi);  
elseif isa(vmi, 'dataset')
    [descLabel, labelIdx] = validateVarDescCol(vmi);
    if  obj.NVariables == size(vmi, 1)
        if ~isempty(descLabel) && labelIdx ~= 1
            vmi.VariableDescription = vmi.(descLabel);
            vmi.(descLabel) = [];
        end
        obj.VarMetaInfo = getDefaultVarMetadata(obj, vmi);
    else
        error(message('bioinfo:MetaData:MetaData:MismatchVarNumber'))
    end
    
elseif iscellstr(vmi) && obj.NVariables == numel(vmi)
    obj.VarMetaInfo = dataset((vmi(:)), 'VarNames', {'VariableDescription'},...
        'ObsNames', get(obj.Data, 'VarNames'));
else
    error(message('bioinfo:MetaData:MetaData:InvalidVarMetaInfoInput'))
end
end

function obj = setVariableNames(obj, names)
    try
        obj.Data = bioma.util.nameList(obj.Data, 'VarNames', [], names, false, 'Var');
        if ~isempty(obj.VarMetaInfo)
            obj.VarMetaInfo = set(obj.VarMetaInfo, 'ObsNames', get(obj.Data, 'VarNames'));
        end
    catch ME
        bioinfoprivate.bioclsrethrow('MetaData','MetaData', ME);
    end
end

function obj = setSampleNames(obj, names)
try
    if ischar(names) && size(names, 1) == 1 && obj.NSamples == 1
        sNames = {names};
    else
        sNames = bioma.util.setDimNames(names,...
            obj.NSamples, 1, 'Sample', false);
    end
catch ME
    bioinfoprivate.bioclsrethrow('MetaData','MetaData', ME);
end

obj.Data = set(obj.Data, 'ObsNames', sNames);
end

function obj = setName(obj, x)
if isempty(x)
    return;
end

if ischar(x)
    obj.Name = x;
elseif iscellstr(x)
    obj.Name = x{1};
else
    error(message('bioinfo:MetaData:MetaData:BadNameInputFormat'))
end
end % setName

function obj = setDescription(obj, desc)
if isempty(desc)
    return;
elseif ischar(desc)
    obj.Description = desc;
elseif iscellstr(desc)
    obj.Description = desc{1};
else
    error(message('bioinfo:MetaData:MetaData:BadDescInputFormat'))
end
end %setDescription

function obj = setDimensionLabels(obj, labels)
if iscellstr(labels) && numel(labels) >= 2
    if isempty(labels{1})
        obj.DimensionLabels{2} = labels{2};
    else
        obj.DimensionLabels = labels(1:2);
    end
elseif ischar(labels)
    obj.DimensionLabels{1} = labels;
end
end

%--------
function varMI = getDefaultVarMetadata(obj, oldVMI)
% If VarMetaInfo is empty, create a dataset with VariableDescription set
% to 'NA'. Check the obj's DATA variable names with the VarMetaInfo's row
% names. If the names matches, not changes. If any names in Data not in
% VarMetaInfo, the missing names are added to VarMetaInfo.

if nargin == 1 || (nargin > 1 && isempty(oldVMI)) 
    oldVMI = obj.VarMetaInfo;
    errFlag = false;
end

if nargin == 2 && ~isempty(oldVMI)
    errFlag = true;
end

varDesc = cellstr(repmat('NA', obj.NVariables, 1));

if isempty(oldVMI) 
    varMI = dataset(varDesc,  'VarNames', 'VariableDescription',...
                              'ObsNames', obj.variableNames');
else    
   varDNames = get(obj.Data, 'VarNames');
   varMINames = get(oldVMI, 'ObsNames');
   varColNames = get(oldVMI, 'VarNames');
   
   [itsVarNames, vdIdx, vMIIdx] = intersect(varDNames, varMINames);
   
   if errFlag
       if isempty(itsVarNames)
          error(message('bioinfo:MetaData:MetaData:NoMatchingVarNames'))
       elseif numel(vMIIdx) < numel(varMINames)
           warnNames = varMINames(~ismember(1:numel(varMINames), vMIIdx));
           warnStr = '';
           for i=1:numel(warnNames)
               warnStr = [mf warnNames{i} ', '];
           end
           warnStr = warnStr(1:end-2);
           warning(message('bioinfo:MetaData:MetaData:IgnoredNotMatchingVarNames', warnStr))
       end
       
   end
   varMI = dataset({cell(obj.NVariables, size(oldVMI, 2)), varColNames{:}}, ...
                        'ObsNames', (get(obj.Data, 'VarNames'))'); %#ok
   varMI  = replacedata(varMI , varDesc, 'VariableDescription');
   
   if ~isempty(itsVarNames)
       for loop=1:size(varMI, 2)
           varMI (vdIdx, :) = replacedata(varMI(vdIdx, :), ...
               oldVMI.(varColNames{loop})(vMIIdx), varColNames{loop});
       end
   end
end
end

function [vrel, vidx] = validateVarDescCol(a)
% Validate the column names in dataset a contains 'VariableDescription',
% 'VarDescription' or 'Description'. Returns the validate result as true or
% false in vrel.
vrel = [];
vidx = 1;
varNames = get(a, 'VarNames');
namesToCheck ={'VariableDescription',...
               'VarDescription',...
               'Description',...
               'LabelDescription'};

for loop = 1: numel(namesToCheck)
    try 
        validatestring(namesToCheck{loop}, varNames);
        vrel = namesToCheck{loop};
        vidx = loop;
        break;
    catch ME %#ok
        %== do nothing
    end
end    
end

%-------------------
function obj = parseInputs(obj, varargin)
%PARSEINPUTS Parse input PV pairs for classes. 
if nargin < 2
    return;
end

% Check for the right number of inputs
if rem(nargin-1, 2)== 1
    error(message('bioinfo:MetaData:MetaData:IncorrectNumberOfArguments'))
end

okargs = [properties(obj);'VariableNames';'SampleNames'];
for i = 1:(nargin-1)/2
    name = varargin{2*i-1};
    value = varargin{2*i};
    k = bioinfoprivate.pvpair(name,[],okargs,'MetaData:MetaData');
    name = okargs{k};
    
    switch name
        case 'Name' % Names
            obj = setName(obj,value);
        case 'Description' % Description
            obj = setDescription(obj,value);
        case 'SampleNames' % SampleNames
            obj = setSampleNames(obj,value);
        case 'VariableNames' % VariableNames
            obj = setVariableNames(obj,value);
        case 'DimensionLabels' % DimemsionLabels
            obj = setDimensionLabels(obj,value);
        case 'NSamples'
            error(message('bioinfo:MetaData:MetaData:IllegalPropertySet', name));
        case 'NVariables'
            error(message('bioinfo:MetaData:MetaData:IllegalPropertySet', name));
        otherwise
            error(message('bioinfo:MetaData:MetaData:UnknownPropName', name));
    end
end
end


