classdef ExptData
    %EXPTDATA Class to contain microarray experiment data values.
    %
    %   The ExptData class contains one or more identically sized matrices
    %   of the same dimension containing measurement values or other data.
    %   All the matrices (elements) in the class are structured so that
    %   rows represent 'features' on the array and columns represent
    %   samples used in the experiment. Validity methods enforce that row
    %   and column names of for the matrices are identical.
    %
    %   ExptData properties:
    %       Name            - Name of the ExptData object.
    %       NFeatures       - Number of features.
    %       NSamples        - Number of samples.
    %       NElements       - Number of elements (DataMatrix objects).
    %       ElementClass    - Class names of the elements.
    %
    %   ExptData methods:
    %       ExptData        - Create an ExptData object.
    %       elementData     - Access and set a DataMatrix element.
    %       elementNames    - Access and set the DataMatrix element names.
    %       featureNames    - Access and set the feature names.
    %       sampleNames     - Access and set the sample names.
    %       dmNames         - Access and set the Name properties of DataMatrix elements.
    %       size            - Size of an ExptData object.
    %       combine         - Combine data from two ExptData objects.
    %       isempty         - True for empty ExptData object.
    %
    %   Examples:
    %     d1 = bioma.data.DataMatrix(rand(5, 4),'Name', 'ExptDataExample',...
    %                                'RowNames','Feature', 'ColNames', 'Sample');
    %     v1 = [11 12 13 14; 21 22 23 24; 31 32 33 34; 41 42 43 44;51 52 53 54]; 
    %     ed = bioma.data.ExptData(d1, v1, 'ElementNames', {'RedData', 'GreenData'})
    %
    %   See also DATAMATRIX, DATASET, EXPRESSIONSET, METADATA, MIAME.
    
    %   Copyright 2009-2012 The MathWorks, Inc.

    
    properties (SetAccess = 'public', GetAccess = 'public')
        %NAME Name of the ExptData object.
        %    The Name property is a string specifying the name of an
        %    ExptData object. 
        %
        %    See also EXPTDATA.
        Name = '';
    end
    
    properties (SetAccess = 'private', GetAccess = 'public')
        %NFEATURES The number of features.
        %   The NFeatures property is a scalar equal to the number of
        %   features (rows).
        %
        %   See also EXPTDATA.
        NFeatures = 0;
        
        %NSAMPLES The number of samples.
        %   The NSamples property is a scalar equal to the number of
        %   samples (columns).
        %
        %   See also EXPTDATA.
        NSamples = 0;
        
        %NELEMENTS The number of DataMatrix objects in ExptData.
        %   The NElements property is a scalar equal to the number of
        %   DataMatrix objects or elements in the ExptData object.
        %
        %   See also EXPTDATA.
        NElements = 0;
        
        %ELEMENTCLASS Class names of the elements.
        %   The ElementClass property is a cell array of class names of the
        %   elements (DataMatrix objects) in the ExptData object.
        %
        %   See also EXPTDATA.
        ElementClass = [];
    end
    
    properties(SetAccess = 'private', GetAccess = 'private', Hidden=true)
        %MATRICES Cell array containing the data arrays
        % The Matrices property is a cell array with equal sized MATLAB
        % numeric or logical arrays as elements.
        Matrices = [];

        %SAMPLENAMES Names of the samples.
        %   The SampleNames property is a cell array of strings containing
        %   the sample names. It is also the column names of the ExptData.
        SampleNames =[];
        
        %FEATURENAMES Names of the features.
        %   The FeatureNames property is a cell array of strings containing
        %   the feature names. It represents the row names of the ExptData.
        FeatureNames =[];
  
        %ELEMENTNAMES Names of the DataMatrix elements.
        %   The ElementNames property is a cell array of strings containing
        %   the DataMatrix names.
        ElementNames =[];
        
        %DMNAMES Names of the name property of DataMatrix objects.
        %   The DMNames property is a cell array of strings of contains
        %   from the Name property of the DataMatrix object.
        DMNames = [];
    end
    
    methods
        %== Constructor
        function obj = ExptData(varargin)
            %EXPTDATA Create an ExptData object.
            %
            %   E = EXPTDATA(DATA1, DATA2, ...) creates an ExptData object
            %   E from variables DATA1, DATA2, ..., . Each variable must be
            %   a DataMatrix object, a numeric matrix or a logical matrix.
            %   All variables must be the same size. The DataMatrix objects
            %   must have identical row names and column names. Any numeric
            %   matrices or logical matrices are converted to a DataMatrix
            %   object with row names and column names.
            %
            %   E = EXPTDATA(..., {DATA, 'NAME'}, ...) creates an ExptData
            %   object E with a DataMatrix element named 'NAME' in E. The
            %   element names must be a valid MATLAB identifiers and
            %   unique.
            %
            %   E = EXPTDATA({DATA1, DATA2,...}, ...) creates an ExptData
            %   object E from a cell array of DataMatrix objects, numeric
            %   matrices or logical matrices. The elements in the cell
            %   array must be the same size. The DataMatrix objects must
            %   have identical row names and column names. Any numeric
            %   matrices or logical matrices are converted to a DataMatrix
            %   object with row names and column names.
            %
            %   E = EXPTDATA(..., 'ELEMENTNAMES', ENAMES) creates
            %   DataMatrix elements with specified names. The ENAMES must
            %   be valid MATLAB identifiers, and unique. If names already
            %   exist for the DataMatrix elements, they are replaced by
            %   ENAMES.
            %
            %   E = EXPTDATA(..., 'FEATURENAMES',FNAMES) creates an
            %   ExptData object E, with feature (row) names specified by
            %   FNAMES. The number of the names in FNAMES must match the
            %   number of rows in the input matrices. The feature names can
            %   be a cell array of strings, a character array, or a numeric
            %   or logical vector. FNAMES can also be a single string as a
            %   prefix for feature names; feature numbers will be appended
            %   to the prefix. FNAMES can also be a logical true or false
            %   (default); if true, default unique names will be assigned
            %   to the features.
            %
            %   E = EXPTDATA(..., 'SAMPLENAMES',SNAMES) creates an ExptData
            %   object E, with sample (column) names specified by SNAMES.
            %   The number of names in SNAMES must match the number of
            %   columns in the input matrices. The sample names can be a
            %   cell array of strings, a character array, or a numeric or
            %   logical vector. SNAMES can also be a single string as a
            %   prefix for sample names; sample numbers will be appended to
            %   the prefix. SNAMES can also be a logical true or false
            %   (default); if true, default unique names will be assigned
            %   to the samples.
            %
            %   See also EXPTDATA.
            
            if nargin == 0
                return;
            end
            obj = createExptData(obj, varargin{:});
        end %end of ExptData constructor
    end
    
    methods
        function data = elementData(obj, idx, dm)
            %ELEMENTDATA Access and set a DataMatrix element of an ExptData object.
            %
            %   D = ELEMENTDATA(E, ELMT) returns a DataMatrix object D of
            %   specified element ELMT from ExptData object E. ELMT can be
            %   a positive integer, or a string specifying an element name.
            %
            %   B = ELEMENTDATA(E, ELMT, DM) returns an ExptData object B
            %   with specified element replaced by a DataMatrix DM. The
            %   column names of DM must matching the existing sample names
            %   and the row names must matching existing feature names in
            %   A.
            %
            %   See also EXPTDATA, DMNAMES, ELEMENTNAMES, FEATURENAMES, SAMPLENAMES.
            
            %== Input check
            bioinfochecknargin(nargin, 2, 'ExptData:elementData')
            
            try
                if nargin ==2
                    data = getElementDM(obj, idx, 'elementData');
                else
                    data = setElementDM(obj, idx, dm, 'elementData');
                end
            catch ME
                bioinfoprivate.bioclsrethrow('ExptData','elementData', ME)
            end
        end
        
        function out = elementNames(obj, idx, names)
            %ELEMENTNAMES Access and set DataMatrix element names.
            %
            %   ENAMES = ELEMENTNAMES(E) returns a cell array of all
            %   DataMatrix element names of ExptData object E.
            %
            %   ENAMES = ELEMENTNAMES(E, I) returns the names of elements
            %   specified by I in ExptData object E. I can be a positive
            %   integer, a vector of positive integers, a string specifying
            %   an element name, a cell array containing one or more
            %   element names, or a logical vector.
            %
            %   B = ELEMENTNAMES(E, I, NAMES) returns an ExptData object B
            %   with the names of specified elements set to NAMES. The
            %   number of names in NAMES must equal the number of elements
            %   specified by I. NAMES can be a cell array of strings, a
            %   character array, or a numeric vector. NAMES can also be a
            %   single string as a prefix for element names; element
            %   numbers will be appended to the prefix. NAMES can also be a
            %   logical true or false (default); if true, default unique
            %   names will be assigned to the elements.
            %
            %   See also EXPTDATA, DMNAMES, FEATURENAMES, SAMPLENAMES.
            
            %== Input check
            bioinfochecknargin(nargin, 1, 'ExptData:elementNames')
            try
                if nargin < 2
                    out = obj.ElementNames;
                elseif nargin < 3
                    out = bioma.util.nameList(obj, 'elementNames',idx);
                elseif nargin ==3
                    if isempty(idx) || bioma.util.isColon(idx)
                        out = setElementNames(obj, names);
                    else
                        out = updateNameProperties(obj,...
                            'elementNames', idx, names);
                    end
                end
            catch ME
                bioinfoprivate.bioclsrethrow('ExptData','elementNames', ME)
            end
        end % elementNames
        
        function out = featureNames(obj, idx, names)
            %FEATURENAMES Access and set the feature names.
            %
            %   FNAMES = FEATURENAMES(E) returns a cell array of feature
            %   (row) names of ExptData object E.
            %
            %   FNAMES = FEATURENAMES(E, I) returns the names of features
            %   specified by I in ExptData object E. I can be a positive
            %   integer, a vector of positive integers, a string specifying
            %   a feature name, a cell array containing one or more feature
            %   names, or a logical vector.
            %
            %   B = FEATURENAMES(E, I, NAMES) returns an ExptData object B
            %   with the names of specified features set to NAMES. The
            %   number of names in NAMES must equal the number of features
            %   specified by I. NAMES can be a cell array of strings, a
            %   character array, or a numeric vector. NAMES can also be a
            %   single string as a prefix for feature names; feature
            %   numbers will be appended to the prefix. NAMES can also be a
            %   logical true or false (default); if true, default unique
            %   names will be assigned to the features.
            %
            %   See also EXPTDATA, DMNAMES, ELEMENTNAMES, SAMPLENAMES.
            
            %== Input check
            bioinfochecknargin(nargin, 1, 'ExptData:featureNames')
            try
                if nargin < 2
                    out = obj.FeatureNames;
                elseif nargin < 3
                    out = bioma.util.nameList(obj, 'featureNames',idx);
                elseif nargin ==3
                    if isempty(idx) || bioma.util.isColon(idx)
                        out = setFeatureNames(obj, names);
                    else
                        out = updateNameProperties(obj,...
                            'featureNames', idx, names);
                    end
                end
            catch ME
                bioinfoprivate.bioclsrethrow('ExptData','featureNames', ME)
            end
        end % featureNames
        
        function out = sampleNames(obj, idx, names)
            %SAMPLENAMES Access and set the sample names.
            %
            %   SNAMES = SAMPLENAMES(E) returns a cell array of sample
            %   (column) names of an ExptData object E.
            %
            %   SNAMES = SAMPLENAMES(E, I) returns the names of samples
            %   specified by I in ExptData object E. I can be a positive
            %   integer, a vector of positive integers, a string specifying
            %   a sample name, a cell array containing one or more sample
            %   names, or a logical vector.
            %
            %   B = SAMPLENAMES(E, I, NAMES) returns an ExptData object B
            %   with the names of specified samples set to NAMES. The
            %   number of names in NAMES must equal the number of samples
            %   specified by I. NAMES can be a cell array of strings, a
            %   character array, or a numeric vector. NAMES can also be a
            %   single string as a prefix for sample names; sample numbers
            %   will be appended to the prefix. NAMES can also be a logical
            %   true or false (default); if true, default unique names will
            %   be assigned to the samples.
            %
            %   See also EXPTDATA, DMNAMES, ELEMENTNAMES, FEATURENAMES.
            
            %== Input check
            bioinfochecknargin(nargin, 1, 'ExptData:sampleNames')
            try
                if nargin < 2
                    out = obj.SampleNames;
                elseif nargin < 3
                    out = bioma.util.nameList(obj, 'sampleNames',idx);
                elseif nargin ==3
                    if isempty(idx) || bioma.util.isColon(idx)
                        out = setSampleNames(obj, names);
                    else
                        out = updateNameProperties(obj,...
                            'sampleNames', idx, names);
                    end
                end
            catch ME
                bioinfoprivate.bioclsrethrow('ExptData','sampleNames', ME)
            end
        end % sampleNames
        
        function out = dmNames(obj, idx, names)
            %DMNAMES Access and set the Name properties of DataMatrix elements.
            %
            %   DNAMES = DMNAMES(E) returns a cell array of names from the
            %   Name properties of the DataMatrix elements in ExptData
            %   object E.
            %
            %   DNAMES = DMNAMES(E, I) returns the names of DataMatrix
            %   elements specified by I of an ExptData object E. I can be a
            %   positive integer, a vector of positive integers, or a
            %   logical vector.
            %
            %   B = DMNAMES(E, I, NAMES) returns an ExptData object B with
            %   the Name properties of specified DataMatrix elements set to
            %   NAMES. The number of names in NAMES must equal the number
            %   of elements specified by I. NAMES can be a cell array of
            %   strings, or a numeric vector. NAMES can also be a single
            %   string as a prefix for DataMatrix names; element numbers
            %   will be appended to the prefix. NAMES can also be a logical
            %   true or false (default); if true, default unique names will
            %   be assigned to the elements.
            %
            %   See also EXPTDATA, ELEMENTNAMES, FEATURENAMES, SAMPLENAMES.
            
            %== Input check
            bioinfochecknargin(nargin, 1, 'ExptData:dmNames')
            try
                if nargin < 2
                    out = obj.DMNames;
                elseif nargin < 3
                    out = bioma.util.nameList(obj, 'dmNames',idx);
                elseif nargin ==3
                    if isempty(idx) || bioma.util.isColon(idx)
                        out = setDMNames(obj, names);
                    else
                        out = updateNameProperties(obj,...
                            'dmNames', idx, names);
                    end
                end
            catch ME
                bioinfoprivate.bioclsrethrow('ExptData','dmNames', ME)
            end
        end % dmNames
        
        function varargout = size(a,dim)
            %SIZE Size of an ExptData object.
            %
            %   D = SIZE(E) returns a two-elment row vector D containing
            %   the number of features and samples in ExptData object E.
            %
            %   [M,N] = SIZE(E) returns the number of features and columns
            %   in ExptData object E as separate output variables.
            %
            %   M = SIZE(E,DIM) returns the length of the dimension
            %   specified by the scalar DIM. For example, SIZE(E,1) returns
            %   the number of features. If DIM > NDIMS(E), M will be 1.
            %
            %   See also ExptData.
            
            if nargin == 1
                if nargout < 2
                    varargout{:} = size(a.Matrices{1});
                elseif nargout == 2
                    [varargout{1}, varargout{2}] = size(a.Matrices{1});
                else
                    [varargout{1}, varargout{2}] = size(a.Matrices{1});
                    varargout(3:nargout) = {1};
                end
                
            else
                if nargout > 1
                    error(message('bioinfo:ExptData:size:TooManyOutputs'));
                end
                
                try
                    varargout{1} = size(a.Matrices{1}, dim);
                catch ME
                    bioinfoprivate.bioclsrethrow('ExptData','size', ME)
                end
            end
        end % size method
        
        function t = isempty(obj)
            %ISEMPTY True for empty ExptData object.
            %
            %   TF = ISEMPTY(E) returns true (1) if E is an empty ExptData
            %   object and false (0) otherwise. An empty ExptData object
            %   has zero elements.
            %
            % See also EXPTDATA.
            t = isempty(obj.Matrices);
        end
        
        function c = combine(a, b)
            %COMBINE Combine data from two ExptData objects.
            %
            %   C = COMBINE(A, B) combines data from ExptData object A with
            %   ExptData object B and returns the result in ExptData object
            %   C. Sample names and feature names in A and B must be
            %   identical. The element names in A and B must be unique.
            %
            %   See also EXPTDATA.
            bioinfochecknargin(nargin, 2, 'ExptData:combine')
            
            if ~isa(a, 'bioma.data.ExptData') || ~isa(b, 'bioma.data.ExptData')
                error(message('bioinfo:ExptData:combine:InvalidInput'))
            end
            %== Find the common variable names in A and B
            if all(ismember(a.FeatureNames, b.FeatureNames)) && ...
                    all(ismember(a.SampleNames, b.SampleNames))
                try
                    data = [a.Matrices, b.Matrices];
                    
                    if ~iscellstr(a.ElementNames)
                        a.ElementNames = {a.ElementNames};
                    end
                    
                    if ~iscellstr(b.ElementNames)
                        b.ElementNames = {b.ElementNames};
                    end
                    elementNames = [a.ElementNames, b.ElementNames];
                    c = bioma.data.ExptData(data, 'ElementNames', elementNames,...
                                                  'SampleNames',  a.SampleNames,...
                                                  'FeatureNames', a.FeatureNames);
                    c = dmNames(c, ':', [a.DMNames, b.DMNames]);
                catch ME
                    bioinfoprivate.bioclsrethrow('ExptData','combine', ME)
                end
            else
                error(message('bioinfo:ExptData:combine:MismatchNames'))
            end
        end
    end % end of method block
    
    methods
        function n = get.NElements(obj)
            n = numel(obj.Matrices);
        end
        
        function elmtTypes = get.ElementClass(obj)
            elmtTypes = cell(1, obj.NElements);
            
            for loop = 1:obj.NElements
               elmtTypes{loop} = class(obj.Matrices{loop}); 
            end 
        end
    end
    
    methods(Hidden=true)
        function disp(obj)
            fprintf('Experiment Data:\n')
            colPad = repmat(' ', 1, 2);
            fprintf('%s%d features,  %d samples\n', colPad, obj.NFeatures, obj.NSamples);
            fprintf('%s%d elements\n', colPad, obj.NElements);
            if iscellstr(obj.ElementNames)
                fprintf('%sElement names: %s\n', colPad,...
                    bioma.util.catCellStrToStr(obj.ElementNames));
            elseif ischar(obj.ElementNames)
                fprintf('%sElement names: %s\n',colPad, obj.ElementNames)
            end
        end % disp
        
        function [varargout] = subsref(a, s)
            %SUBSREF Subscripted reference for an ExptData object.
            %
            %   B = A(I) returns a DataMatrix of Ith element of the
            %   ExptData A. I can be a positive integer, a string of
            %   element name. B is a DataMatrix with feature names as its
            %   row names and samples as its column names.
            %
            %   B = A(I, J) returns an subset of ExptData with Ith feature
            %   and Jth samples of the ExptData A. I and J can be a
            %   positive integer, a string of element name.
            %
            %   Note: Cell indexing A{I} is not supported.
            %
            %   P = A.PROPERTYNAME returns an ExptData property.
            switch s(1).type
                case '()'
                    %==Only one-dimension indexing
                    if numel(s.subs) > 2
                        error(message('bioinfo:ExptData:subsref:IncorrectIndexing'))
                    elseif ~isscalar(s)
                        error(message('bioinfo:ExptData:subsref:InvalidSubsExpr'));
                    end
                    
                    if numel(s.subs) == 1
                        try
                           [varargout{1:nargout}] = getElementDM(a, s.subs{:},'subsref'); 
                        catch ME
                            if strfind(ME.identifier, 'IndexOutOfBounds')
                                error(message('bioinfo:ExptData:subsref:IndexOutOfBounds', a.NElements));
                            else
                                bioinfoprivate.bioclsrethrow('ExptData','subsref',ME)
                            end
                        end
                    elseif numel(s.subs) == 2
                        try
                            rowIndices = bioma.util.findLabelIndices(a.FeatureNames, s(1).subs{1});
                            if isempty(rowIndices)
                                error(message('bioinfo:ExptData:subsref:NotSuchFeatures'))
                            elseif any(rowIndices == 0) && ~islogical(rowIndices)
                                error(message('bioinfo:ExptData:subsref:NotSuchFeature'))
                            end
                            
                            colIndices = bioma.util.findLabelIndices(a.SampleNames, s(1).subs{2});
                            if isempty(colIndices)
                                error(message('bioinfo:ExptData:subsref:NotSuchSamples'))
                            elseif any(colIndices == 0) && ~islogical(colIndices)
                                error(message('bioinfo:ExptData:subsref:NotSuchSample'))
                            end
                            
                            tmpCell = cell(1,a.NElements);
                            
                            for iloop = 1:a.NElements
                                tmpCell{iloop} = a.Matrices{iloop}(rowIndices, colIndices);
                            end
                            [varargout{1:nargout}] = bioma.data.ExptData(tmpCell,...
                                'ElementNames', a.ElementNames,...
                                'FeatureNames', a.FeatureNames(rowIndices),...
                                'SampleNames', a.SampleNames(colIndices),...
                                'DMNames', a.DMNames);
                        catch ME
                            if strfind(ME.identifier, 'IndexOutOfBounds')
                                error(message('bioinfo:ExptData:subsref:IndexOutOfBounds2D', a.NFeaturs, a.NSamples));
                            else
                                bioinfoprivate.bioclsrethrow('ExptData','subsref',ME)
                            end
                        end
                    end   
                case '{}'
                    error(message('bioinfo:ExptData:subsref:CellSubscript'));
                case '.'
                    % Intercept only when it is a valid property name,
                    % otherwise error out.
                    queryName = s(1).subs;
                    if any(strcmp(properties(a),queryName))
                        [varargout{1:nargout}] = builtin('subsref', a, s);
                    elseif any(strcmp(methods(a),queryName))
                        [varargout{1:nargout}] = builtin('subsref', a, s);
                    else
                        error(message('bioinfo:ExptData:subsref:NotAProperty'));
                    end                    
            end
        end %subsref
        
        function obj = subsasgn(obj,s,b)
            % Subscripted assignment to an ExptData object.
            %
            %   A = SUBSASGN(A,S,B) is called for the syntax A(I)=B.
            %
            %   A(I) = B assigns the values of the DataMatrix object or a
            %   valid MATLAB numeric or logical array, B into the Ith
            %   element of the ExptData object A. I can be a positive
            %   integer, or a string of name. The size of B must match the
            %   number of features and the number of samples in ExptData A.
            %   If B is a DataMatrix, its row names and column names must
            %   matching the feature names and sample names of ExptData B.
            %
            %   Note: Cell indexing A{I}= B is not supported.
            %
            %   A.PROPERTYNAME = P assigns to an ExptData property.
            fcnName = 'subsasgn';
            switch s(1).type
                case '()'
                    if numel(s.subs) ~= 1
                        error(message('bioinfo:ExptData:subsasgn:IncorrectIndexing'))
                    elseif ~isscalar(s)
                        error(message('bioinfo:ExptData:subsasgn:InvalidSubsExpr'));
                    end
                    obj = setElementDM(obj, s.subs{:}, b, fcnName);
                case '{}'
                    error(message('bioinfo:ExptData:subsasgn:CellSubscript'));
                case '.'
                    % Intercept only when it is a valid property name,
                    % otherwise error out.
                    propName = s(1).subs;
                    if any(strcmp(properties(obj),propName))
                        switch propName
                            case 'Name' % Names is the only property can be set with dot notation
                                obj = setName(obj,b);
                            otherwise
                                error(message('bioinfo:ExptData:subsasgn:SetProhibited', propName));
                        end
                    else
                        error(message('bioinfo:ExptData:subsasgn:NotAProperty'));
                    end
            end
        end %subsasgn
    end
end % class ExptData block

%== Helper functions
function obj = createExptData(obj, varargin)
%Create an ExptData object.
argCount = 1;
nInput = nargin -1;
elementNames = repmat({''}, 1, nInput);
dmNames = repmat({''}, 1, nInput);
obj.Matrices = cell(1, nInput);

%== Processing individual input arguments
matCount = 1;

while argCount <= nInput
    arg = varargin{argCount};
    %== Guess if the input is param name/value pairs
    if bioma.util.isString(arg)
        %== Start of param name/value pairs
        break;
    elseif iscell(arg) && ~isscalar(arg) && isvector(arg)
        if numel(arg) == 2 && bioma.util.isString(arg{2})
            %=={element, name}
            obj = validateSingleArg(arg{1}, obj);
            elementNames{matCount} = arg{2};
            
            if isa(arg{1}, 'bioma.data.DataMatrix')
                obj = transferDimNames(obj, arg{1});
                obj.Matrices{matCount} = arg{1}.(':')(':');
                dmNames{matCount} = arg{1}.Name;
            else
                obj.Matrices{matCount} = arg{1};
            end
            matCount = matCount + 1;
        else
            if argCount == 1
                obj.Matrices = [];
                elementNames = [];
                dmNames = [];
            end
            dmnames = repmat({''}, 1, numel(arg));
            
            for iloop = 1:numel(arg)
                obj =  validateSingleArg(arg{iloop}, obj);
                if isa(arg{iloop}, 'bioma.data.DataMatrix')
                    obj = transferDimNames(obj, arg{1});
                    dmnames{iloop} = arg{iloop}.Name;
                    arg{iloop} = arg{iloop}.(':')(':');
                end
            end
            
            elmtnames = bioma.util.appendUniqueNumToNames('Elmt',...
                                           numel(arg)+matCount-1, matCount);
            
            obj.Matrices = [obj.Matrices arg];
            matCount = matCount + numel(arg);
            
            elementNames = [elementNames elmtnames]; %#ok
            dmNames = [dmNames dmnames]; %#ok
        end
    else
        if iscell(arg) && isscalar(arg)
            arg = arg{:};
        end
        obj = validateSingleArg(arg, obj);
        elementNames(matCount) = bioma.util.appendUniqueNumToNames('Elmt',...
                                    matCount, matCount);
        
        if isa(arg, 'bioma.data.DataMatrix')
            if isempty(obj.FeatureNames)
                obj.FeatureNames = rownames(arg);
            end
            
            if isempty(obj.SampleNames)
                obj.SampleNames = colnames(arg);
            end
            obj.Matrices{matCount} = arg.(':')(':');
            dmNames{matCount} = arg.Name;
        else
            obj.Matrices{matCount} = arg;
        end
        matCount = matCount + 1;
    end
    argCount = argCount + 1;
end % while argCount < nargin processing individual input argument

matCount = matCount - 1;
elementNames = elementNames(1:matCount);
%== fix the duplicated element names and warn about it.
[~,k1,k2] = unique(elementNames(1:matCount));
if (numel(k2)~=numel(elementNames(1:matCount)))
    idx = sort(k1(accumarray(k2(:),1)==1));
    repNames = cellstr(strcat('Elmt', num2str(idx', '%d')));
    elementNames(idx) = repNames;
    warning(message('bioinfo:ExptData:DuplicatedElmtNames'));
end

obj.Matrices = obj.Matrices(1:matCount);
obj.ElementNames = elementNames;
obj.DMNames = dmNames(1:matCount);

if isempty(obj.FeatureNames)
    obj = obj.featureNames(':', true);
end

if isempty(obj.SampleNames)
    obj = obj.sampleNames(':', true);
end

%==Processing param name/value pair
if argCount-matCount < nInput
    obj = parseInputs(obj, varargin{argCount:end});
end

end

%--------------
function obj = validateSingleArg(arg, obj)
if all([obj.NFeatures, obj.NSamples]==0)
    [obj.NFeatures, obj.NSamples] = size(arg);
elseif ~isequal([obj.NFeatures, obj.NSamples], size(arg))
    error(message('bioinfo:ExptData:ExptData:MismatchSize'));
end

if isa(arg, 'bioma.data.DataMatrix')
    if  ~isempty(obj.FeatureNames) && ~isempty(obj.SampleNames)&& ~isMetaEqual(arg, obj)
        error(message('bioinfo:ExptData:ExptData:InvalidDataMatrixInput', arg.Name));
    end
else
    try
        validateattributes(arg, {'numeric', 'logical'}, {'real', '2d'});
    catch ME 
        error(message('bioinfo:ExptData:ExptData:InvalidArgumentMatrix'))
    end
    
    if ~isMetaEqual(arg, obj)
        error(message('bioinfo:ExptData:ExptData:WrongMatrixSize'));
    end
end
end

%---------
function t = isMetaEqual(d, dm)
% TRUE if DataMatrix d equal in size, row names and column names with
% DataMatrix dm.
% TRUE if a numeric matrix is equal in size with DataMatrix dm.
if isa(d, 'bioma.data.DataMatrix')
    t = isequal(dm.NFeatures, d.NRows) && ...
        isequal(dm.NSamples, d.NCols) && ...
        isequal(dm.FeatureNames, d.RowNames) && ...
        isequal(dm.SampleNames, d.ColNames);
else
    t = isequal(dm.NFeatures, size(d, 1)) && ...
        isequal(dm.NSamples, size(d, 2));
end
end

%------------
function obj = transferDimNames(obj, dM)
% Set ExptData Feature and sample names to the row names and column names
% of the DataMatrix object dM
if isempty(obj.FeatureNames) && ~isempty(dM.RowNames)
    obj.FeatureNames = dM.RowNames;
end

if isempty(obj.SampleNames) && ~isempty(dM.ColNames)
    obj.SampleNames = dM.ColNames;
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
    error(message('bioinfo:ExptData:ExptData:IncorrectNumberOfArguments'))
end

okargs = {'Name','ElementNames','SampleNames','FeatureNames','DMNames'};
for i = 1:(nargin-1)/2
    name = varargin{2*i-1};
    value = varargin{2*i};
    k = bioinfoprivate.pvpair(name,[],okargs,'ExptData:ExptData');
    name = okargs{k};
    
    switch name
        case 'Name' % Names
            obj = setName(obj,value);
        case 'ElementNames' % Names
            obj = setElementNames(obj,value);
        case 'SampleNames' % SampleNames
            obj = setSampleNames(obj, value);
        case 'FeatureNames' % VariableNames
            obj = setFeatureNames(obj, value);
        case 'DMNames' % VariableNames
            obj = setDMNames(obj,value);
        otherwise
            error(message('bioinfo:ExptData:ExptData:SetProhibited', name));
    end
end
end

function obj = setName(obj, name)
if iscellstr(name)
    name = char(name);
end
if bioma.util.isString(name)
    obj.Name = strtrim(name(:)');
else
    error(message('bioinfo:ExptData:ExptData:InvalidInputFormat'))
end
end % Name

function obj = setElementNames(obj, names)
try    
    if ischar(names) && size(names, 1) == 1 && obj.NElements == 1
        obj.ElementNames = {names};
    else
        obj.ElementNames = bioma.util.setDimNames(names,...
            obj.NElements, 2, 'Element', false);
    end
catch ME
    bioinfoprivate.bioclsrethrow('ExptData','ExptData', ME);
end
end % ElementNames

function obj = setFeatureNames(obj, names)
try
    if ischar(names) && size(names, 1) == 1 && obj.NFeatures == 1
        obj.FeatureNames = {names};
    else
        obj.FeatureNames = bioma.util.setDimNames(names,...
            obj.NFeatures, 1, 'Feature', false);
    end
catch ME
    bioinfoprivate.bioclsrethrow('ExptData','ExptData', ME);
end
end % FeatureName

function obj = setSampleNames(obj, names)
try
    if ischar(names) && size(names, 1) == 1 && obj.NSamples == 1
        obj.SampleNames = {names};
    else
        obj.SampleNames = bioma.util.setDimNames(names,...
            obj.NSamples, 2, 'Sample', false);
    end
catch ME
    bioinfoprivate.bioclsrethrow('ExptData','ExptData', ME);
end
end % SampleName

function obj = setDMNames(obj, names)
try
    if ischar(names) && size(names, 1) == 1 && obj.NElements == 1
        obj.DMNames = {names};
    else
        obj.DMNames = bioma.util.setDimNames(names,...
            obj.NElements, 2, 'DMName', true);
    end
catch ME
    bioinfoprivate.bioclsrethrow('ExptData','ExptData', ME);
end
end % DMName

function dm = getElementDM(obj, idx, fcnName)
% Return a DataMatrix of a specified element.
elmtIdx = bioma.util.findLabelIndices(obj.ElementNames,idx);

try
    if isempty(elmtIdx)
        error(message('bioinfo:ExptData:ExptData:NotSuchElement')) 
    elseif numel(elmtIdx) > 1
        error(message('bioinfo:ExptData:ExptData:TooManyIndices'))
    end
catch ME
    bioinfoprivate.bioclsrethrow('ExptData',fcnName, ME);
end

dm = bioma.data.DataMatrix(obj.Matrices{elmtIdx},...
    obj.FeatureNames, obj.SampleNames, 'Name', obj.DMNames{elmtIdx});

end

function obj = setElementDM(obj, idx, b, fcnName)
% Assign a DataMatrix object b to the specified element in obj.

elmtIdx = bioma.util.findLabelIndices(obj.ElementNames,idx);

if numel(elmtIdx) > 1
    try
        error(message('bioinfo:ExptData:ExptData:TooManyIndices'))
    catch ME
        bioinfoprivate.bioclsrethrow('ExptData',fcnName, ME);
    end
elseif isempty(elmtIdx)
    if ischar(idx) % Assign a new element
        obj = validateSingleArg(b, obj);
        
        if isa(b, 'bioma.data.DataMatrix')
            obj.Matrices{end+1} = b.(':')(':');
            obj.DMNames{end+1} = b.Name;
        else
            obj.Matrices{end+1} = b;
        end
        
        obj.ElementNames{end+1} = idx;
    else
        try
           error(message('bioinfo:ExptData:ExptData:NotSuchElement')) 
        catch ME
           bioinfoprivate.bioclsrethrow('ExptData',fcnName, ME);
        end
    end
else
    if isnumeric(b)&&isempty(b)&&~isobject(b)&&~isstruct(b) % is Square Bracket ?
        obj.Matrices(elmtIdx) = [];
        obj.ElementNames(elmtIdx) = [];
        obj.DMNames(elmtIdx) = [];
    elseif isa(b, 'bioma.data.DataMatrix')
        obj = validateSingleArg(b, obj);
        obj.Matrices{elmtIdx} = b.(':')(':');
        obj.DMNames{elmtIdx} = b.Name;
    elseif isnumeric(b) || islogical(b)
        obj = validateSingleArg(b, obj);
        obj.Matrices{elmtIdx} = b;
        obj.DMNames{elmtIdx} = '';
    else
        try
            error(message('bioinfo:ExptData:ExptData:InvalidAssignment'));
        catch ME
            bioinfoprivate.bioclsrethrow('ExptData',fcnName, ME);
        end
    end
    
    if isempty(obj.FeatureNames) || isempty(obj.FeatureNames{1})
        obj = obj.featureNames(':', true);
    end
    
    if isempty(obj.SampleNames) || isempty(obj.SampleNames{1})
        obj = obj.sampleNames(':', true);
    end
end
end

function obj = updateNameProperties(obj, propLabel, idx, newNames)
%UPDATENAMEPROPERTIES Update the indexed names in a name property.
% 
%   NEWOBJ = UPDATENAMEPROPERTIES(OBJ, PROPLABEL, IDX, NEWNAMES) update the
%   names in a name property PROPLABEL indexed by IDX and return a new
%   object.  

idx = bioma.util.findLabelIndices(obj.(propLabel), idx);

if ischar(newNames)
    newNames = cellstr(newNames);
end

N = numel(idx);

if iscellstr(newNames)
    if isempty(newNames) || numel(newNames) ~= N || (numel(unique(newNames))~=numel(newNames))
        error(message('bioinfo:ExptData:updateNameProperties:InvalidInput', propLabel));
    end
else
    error(message('bioinfo:ExptData:updateNameProperties:InvalidInputType', propLabel));
end

tNames = obj.(propLabel);
tNames(idx) = newNames;

obj = obj.(propLabel)(':', tNames);

end
