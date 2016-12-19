classdef ExpressionSet < bioma.AExperiment
    % ExpressionSet Class to contain microarray gene expression experiment data.
    %
    %    The ExpressionSet class is a data container for data from a
    %    microarray gene expression experiment. It contains at least one
    %    DataMatrix object named Expressions for the gene expression
    %    values. This class is derived from the bioma.AExperiment abstract
    %    class.
    %
    %    ExpressionSet properties:
    %       NFeatures       - The number of features in the experiment.
    %       NSamples        - The number of samples in the experiment.
    %       NElements       - The number of elements (DataMatrix objects) in experiment data.
    %
    %   ExpressionSet methods:
    %       ExpressionSet       - Create an ExpressionSet object.
    %       expressions         - Access and set the Expressions DataMatrix object.
    %       exptData            - Access and set experiment data.
    %       elementData         - Access and set a DataMatrix element in experiment data.
    %       elementNames        - Access and set element names of experiment data.
    %       featureNames        - Access and set feature names.
    %       featureVarNames     - Access and set variable names describing the features.
    %       featureData         - Access and set the feature metadata.
    %       featureVarValues    - Access and set feature variable values.
    %       featureVarDesc      - Access and set the metadata describing variables in feature data.
    %       sampleNames         - Access and set sample names.
    %       sampleVarNames      - Access and set variable names describing the samples.
    %       sampleData          - Access and set sample metadata.
    %       sampleVarValues     - Access and set sample variable values.
    %       sampleVarDesc       - Access and set the metadata describing variables in sample data.   
    %       exptInfo            - Access and set the experiment information.
    %       abstract            - Access and set abstract in experiment information.
    %       pubMedID            - Access and set pubMed IDs in experiment information.
    %       exprWrite           - Write values from Expressions element to a text file.
    %       size                - Return the size of experiment data.
    %
    %   See also BIOMA.DATA.DATAMATRIX, DATASET, BIOMA.DATA.EXPTDATA,
    %   BIOMA.DATA.METADATA, BIOMA.DATA.MIAME. 
    
    %   Copyright 2009-2012 The MathWorks, Inc.

    
    
    methods
        function obj = ExpressionSet(data, varargin)
            %EXPRESSIONSET Create an object of ExpressionSet class. 
            %
            %   E = EXPRESSIONSET(DATA) creates an object E of
            %   ExpressionSet class from a DataMatrix or numeric array
            %   DATA. An ExptData object is implicitly created to contain
            %   DATA with an element name of 'Expressions'. DATA can also
            %   be an instance of the ExptData class. If there is no
            %   element named 'Expressions' in DATA, the first element in
            %   ExptData will be renamed to be 'Expressions'.
            %
            %   E = EXPRESSIONSET(DATA, {DM1, 'NAME'}, ...) creates an
            %   ExpressionSet object, E, with other named DataMatrix
            %   arguments, DMx, with the same dimensions as the
            %   'Expresssions' DataMatrix.
            %
            %   E = EXPRESSIONSET(..., 'SDATA', SD) creates an
            %   ExpressionSet object E, that contains SD, an object of the
            %   MetaData class, containing the metadata about the samples
            %   used in the experiment.
            %
            %   E = EXPRESSIONSET(..., 'FDATA', FD) creates an
            %   ExpressionSet object E, that contains FD, an object of the
            %   MetaData class, containing the metadata about the array
            %   features unique to the experiment.
            %
            %   E = EXPRESSIONSET(..., 'EINFO', INFO) creates an
            %   ExpressionSet object E, that contains details of experiment
            %   methods in, INFO, an object of the MIAME class.
            %
            %   See also EXPRESSIONSET.
            
            exprsName = 'Expressions';
            if nargin == 0
                superData = bioma.data.ExptData({bioma.data.DataMatrix,exprsName});
                pvPairs = {};
            else
                %== Handle input expression Data, must have a DataMatrix
                %   element named 'Expressions'
                if isa(data, 'bioma.data.DataMatrix')
                    superData = bioma.data.ExptData({data,exprsName});
                elseif isa(data, 'bioma.data.ExptData')
                    if isempty(find(ismember(elementNames(data), exprsName), 1))
                        superData = data.elementNames(1, exprsName);
                    else
                        superData = data;
                    end
                else
                    try
                        validateattributes(data, {'numeric', 'logical'}, {'real', '2d'});
                        superData = bioma.data.ExptData({data,exprsName});
                    catch ME 
                        error(message('bioinfo:ExpressionSet:ExpressionSet:InvalidExprsInput'))
                    end
                end
                
                %== Handle extra data input
                extraDataCount = 0;
                argCount = 1;
                nInput = nargin - 1;
                
                while argCount <= nInput
                    arg = varargin{argCount};
                    %== Guess if the input is param name/value pairs
                    if bioma.util.isString(arg)
                        break;
                    else
                        extraDataCount = extraDataCount +1;
                    end
                    argCount = argCount + 1;
                end
                
                if extraDataCount > 0
                    try
                        extraData = bioma.data.ExptData(varargin{1:extraDataCount});
                        superData = combine(superData, extraData);    
                    catch ME
                        bioinfoprivate.bioclsrethrow('ExpressionSet','ExpressionSet', ME)
                    end
                end
                pvPairs = varargin(extraDataCount+1 : end);
            end
            obj = obj@bioma.AExperiment(superData, pvPairs{:});
        end
        
    end
    
    methods
        function exprs = expressions(obj, dmatrix)
            %EXPRESSIONS Access and set the Expressions DataMatrix object.
            %
            %   EXPRS = EXPRESSIONS(E) returns the DataMatrix named
            %   Expressions of ExpressionSet object E.
            %
            %   B = EXPRESSIONS(E, EXPRS) returns an ExpressionSet object
            %   B, with the DataMatrix element named Expressions replaced
            %   by DataMatrix object EXPRS.
            %
            %   See also EXPRESSIONSET.
            try
                if nargin < 2
                    exprs = obj.elementData('Expressions');
                else
                    exprs = obj.elementData('Expressions', dmatrix);  
                end
            catch ME
                bioinfoprivate.bioclsrethrow('ExpressionSet', 'expressions', ME)
            end
        end
        
        function out = featureNames(obj, varargin)
            %FEATURENAMES Access and set the feature names of an ExpressionSet.
            %
            %   FNAMES = FEATURENAMES(E) returns a cell array of all
            %   feature names in ExpressionSet object E.
            %
            %   FNAMES = FEATURENAMES(E, I) returns the names of features
            %   specified by I of ExpressionSet object E. I can be a
            %   positive integer, a vector of positive integers, a string
            %   specifying a feature name, a cell array containing one or
            %   more feature names, or a logical vector.
            %
            %   B = FEATURENAMES(E, I, NAMES) returns an ExpressionSet
            %   object B with the names of specified features set to NAMES.
            %   The number of names in NAMES must equal the number of
            %   features specified by I. NAMES can be a cell array of
            %   strings, a character array, or a numeric vector. NAMES can
            %   also be a single string as a prefix for feature names;
            %   numbers will be appended to the prefix. NAMES can also be a
            %   logical true or false (default); if true, default unique
            %   names will be assigned to the features.
            %
            %   See also EXPRESSIONSET, FEATUREVARNAMES, SAMPLENAMES,
            %   SAMPLEVARNAMES. 
            
            try
                rel = featureNames(obj.Data, varargin{:});
                if isa(rel, 'bioma.data.ExptData')
                    eFeatureData = obj.FeatureData.sampleNames(':', rel.featureNames);
                    out = bioma.ExpressionSet(rel,...
                        'SData', obj.SampleData,...
                        'FData', eFeatureData,...
                        'Einfo', obj.ExptInfo);
                else
                    out = rel;
                end
            catch ME
                bioinfoprivate.bioclsrethrow('ExpressionSet', 'featureNames', ME)
            end
        end
        
        function out = sampleNames(obj, varargin)
            %SAMPLENAMES Access and set the sample names of an ExpressionSet.
            %
            %   SNAMES = SAMPLENAMES(E) returns a cell array of all sample
            %   names in ExpressionSet object E.
            %
            %   SNAMES = SAMPLENAMES(E, I) returns the names of samples
            %   specified by I of ExpressionSet object E. I can be a
            %   positive integer, a vector of positive integers, a string
            %   specifying a sample name, a cell array containing one or
            %   more sample names, or a logical vector.
            %
            %   B = SAMPLENAMES(E, I, NAMES) returns an ExpressionSet
            %   object B with the names of specified samples set to NAMES.
            %   The number of names in NAMES must equal the number of
            %   samples specified by I. NAMES can be a cell array of
            %   strings, a character array, or a numeric vector. NAMES can
            %   also be a single string as a prefix for sample names;
            %   numbers will be appended to the prefix. NAMES can also be a
            %   logical true or false (default); if true, default unique
            %   names will be assigned to the samples.
            %
            %   See also EXPRESSIONSET, FEATURENAMES, FEATUREVARNAMES,
            %   SAMPLEVARNAMES. 
            try
                rel = sampleNames(obj.Data, varargin{:});
                if isa(rel, 'bioma.data.ExptData')
                    eSampleData = obj.SampleData.sampleNames(':', rel.sampleNames);
                    out = bioma.ExpressionSet(rel,...
                        'SData', eSampleData,...
                        'FData', obj.FeatureData,...
                        'Einfo', obj.ExptInfo);
                else
                    out = rel;
                end
            catch ME
                bioinfoprivate.bioclsrethrow('ExpressionSet', 'sampleNames', ME)
            end
        end
        
        function out = elementNames(obj,varargin)
            %ELEMENTNAMES Access and set element names of experiment data.
            %
            %   ENAMES = ELEMENTNAMES(E) returns a cell array of all
            %   element names in the experiment data in ExpressionSet
            %   object E.  
            %
            %   ENAMES = ELEMENTNAMES(E, I) returns the names of elements
            %   specified by I experiment data in ExpressionSet object E. I
            %   can be a positive integer, a vector of positive integers, a
            %   logical vector, a string, or a cell array of strings.
            %
            %   B = ELEMENTNAMES(E, I, NAMES) returns an ExpressionSet
            %   object B with the names of specified elements in the
            %   experiment data set to NAMES. The number of names in NAMES
            %   must equal the number of elements specified by I. NAMES can
            %   be a cell array of strings, a character array, or a numeric
            %   vector. NAMES can also be a single string as a prefix for
            %   element names; numbers will be appended to the prefix.
            %   NAMES can also be a logical true or false (default); if
            %   true, default unique names will be assigned to the
            %   elements.
            %
            %   See also EXPRESSIONSET, EXPTDATA.
            try
                rel = elementNames(obj.Data, varargin{:});
              
                if isa(rel, 'bioma.data.ExptData')                    
                    out = obj;
                    out.Data = elementNames(rel, 1, 'Expressions');
                else
                    out = rel;
                end
            catch ME
                bioinfoprivate.bioclsrethrow('ExpressionSet', 'elementNames', ME)
            end
        end
        
        function out = featureVarNames(obj, varargin)
            %FEATUREVARNAMES Access and set the variable names of feature data.
            %
            %   VNAMES = FEATUREVARNAMES(E) returns a cell array of all
            %   variable names of feature data in ExpressionSet object E.
            %
            %   VNAMES = FEATUREVARNAMES(E, I) returns the names of
            %   variables specified by I in ExpressionSet object E. I can
            %   be a positive integer, a vector of positive integers, a
            %   string specifying a variable name, a cell array containing
            %   one or more variable names, or a logical vector.
            %
            %   B = FEATUREVARNAMES(E, I, NAMES) returns an ExpressionSet
            %   object B with the names of specified feature variables set
            %   to NAMES. The number of names in NAMES must equal the
            %   number of variables specified by I. NAMES can be a cell
            %   array of strings, a character array, or a numeric vector.
            %   NAMES can also be a single string as a prefix for variable
            %   names; variable numbers will be appended to the prefix.
            %   NAMES can also be a logical true or false (default); if
            %   true, default unique names will be assigned to the
            %   variables. The variable names (as rows) in the variable
            %   meta information dataset are updated automatically.
            %
            %   See also EXPRESSIONSET, FEATURENAMES, SAMPLEVARNAMES.
            try
                rel = variableNames(obj.FeatureData, varargin{:});
                if isa(rel, 'bioma.data.MetaData')
                    out = obj;
                    out.FeatureData = rel;
                else
                    out = rel;
                end
            catch ME
                bioinfoprivate.bioclsrethrow('ExpressionSet', 'featureVarNames', ME)
            end
        end
        
        function out = sampleVarNames(obj, varargin)
            %SAMPLEVARNAMES Access and set the variable names of sample data.
            %
            %   VNAMES = SAMPLEVARNAMES(E) returns a cell array ofall  variable
            %   names of sample data in an ExpressionSet object E.
            %
            %   VNAMES = SAMPLEVARNAMES(E, I) returns the names of sample
            %   variables specified by I in ExpressionSet E. I can be a
            %   positive integer, a vector of positive integers, a string
            %   specifying a variable name, a cell array containing one or
            %   more variable names, or a logical vector.
            %
            %   B = SAMPLEVARNAMES(E, I, NAMES) returns an ExpressionSet
            %   object B with the names of specified sample variables set
            %   to NAMES. The number of names in NAMES must equal the
            %   number of variables specified by I. NAMES can be a cell
            %   array of strings, a character array, or a numeric vector.
            %   NAMES can also be a single string as a prefix for variable
            %   names; variable numbers will be appended to the prefix.
            %   NAMES can also be a logical true or false (default); if
            %   true, default unique names will be assigned to the
            %   variables. The variable names (as rows) in the variable
            %   meta information dataset are updated automatically.
            %
            %   See also FEATURENAMES, FEATUREVARNAMES, SAMPLENAMES.
            try
                rel = variableNames(obj.SampleData, varargin{:});
                if isa(rel, 'bioma.data.MetaData')
                    out = obj;
                    out.SampleData = rel;
                else
                    out = rel;
                end
            catch ME
                bioinfoprivate.bioclsrethrow('ExpressionSet', 'sampleVarNames', ME)
            end
        end
        
        function out = exptData(obj,varargin)
            %EXPTDATA Access and set the experiment data in an ExpressionSet object.
            %
            %   D = EXPTDATA(E) returns an ExptData object D stored in
            %   ExpressionSet object E.
            %
            %   B = EXPTDATA(E, D) returns an ExpressionSet B with its
            %   experiment data replaced by ExptData object D. The sample
            %   names and feature names in D must match the existing sample
            %   names and feature names in ExpressionSet object E.
            %
            %   See also EXPRESSIONSET, FEATUREDATA, SAMPLEDATA.
            try
                if nargin < 2
                    out = obj.Data;
                else
                    out = obj;
                    out.Data = varargin{1};
                    out = elementNames(out, 1, 'Expressions');
                end
            catch ME
                bioinfoprivate.bioclsrethrow('ExpressionSet', 'exptData', ME)
            end
        end
        
        function varargout = size(a,varargin)
            %SIZE Size of an ExpressionSet object.
            %
            %   D = SIZE(E) returns a two-elment row vector D containing
            %   the number of features and samples in ExpressionSet object
            %   E.
            %
            %   [M,N] = SIZE(E) returns the number of features and columns
            %   in ExpressionSet object E as separate output variables.
            %
            %   M = SIZE(E,DIM) returns the length of the dimension
            %   specified by the scalar DIM. For example, SIZE(E,1) returns
            %   the number of features. If DIM > NDIMS(E), M will be 1.
            %
            %   See also EXPRESSIONSET.
            
            try
                varargout{1:nargout} = a.Data.size(varargin{:});
            catch ME
                bioinfoprivate.bioclsrethrow('ExpressionSet', 'size', ME)
            end
        end % size method
        
        function exprWrite(obj, filename, varargin)
            %EXPRWRITE Write values from Expressions element to a text file.
            %
            %   EXPRWRITE(E, FILENAME) writes expression values in the
            %   Expressions element (a DataMatrix object) of ExpressionSet
            %   object E to a text file, FILENAME, using default delimiter
            %   (\t) to separate columns in the DataMatrix object. The data
            %   is written starting at the first column of the first row in
            %   the destination file. FILENAME must be a string.
            %
            %   EXPRWRITE(..., 'DELIMITER', DELIMITER) specifies delimiter
            %   DELIMITER to separate matrix columns. '\t' is the default
            %   delimiter.
            %
            %   EXPRWRITE(..., 'PRECISION',PS) specifies the numeric
            %   precision PS used in writing data to the file. PS can be
            %   the number of significant digits or a C-style format string
            %   starting in %, such as '%6.5f'. Default is 5.
            %
            %   EXPRWRITE(..., 'HEADER', HEADERTEXT) specifies the first
            %   line of the file. The default is the name in the Name
            %   property of the Expressions element.
            %
            %   EXPRWRITE(..., 'ANNOTATED',FALSE) does not write row and
            %   column names to the file. Default is TRUE.
            %
            %   EXPRWRITE(..., 'APPEND',TRUE) appends the Expressions
            %   DataMatrix object to the end of the existing file. Default
            %   is FALSE.
            %
            %   See also EXPRESSIONSET, BIOMA.DATA.DATAMATRIX.DMWRITE.
            
            exprs = expressions(obj);
            try 
                dmwrite(exprs, filename, varargin{:});
            catch ME
                 bioinfoprivate.bioclsrethrow('ExpressionSet', 'exprWrite', ME);
            end
        end 
    end %end of method

end
