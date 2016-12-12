classdef AExperiment
    % AEXPERIMENT Abstract class to contain for microarray experiment data.
    %
    %    AExperiment is an abstract class for coordinate and contain
    %    microarray experiment data. It contains data measured from a
    %    microarray experiment as an ExptData object, metadata about the
    %    samples as an MetaData object, data about the features the arrays
    %    as a MetaData object, and experiment description as a MIAME
    %    object.
    %
    %    Classes derived from AExperiment contain one or more
    %    identical-size DataMatrix objects as ExptData elements with column
    %    number equals to number of rows in sample data, and row number
    %    equals to the number of rows in feature data. 
    %
    %    AExperiment properties:
    %       NFeatures       - The number of features in the experiment.
    %       NSamples        - The number of samples in the experiment.
    %       NElements       - The number of DataMatrix objects in experiment data.
    %
    %   AExperiment methods:
    %       exptData            - Abstract method to access and set experiment data.
    %       elementNames        - Abstract method to access and set element names of experiment data.
    %       featureNames        - Abstract method to access and set feature names.
    %       featureVarNames     - Abstract method to access and set variable names of feature data.
    %       sampleNames         - Abstract method to access and set sample names.
    %       sampleVarNames      - Abstract method to access and set variable names of sample data.
    %       size                - Abstract method to return the size of experiment data.  
    %       elementData         - Access and set a DataMatrix element in experiment data.
    %       featureData         - Access and set the feature metadata.
    %       featureVarValues    - Access and set feature variable values.
    %       featureVarDesc      - Access and set meta data about variables in feature data.
    %       sampleData          - Access and set sample metadata.
    %       sampleVarValues     - Access and set sample variable values.
    %       sampleVarDesc       - Access and set meta data about variables in sample data.
    %       exptInfo            - Access and set the experiment information.
    %       abstract            - Access and set abstract in experiment information.
    %       pubMedID            - Access and set pubMed IDs in experiment information.
    %
    %   See also BIOMA.DATA.DATAMATRIX, DATASET, EXPRESSIONSET,
    %   BIOMA.DATA.EXPTDATA, BIOMA.DATA.METADATA, BIOMA.DATA.MIAME.
    
    %   Copyright 2009-2012 The MathWorks, Inc.

    
    properties (SetAccess = 'private', GetAccess = 'public')
        %NFEATURES The number of features in the experiment.
        %   The NFeatures property is a scalar equal to the number of
        %   features.
        %
        %   See also BIOMA.AEXPERIMENT.
        NFeatures = 0;
        
        %NSAMPLES The number of samples in the experiment.
        %   The NSamples property is a scalar equal to the number of
        %   samples.
        %
        %   See also BIOMA.AEXPERIMENT.
        NSamples = 0;
        
        %NELEMENTS The number of elements (DataMatrix objects) in experiment data.
        %   The NElements property is a scalar equal to the number of
        %   DataMatrix objects or elements in the experiment data.
        %
        %   See also BIOMA.AEXPERIMENT.
        NElements = 0;
    end
    
    properties (SetAccess = 'protected', GetAccess = 'protected', Hidden=true)
        %DATA Contains DataMatrix objects with equal dimensions.
        %    The Data property is an ExptData containing DataMatrix objects
        %    with equal dimension properties. One of the element name in
        %    Data is 'Expressions'.
        %
        %    See also BIOMA.AEXPERIMENT, BIOMA.DATA.EXPTDATA.
        Data = bioma.data.ExptData;
        
        %SAMPLEDATA Contains values of variables describing samples.
        %    The SampleData property is a MetaData containing values of
        %    variables  measured on the samples. 
        %
        %    See also BIOMA.AEXPERIMENT, BIOMA.DATA.METADATA.
        SampleData = bioma.data.MetaData;
        
        %FEATUREDATA Contains values of variables describing features.
        %    The FeatureData property is a MetaData containing values of
        %    variables measured for the features. 
        %
        %    See also BIOMA.AEXPERIMENT, BIOMA.DATA.METADATA.
        FeatureData = bioma.data.MetaData;
        
        %EXPTINFO Contains description about the experiment.
        %    The ExptInfo property is a MIAME containing details of
        %    experiment methods and conditions.
        %
        %    See also BIOMA.AEXPERIMENT, BIOMA.DATA.MIAME.
        ExptInfo = bioma.data.MIAME;
    end
    
    methods(Abstract)
        %EXPTDATA Abstract method to access and set experiment data.
        %
        %   See also BIOMA.AEXPERIMENT.
        out = exptData(obj,varargin);
        
        %ELEMENTNAMES Abstract method to access and set element names of experiment data.
        %
        %   See also BIOMA.AEXPERIMENT.
        out = elementNames(obj,varargin)
            
        
        %FEATURENAMES Abstract method to access and set feature names.
        %
        %   See also BIOMA.AEXPERIMENT.
        out = featureNames(obj, varargin);

        %SAMPLENAMES Abstract method to access and set sample names.
        %
        %   See also BIOMA.AEXPERIMENT.
        out = sampleNames(obj, varargin);
        
        %FEATUREVARNAMES Abstract method to access and set variable names of feature data.
        %
        %   See also BIOMA.AEXPERIMENT.
        out = featureVarNames(obj, varargin);
        
        %SAMPLEVARNAMES Abstract method to access and set variable names of sample data.
        %
        %   See also BIOMA.AEXPERIMENT.
        out = sampleVarNames(obj, varargin);

        %SIZE Abstract method to return the size of experiment data.
        %
        %   See also BIOMA.AEXPERIMENT.    
        varargout = size(a,varargin);
    end
    
    methods
        function obj = AExperiment(data, varargin)
            %AEXPERIMENT Create an object for containing experiment values and meta data.
            %
            %   A = AEXPERIMENT(DATA) creates an AExperiment A from an
            %   ExptData object DATA.
            %
            %   A = AEXPERIMENT(..., 'SDATA', SD) creates an AExperiment A
            %   that contains variables describing samples, SD of MetaData
            %   class.
            %
            %   A = AEXPERIMENT(..., 'FDATA', FD) creates an AExperiment A
            %   that contains variables describing features unique to the
            %   experiment, FD of MetaData class.
            %
            %   A = AEXPERIMENT(..., 'EInfo', ED) creates an AExperiment A
            %   that contains details of experiment methods and conditions,
            %   ED of MIAME class.
            %
            %   See also BIOMA.AEXPERIMENT.
            
            if nargin == 0
                return;
            end
            
            if isa(data, 'bioma.data.ExptData')
                obj.Data = data;
            else
                error(message('bioinfo:AExperiment:AExperiment:InvalidExptDataInput'))
            end
            
            obj = parsePVInputs(obj, varargin{:});
        end
        
        function out = elementData(obj,varargin)
            %ELEMENTDATA Access and set a DataMatrix element of the experiment data.
            %
            %   D = ELEMENTDATA(A, ELMT) returns a DataMatrix object D that
            %   is the specified element ELMT from experiment data in
            %   object A of AExperiment class.
            %
            %   B = ELEMENTDATA(A, ELMT, DM) returns an object B of
            %   AExperiment class with specified element of experiment data
            %   replaced by a DataMatrix object DM. The sample names and
            %   feature names of DM must match the existing sample names
            %   and feature names in A.
            %
            %   See also BIOMA.AEXPERIMENT, BIOMA.DATA.EXPTDATA.
            bioinfochecknargin(nargin, 2, ['AExperiment', ':elementData'])
            
            try
                if nargin < 3
                    out = obj.Data(varargin{1});
                else
                    out = obj;
                    out.Data(varargin{1}) = varargin{2};
                end
            catch ME
                bioinfoprivate.bioclsrethrow('AExperiment','elementData',ME)
            end
        end
        
        function out = sampleData(obj,varargin)
            %SAMPLEDATA Access and set sample metadata.
            %
            %   D = SAMPLEDATA(A) returns a MetaData object D containing
            %   sample information about object A of AExperiment class.
            %
            %   B = SAMPLEDATA(A, D) returns an object B of AExperiment
            %   class with its sample data replaced by a MetaData object or
            %   a dataset array D. The sample or row names in D must
            %   match the sample names in A.
            %
            %   See also BIOMA.AEXPERIMENT, FEATUREDATA, SAMPLENAMES.
            try
                if nargin < 2
                    out = obj.SampleData;
                else
                    out = obj;
                     if isa(varargin{1}, 'bioma.data.MetaData')
                        out.SampleData = varargin{1};
                     else
                         out.SampleData = variableValues(obj.SampleData, varargin{1});
                     end
                end
            catch ME
                bioinfoprivate.bioclsrethrow('AExperiment', 'sampleData', ME)
            end
        end
        
        function out = featureData(obj,varargin)
            %FEATUREDATA Access and set feature metadata.
            %
            %   D = FEATUREDATA(A) returns a MetaData object D containing
            %   feature information about object A of AExperiment class.
            %
            %   B = FEATUREDATA(A, D) returns an object B of the
            %   AExperiment class with its feature data replaced by
            %   MetaData object or a dataset array D. The feature or row
            %   names in D must match the feature names in A.
            %
            %   See also BIOMA.AEXPERIMENT, FEATURENAMES, SAMPLEDATA.
            try
                if nargin < 2
                    out = obj.FeatureData;
                else
                    out = obj;
                    if isa(varargin{1}, 'bioma.data.MetaData')
                        out.FeatureData = varargin{1};
                    else
                        out.FeatureData = variableValues(obj.FeatureData, varargin{1});
                    end
                end
            catch ME
                bioinfoprivate.bioclsrethrow('AExperiment', 'featureData', ME)
            end
        end
        
        function out = sampleVarValues(obj,varargin)
            %SAMPLEVARVALUES Access and set sample variable values.
            %
            %   DS = SAMPLEVARVALUES(A) returns a dataset array, DS,
            %   containing the measured values of variable data about
            %   samples in object A of AExperiment class.
            %
            %   B = SAMPLEVARVALUES(A, D) returns an object B of
            %   AExperiment class with its sample variable data replaced
            %   with dataset array D. The sample or row names in D must
            %   match the sample names in A.
            %
            %   See also BIOMA.AEXPERIMENT,
            %   BIOMA.DATA.METADATA.VARIABLEVALUES. 
            
            try
                rel = variableValues(obj.SampleData, varargin{:});
                if isa(rel, 'bioma.data.MetaData')
                    out = obj;
                    out.SampleData = rel;
                else
                    out = rel;
                end
            catch ME
                bioinfoprivate.bioclsrethrow('AExperiment', 'sampleVarValues', ME)
            end
        end
        
        function out = featureVarValues(obj,varargin)
            %FEATUREVARVALUES Access and set feature variable values.
            %
            %   DS = FEATUREVARVALUES(A) returns a dataset array, DS,
            %   containing the measured values of variable data about
            %   features in object A of AExperiment class.
            %
            %   B = FEATUREVARVALUES(A, D) returns an object B of
            %   AExperiment class with its feature variable data replaced
            %   with dataset array D. The feature or row names in D must
            %   match the feature names in A.
            %
            %   See also BIOMA.AEXPERIMENT,
            %   BIOMA.DATA.METADATA.VARIABLEVALUES.
            
            try
                rel = variableValues(obj.FeatureData, varargin{:});
                if isa(rel, 'bioma.data.MetaData')
                    out = obj;
                    out.FeatureData = rel;
                else
                    out = rel;
                end
            catch ME
                bioinfoprivate.bioclsrethrow('AExperiment', 'featureVarValues', ME)
            end
        end
 
        function out = sampleVarDesc(obj,varargin)
            %SAMPLEVARDESC Access and set meta information about sample variables.
            %
            %   DS = SAMPLEVARDESC(A) returns a dataset array, DS, with
            %   variable names as rows, and description labels as columns
            %   for sample data in object A of AExperiment class.
            %
            %   B = SAMPLEVARDESC(A, DS) returns an object B of AExperiment
            %   class with its sample variable metadata set to dataset DS.
            %   DS can also be a cell array of strings containing
            %   descriptions of the variables.
            %
            %   See also BIOMA.AEXPERIMENT,
            %   BIOMA.DATA.METADATA.VARIABLEDESC. 
            
            try
                rel = variableDesc(obj.SampleData, varargin{:});
                if isa(rel, 'bioma.data.MetaData')
                    out = obj;
                    out.SampleData = rel;
                else
                    out = rel;
                end
            catch ME
                bioinfoprivate.bioclsrethrow('AExperiment', 'sampleVarDesc', ME)
            end
        end
        
        function out = featureVarDesc(obj,varargin)
            %FEATUREVARDESC Access and set meta information about feature variables.
            %
            %   DS = FEATUREVARDESC(A) returns a dataset array, DS, with
            %   variable names as rows, description labels as columns of
            %   feature data in object A of AExperiment class.
            %
            %   B = FEATUREVARDESC(A, DS) returns an object B of
            %   AExperiment class with its feature variable meta
            %   information data set to a dataset DS. DS can also be a cell
            %   array of strings containing descriptions of the variables.
            %
            %   See also BIOMA.AEXPERIMENT,
            %   BIOMA.DATA.METADATA.VARIABLEDESC.
            
            try
                rel = variableDesc(obj.FeatureData, varargin{:});
                if isa(rel, 'bioma.data.MetaData')
                    out = obj;
                    out.FeatureData = rel;
                else
                    out = rel;
                end
            catch ME
                bioinfoprivate.bioclsrethrow('AExperiment', 'featureVarDesc', ME)
            end
        end
        
        function out = exptInfo(obj, mv)
            %EXPTINFO Access and set description about the experiment.
            %
            %   D = EXPTINFO(A) returns a MIAME object D from object A of
            %   AExperiment class. 
            %
            %   B = EXPTINFO(A, M) returns an object B of AExperiment class
            %   with its experiment information updated with information
            %   from MIAME object M.
            %
            %   See also BIOMA.AEXPERIMENT, BIOMA.DATA.MIAME.
            try
                if nargin < 2
                    out = obj.ExptInfo;
                else
                    out = obj;
                    out.ExptInfo = mv;
                end
            catch ME
                bioinfoprivate.bioclsrethrow('AExperiment', 'exptInfo', ME)
            end
        end
        
        function out = pubMedID(obj, ids)
            %PUBMEDID Access and set PubMed IDs in experiment information.
            %
            %   D = PUBMEDID(A) returns PubMed IDs in the experiment
            %   information from object A of AExperiment class.
            %
            %   B = PUBMEDID(A, IDs) returns an object B of AExperiment
            %   class with its PubMed IDs in the experiment information
            %   updated to IDS.
            %
            %   See also BIOMA.AEXPERIMENT, BIOMA.DATA.MIAME.
            
            try
                if nargin < 2
                    out = obj.ExptInfo.PubMedID;
                else
                    out = obj;
                    out.ExptInfo.PubMedID = ids;
                end
            catch ME
                bioinfoprivate.bioclsrethrow('AExperiment', 'pubMedID', ME)
            end
        end
        
        function out = abstract(obj, abs)
            %ABSTRACT Access and set the abstract in experiment information.
            %
            %   D = ABSTRACT(A) returns the abstract in the experiment
            %   information from object A of AExperiment class.
            %
            %   B = ABSTRACT(A, ABS) returns an object B of AExperiment
            %   class with its abstract in the experiment information
            %   updated to string ABS.
            %
            %   See also BIOMA.AEXPERIMENT, BIOMA.DATA.MIAME.
            
            try
                if nargin < 2
                    out = obj.ExptInfo.Abstract;
                else
                    out = obj;
                    out.ExptInfo.Abstract = abs;
                end
            catch ME
                bioinfoprivate.bioclsrethrow('AExperiment', 'abstract', ME)
            end
        end
    end
    
    methods %set methods
        function n = get.NSamples(obj)
            n = obj.Data.NSamples;
        end
        
        function n = get.NFeatures(obj)
            n = obj.Data.NFeatures;
        end
        
        function n = get.NElements(obj)
            n = obj.Data.NElements;
        end
        
        function obj = set.Data(obj, data)
            if isa(data, 'bioma.data.ExptData')
                if ~isempty(obj.Data)
                    try
                        if ~isempty(obj.SampleData) %#ok
                            validateNames(sampleNames(obj.Data), data.sampleNames,...
                                'Input variable', 'existing samples');
                        end
                        if ~isempty(obj.FeatureData) %#ok
                            validateNames(featureNames(obj.Data), data.featureNames,...
                                'Input variable', 'existing features');
                        end
                    catch ME
                        bioinfoprivate.bioclsrethrow('AExperiment','AExperiment', ME);
                    end
                end
                obj.Data = data;
            else
                error(message('bioinfo:AExperiment:AExperiment:InvalidExprDataInput'))
            end
        end
        
        function obj = set.SampleData(obj, sMData)
            if isempty(sMData)
                obj.SampleData = bioma.data.MetaData;
            elseif isa(sMData, 'bioma.data.MetaData')
                try
                    validateNames(sampleNames(obj.Data), sMData.sampleNames,...
                        'Input variable', 'existing samples'); %#ok
                    obj.SampleData = sMData;
                catch ME
                    bioinfoprivate.bioclsrethrow('AExperiment','AExperiment', ME);
                end
            else
                error(message('bioinfo:AExperiment:AExperiment:InvalidSampleDataInput'))
            end
        end
        
        function obj = set.FeatureData(obj, fMData)
            if isempty(fMData)
                obj.FeatureData = bioma.data.MetaData;
            elseif isa(fMData, 'bioma.data.MetaData')
                try
                    validateNames(featureNames(obj.Data), fMData.sampleNames,...
                        'Input variable', 'existing features'); %#ok
                    obj.FeatureData = fMData;
                catch ME
                    bioinfoprivate.bioclsrethrow('AExperiment','AExperiment', ME);
                end
            else
                error(message('bioinfo:AExperiment:AExperiment:InvalidFeatureDataInput'))
            end
        end
        
        function obj = set.ExptInfo(obj, miame)
            if isa(miame, 'bioma.data.MIAME')
                obj.ExptInfo = miame;
            else
                error(message('bioinfo:AExperiment:AExperiment:InvalidExptInfoInput'))
            end
        end
    end
    
    methods(Hidden=true)
        function obj = updateExpt(obj, data, type)
            %Update all the properties, if there is new data changes
            type = upper(type);
            switch(type)
                case 'SAMPLENAMES' %SampleNames
                    obj = sampleNames(obj.Data, data);
                    obj = obj.SampleData.rowNames(':', data);
                case 'FEATURENAMES' %FeatureNames
                    obj = featureNames(obj.Data, data);
                    obj = obj.FeatureData.rowNames(':', data);
            end
        end %updateExpt;
        
        function disp(obj)
            clsName = class(obj);
            clsIdx = strfind(clsName, '.');
            clsName = clsName(clsIdx(end)+1 : end);
            colPad = repmat(' ', 1, 2);
            
            fprintf('%s\n', clsName)
            
            fprintf('Experiment Data: %d features, %d samples\n',...
                obj.Data.NFeatures, obj.Data.NSamples);
            if obj.NElements == 0 || (obj.NElements == 1 && isempty(obj.Data(1)))
                fprintf('%s Element names: %s\n',...
                    colPad, 'none');
            else
                fprintf('%sElement names: %s\n',...
                    colPad, bioma.util.catCellStrToStr(obj.Data.elementNames));
                if isempty(obj.SampleData)
                    fprintf('Sample Data: none\n')
                else
                    fprintf('Sample Data:\n')
                    dispMetadata(obj.SampleData, 'Sample', 'MetaData');
                end
                
                if isempty(obj.FeatureData)
                    fprintf('Feature Data: none\n')
                else
                    fprintf('Feature Data:\n')
                    dispMetadata(obj.FeatureData, 'Feature', 'MetaData')
                end
            end
            
            if isempty(obj.ExptInfo)
                fprintf('Experiment Information: none\n');
            else
                fprintf('Experiment Information: use ''exptInfo(obj)''\n');
            end
        end % disp
        
        function [varargout] = subsref(a, s)
            %SUBSREF Subscripted reference for an AExperiment object.
            %
            %   B = A(I,J) returns an AExperiment B that contains a subset
            %   of the features and samples in the AExperiment object A. I
            %   and J can be a positive integer, a vector of positive
            %   integers, a string of names, a cell array containing one or
            %   more unique names, or a logical vector. B contains the same
            %   property values as A, subsetted for rows or columns where
            %   appropriate. () subscribe type does not allow linear
            %   indexing.
            %
            %   Note: Cell indexing A{I,J} is not supported.
            %
            %   P = A.PROPERTYNAME returns an ExptData property.
            switch s(1).type
                case '()'
                    %==Only one-dimension indexing
                    if nargout > 1
                        error(message('bioinfo:AExperiment:subsref:TooManyOutputs'));
                    elseif ~isscalar(s)
                        error(message('bioinfo:AExperiment:subsref:InvalidSubsExpr'));
                    elseif numel(s(1).subs) ~= 2
                        error(message('bioinfo:AExperiment:subsref:NDSubscript'));
                    end
                    
                    try
                        %Subscript for sample data
                        s_spl = s;
                        s_spl.subs = {s.subs{2},':'};
                        %subscript for feature data
                        s_ftr = s;
                        s_ftr.subs = {s.subs{1},':'};
                        b = a;
                        b.SampleData = [];
                        b.FeatureData = [];
                        
                        b.Data = subsref(a.Data,s);
                        if ~isempty(a.SampleData)
                            b.SampleData = subsref(a.SampleData, s_spl);
                        end
                        
                        if ~isempty(a.FeatureData)
                            b.FeatureData = subsref(a.FeatureData, s_ftr);
                        end
                        [varargout{1:nargout}] = b;
                    catch ME
                        if strfind(ME.identifier, 'IndexOutOfBounds')
                            error(message('bioinfo:AExperiment:subsref:IndexOutOfBounds', a.NElements));
                        else
                            bioinfoprivate.bioclsrethrow('AExperiment','subsref', ME)
                        end
                    end
                case '{}'
                    error(message('bioinfo:AExperiment:subsref:CellSubscript'));
                case '.'
                    varName = s(1).subs;
                    if ischar(varName) && size(varName, 1) == 1
                        try
                            if isempty(find(strcmp(varName, a.sampleVarNames), 1))
                                bioma.util.validateMethodProps(a, varName);
                                [varargout{1:nargout}] = builtin('subsref',a,s);
                            else
                                [varargout{1:nargout}] = subsref(a.SampleData, s);
                            end
                        catch ME
                            bioinfoprivate.bioclsrethrow('AExperiment','subsref', ME)
                        end
                    else
                        error(message('bioinfo:AExperiment:subsref:IllegalVarSubscript'));
                    end
            end
        end %subsref
    end
end %AExperiment

%------Helper functions-----------
function validateNames(validNames, namesToCheck, vldName, chkName)
if ischar(validNames)
    validNames = cellstr(validNames);
end
nameIdx = bioma.util.findLabelIndices(validNames, namesToCheck);
if numel(nameIdx) ~= numel(validNames)
    error(message('bioinfo:AExperiment:AExperiment:NumberOfNamesNotMatch', vldName, chkName));
end
end

function obj = parsePVInputs(obj, varargin)
% Parse input PV pairs.

% Check for the right number of inputs
if rem(nargin-1, 2)== 1
    error(message('bioinfo:AExperiment:AExperiment:IncorrectNumberOfArguments'))
end
validNames = {'sdata', 'fdata', 'einfo'};
for j=1:2:nargin-1
    k = bioinfoprivate.pvpair(varargin{j},[],validNames,'AExperiment:AExperiment');
    switch(k)
        case 1 % sample data
            obj.SampleData = varargin{j+1};
        case 2 % feature data
            obj.FeatureData = varargin{j+1};
        case 3 % experiment information
            obj.ExptInfo = varargin{j+1};
    end
end
end %parse_input

function dispMetadata(a, propName, clsName)
%DISP Display AExperiment object.
%
%   DISP(A) displays the AExperiment A, including the first two and last row
%   names, a table formatted variable labels and the variable metadata,
%   without printing the MetaData object name.

colPad = repmat(' ', 1, 4);
rowNames = a.sampleNames;

amChars = bioma.util.printOneLineMultiNames(rowNames, colPad,...
    [a.NSamples, a.NVariables], clsName);

%== Display row names
if isempty(rowNames)
    amChars = '';
end

fprintf('%s%s names: %s\n', colPad, propName, amChars)

varMetaInfo = a.variableDesc;
if isempty(varMetaInfo)
    fprintf('%s%s variable names and meta information: %s\n',...
            colPad,propName, 'none')
else
    fprintf('%s%s variable names and meta information: \n',...
            colPad, propName)
    varNames = get(a.variableValues, 'VarNames');
    
    for i = 1:a.NVariables
        if isempty(varMetaInfo) || isempty(varMetaInfo(i, 'VariableDescription'))
            varDesc = 'NA';
        else
            varDesc = varMetaInfo{i, 'VariableDescription'};
        end
        fprintf('%s%s%s: %s\n',colPad, colPad,varNames{i}, varDesc)
    end
end
end % disp method
