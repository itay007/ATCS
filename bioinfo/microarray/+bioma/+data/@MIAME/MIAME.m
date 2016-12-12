classdef MIAME
    %MIAME Class for storing information about a microarray experiment.
    %
    %   The MIAME class is a collection of microarray experiment
    %   information loosely following the Minimum Information About a
    %   Microarray Experiment (MIAME) specification. The minimum experiment
    %   information including:
    %       1. Experiment design
    %       2. Samples used, extract preparation and labeling
    %       3. Hybridization procedures and parameters
    %       4. Measurement data and specifications of data processing
    %   The MIAME class stores information about experiment design,
    %   samples, hybridyzations, normalization controls, and pre-processing
    %   information.
    % 
    %   MIAME properties:
    %       Investigator    - Name of the investigator of the experiment.        
    %       Laboratory      - The laboratory where the experiment was conducted.
    %       Contact         - Contact information for laboratory or investigator.
    %       URL             - A URL for the experiment.
    %       Title           - An experiment title.
    %       Abstract        - An abstract describing the experiment.
    %       ExptDesign      - A summary about the experiment design.
    %       Arrays          - Information about the arrays.
    %       Samples         - Information about the samples.
    %       Hybridization   - Information about the hybridization.
    %       QualityControl  - Information about quality controls.
    %       PreProcessing   - Information about pre-processing steps. 
    %       PubMedID        - PubMed identifiers of relevant publications.
    %       Other           - Other information about the experiment.
    %
    %   MIAME methods:
    %       MIAME       - Create an MIAME object. 
    %       combine     - Combine two MIAME objects.
    %       isempty     - True for empty MIAME object.
    %
    %   Examples:
    %       expDesc = bioma.data.MIAME('investigator', 'Jane OneName',...
    %                       'lab',      'One Bioinformatics Laboratory',...
    %                       'contact',  'jonename@lab.not.exist',...
    %                       'url',      'www.lab.not.exist',...
    %                       'title',    'Normal vs. Diseased Experiment',...
    %                       'abstract', 'An example for using ExpressionSet',...
    %                       'other',    {'Notes: Created from a text files.'});
    %
    %   See also DATAMATRIX, DATASET, EXPRESSIONSET, EXPTDATA, METADATA.
    
    %   Copyright 2009-2012 The MathWorks, Inc.

    
    % References:
    % http://www.mged.org/Workgroups/MIAME/miame_1.1.html
    
    properties(SetAccess = 'public', GetAccess='public', Hidden=false)
        %INVESTIGATOR Name of the experiment investigator.
        %   The Investigator property is a character array containing the
        %   investigator name.
        %
        %   See also MIAME.
        Investigator = '';
        
        %LABORATORY  The laboratory where the experiment was conducted.
        %   The Laboratory property is a character array containing the
        %   laboratory where the experiment was conducted.
        %
        %   See also MIAME.
        Laboratory = '';
        
        %CONTACT  Contact information for laboratory or investigator.
        %   The Contact property is a character array containing the
        %   contact information of the laboratory and/or the experiment
        %   investigator.
        %
        %   See also MIAME.
        Contact = '';
        
        %URL  A URL for the experiment.
        %   The URL property is a character array containing a URL for the
        %   experiment. 
        %
        %   See also MIAME.
        URL = '';
        
        %TITLE  An experiment title.
        %   The Title property is a character array containing a
        %   single-sentence experiment title. 
        %
        %   See also MIAME.
        Title = '';

        %ABSTRACT  An abstract describing the experiment.
        %   The Abstract property is a character array containing an
        %   abstract describing the experiment. 
        %
        %   See also MIAME.
        Abstract = '';
        
        %EXPTDESIGN A brief description of the experiment design.
        %   The ExptDesign property is a character array containing a
        %   summary describing the experiment overall design. 
        %
        %   See also MIAME.
        ExptDesign = '';

        %ARRAYS  Information about the arrays.
        %   The Arrays property is a cell array containing array-related
        %   information in the experiment. The design information includes,
        %   for example, array design name, platform type or number of
        %   features on the array, etc.
        %
        %   See also MIAME.
        Arrays = {};
        
        %SAMPLES  Cell array of information about the samples.
        %   The Samples property is a cell array containing information
        %   about the samples. Sample information can be, for example, the
        %   sample source, organism, treatment type, compound, labeling
        %   protocol and external control, etc.
        %
        %   See also MIAME.
        Samples = {};

        %HYBRIDIZATION  Cell array of information about the hybridization.
        %   The Hybridization property is a cell array containing
        %   information about the hybridization. Hybridization description
        %   may include, for example, the hybridization protocol including
        %   time, concentration, volume and temperature, etc.
        %
        %   See also MIAME.
        Hybridization = {};
        
        %QUALITYCONTROL Information about the quality controls. 
        %   The QualityControl property is a cell array of information
        %   about the quality control steps taken, for example, replicates
        %   or dye swap, etc.
        %
        %   See also MIAME.
        QualityControl = {};
        
        %PREPROCESSING  Information about the pre-processing steps
        %   The PreProcessing property is a cell array of information about
        %   the pre-processing steps used on the raw data from the
        %   experiment.
        %
        %   See also MIAME.
        PreProcessing = {};
        
        %PUBMEDID  PubMed identifiers of relevant publications. 
        %   The PubMedID property is a character array listing PubMed
        %   identifiers of the papers relevant to the data set from the
        %   experiment.
        %
        %   See also MIAME.
        PubMedID = '';
        
        %OTHER  Other information not covered in other properties.
        %   The Other property is a cell array containing other information
        %   not covered in other properties describing the experiment.
        %
        %   See also MIAME.
        Other = {};
    end
    
    methods
        function obj = MIAME(varargin)
            %MIAME Create an MIAME object.
            %
            %   E = MIAME() creates an object E of MIAME class for storing
            %   information describing a microarray experiment. Namely,
            %   descriptions about experiment design, arrays, samples,
            %   hybridizations, quality controls, and preprocessing steps.
            %   By default all the properties are empty.
            %
            %   E = MIAME(GEOSERSTRUCT) creates an MIAME object E from a
            %   structure containing GEO series data. GEOSERSTRUCT must be
            %   a structure returned by the GETGEODATA function, or its
            %   sub-structure, GEOSERSTRUCT.Header.Series.
            %
            %   E = MIAME(..., 'INVESTIGATOR', NAME) specifies the name of
            %   the experiment investigator. NAME is a character array.
            %
            %   E = MIAME(..., 'LAB', LAB) specifies the laboratory where
            %   the experiment was conducted. LAB is a character array.
            %
            %   E = MIAME(..., 'CONTACT', CONT) specifies the contact
            %   information of the laboratory and/or the experiment
            %   investigator. CONT is a character array.
            %   
            %   E = MIAME(..., 'URL', URL) specifies the URL for the
            %   experiment. URL is a character array.
            
            %   E = MIAME(..., 'TITLE', TLE) specifies the experiment
            %   title. TLE is a character array.
            %
            %   E = MIAME(..., 'ABSTRACT', ABS) specifies the abstract
            %   describing the experiment. ABS is a character array.
            %
            %   E = MIAME(..., 'DESIGN', ED) specifies the brief summary
            %   describing the experiment design. ED is a character
            %   array.
            %
            %   E = MIAME(..., 'ARRAYS', ARR) specifies the array related
            %   information of the experiment. ARR is a cell array.
            %
            %   E = MIAME(..., 'SAMPLES', SAM) specifies the information
            %   about the samples. SAM is a cell array.
            %
            %   E = MIAME(..., 'HYBRID', HYB) specifies the information
            %   about the hybridization. HYB is a cell array.  
            %
            %   E = MIAME(..., 'CONTROLS', QC) specifies the quality
            %   control steps. QC is a cell array.
            %
            %   E = MIAME(..., 'PREPROCESS', PREP) specifies the
            %   information about pre-processing steps used on the raw data
            %   from the experiment. PREP is a cell array.
            %
            %   E = MIAME(..., 'PUBMED', ID) specifies the PubMed
            %   identifications of the papers relevant to the data set. ID
            %   is character array.
            %
            %   E = MIAME(..., 'OTHER', OTH) specifies the other
            %   information not covered in other properties. OTH is a
            %   cell array.
            %
            %   See also MIAME.
            
            %== Process inputs
            if nargin == 0
                return;
            end
            
            %= Check if the first argument is a structure
            if isstruct(varargin{1})
                obj = convertGEOSeries(obj, varargin{1});
                
                obj = createMIAME(obj, varargin{2:end});
            else
                obj = createMIAME(obj, varargin{:});
            end
        end
        
        function c = combine(a, b)
            %COMBINE Combine data from two MIAME objects.
            %
            %   C = COMBINE(A, B) combines data from MIAME object A with
            %   MIAME object B and returns the result in a MIAME object C.
            %   The properties are concatenated together.
            %
            %   See also MIAME
            bioinfochecknargin(nargin, 2, 'MIAME:combine')
            
            if ~isa(a, 'bioma.data.MIAME') || ~isa(b, 'bioma.data.MIAME')
                error(message('bioinfo:MIAME:combine:InvalidInput'))
            end
            %== Find the common variable names in A and B
            try
                c = bioma.data.MIAME;
                c.Investigator = [a.Investigator ' ', b.Investigator];
                c.Laboratory = [a.Laboratory  ' ', b.Laboratory ];
                c.Contact = [a.Contact ' ', b.Contact];
                c.URL = [a.URL ' ', b.URL];
                c.Title = [a.Title ' ', b.Title];
                c.Abstract = [a.Abstract ' ', b.Abstract];
                c.PubMedID = [a.PubMedID ' ', b.PubMedID];
                
                c.ExptDesign = [a.ExptDesign b.ExptDesign];
                c.Arrays = [a.Arrays b.Arrays];
                c.Samples = [a.Samples b.Samples];
                c.Hybridization = [a.Hybridization b.Hybridization];
                c.QualityControl = [a.QualityControl b.QualityControl];
                c.PreProcessing = [a.PreProcessing b.PreProcessing];
                c.Other = [a.Other b.Other];
            catch ME
                bioinfoprivate.bioclsrethrow('MIAME', 'combine', ME)
            end
        end
        
        function t = isempty(M)
            %ISEMPTY True for empty MIAME object.
            %
            %   TF = ISEMPTY(M) returns true (1) if M is an empty MIAME
            %   object and false (0) otherwise. An empty MIAME object has 
            %   empty properties. 
            %
            %   See also MIAME.
            
            t = isempty(M.Investigator) && ...
                isempty(M.Laboratory) && ...
                isempty(M.Contact) && ...
                isempty(M.URL) && ...
                isempty(M.Title) && ...
                isempty(M.Abstract) &&...
                isempty(M.ExptDesign) && ...
                isempty(M.Arrays) && ...
                isempty(M.Samples) && ...
                isempty(M.Hybridization) && ...
                isempty(M.QualityControl) && ...
                isempty(M.PreProcessing) && ...
                isempty(M.PubMedID) && ...
                isempty(M.Other);
        end
    end
    
    methods
        function obj = set.Investigator(obj, name)
            obj.Investigator = setString(name, 'INVESTIGATOR');
        end % Investigator
        
        function obj = set.Laboratory(obj, lab)
           obj.Laboratory = setString(lab, 'LABORATORY');
        end % Laboratory
        
        function obj = set.Contact(obj, cnt)
           obj.Contact = setString(cnt, 'CONTACT');
        end % Contact
        
        function obj = set.URL(obj, url)
           obj.URL = setString(url, 'URL');
        end % URL
        
        function obj = set.Title(obj, title)
           obj.Title = setString(title, 'TITLE');
        end % Title
        
        function obj = set.Abstract(obj, abstract)
           obj.Abstract = setString(abstract, 'ABSTRACT');
        end % Abstract
        
        function obj = set.ExptDesign(obj, design)
           obj.ExptDesign = setString(design, 'EXPDESIGN');
        end % ExptDesign
        
        function obj = set.Arrays(obj, arr)
            obj.Arrays = setCellArray(arr, 'ARRAYS');
        end % Arrays
        
        function obj = set.Samples(obj, sam)
            obj.Samples = setCellArray(sam, 'SAMPLES');
        end % Samples
        
        function obj = set.Hybridization(obj, hyb)
            obj.Hybridization = setCellArray(hyb, 'HYBRIDIZATIONS');
        end % Hybridization
        
        function obj = set.QualityControl(obj, qc)
            obj.QualityControl = setCellArray(qc, 'QUALITYCONTROLS');
        end % QualityControl
        
        function obj = set.PreProcessing(obj, prep)
            obj.PreProcessing = setCellArray(prep, 'PREPROCESSING');
        end % PreProcessing
        
        function obj = set.PubMedID(obj, ids)
            if isnumeric(ids)
                ids = num2str(ids);
            end
            obj.PubMedID = setString(ids, 'PUBMEDID');
        end % PubMedIds
        
        function obj = set.Other(obj, oth)
            obj.Other = setCellArray(oth, 'OTHER');
        end % Other
    end
    
    methods(Hidden=true)
        function [varargout] = subsref(obj,s)
            %SUBSREF Subscripted reference for a MIAME object.
            %
            %   Indexing A(I,J) and cell indexing A{I,J} are not supported.
            %
            %   P = A.PROPERTYNAME returns a MIAME property.     
            switch s(1).type
                case '()'
                    error(message('bioinfo:MIAME:subsref:SubscriptNotSupport'));
                case '{}'
                    error(message('bioinfo:MIAME:subsref:CellSubscript'));
                case '.'
                    try
                        bioma.util.validateMethodProps(obj, s(1).subs);
                        [varargout{1:nargout}] = builtin('subsref', obj, s);
                    catch ME
                        bioinfoprivate.bioclsrethrow('MIAME', 'subsref', ME)
                    end
            end
        end
        
        function [varargout] = subsasgn(obj,s,b)
            %SUBSASGN Subscripted assignment to a MIAME object.
            %
            %   Subscript indexing A(I,J) = B and cell indexing A{I,J} = B
            %   is not supported.
            %
            %   A.PROPERTYNAME = P assigns to a MIAME property.
            
            %== Empty MIAME
            creating = isequal(obj,[]);
            if creating
                obj = bioma.data.MIAME;
            end
            
            switch s(1).type
                case '()'
                    error(message('bioinfo:MIAME:subsasgn:SubscriptNotSupport'));
                case '{}'
                    error(message('bioinfo:MIAME:subsasgn:CellSubscript'));
                case '.'
                    % Return default subsasgn to this object
                    try
                        bioma.util.validateMethodProps(obj, s(1).subs);
                        [varargout{1:nargout}] = builtin('subsasgn',obj,s,b);
                    catch ME
                        bioinfoprivate.bioclsrethrow('MIAME', 'subsasgn', ME)
                    end
            end %switch
        end %subsasgn
        
        function disp(obj)
            fprintf('Experiment Description:\n')
            colPad = repmat(' ', 1, 2);
            fprintf('%sAuthor name: %s\n', colPad, obj.Investigator);
            fprintf('%sLaboratory: %s\n', colPad, obj.Laboratory);
            fprintf('%sContact information: %s\n', colPad, obj.Contact);
            fprintf('%sURL: %s\n', colPad, obj.URL);
            fprintf('%sPubMedIDs: %s\n', colPad, obj.PubMedID);
            %== Abstract
            if isempty(obj.Abstract)
                fprintf('%sNo abstract available.\n', colPad);
            else
                numWords = numel(regexp(obj.Abstract, '\s+'));
                fprintf(['%sAbstract: A %d word abstract is available. '... 
                            'Use the Abstract property.\n'], colPad, numWords);
            end
            
            %== ExptDesign
            if isempty(obj.ExptDesign)
                fprintf('%sNo experiment design summary available.\n', colPad);
            else
                numWords = numel(regexp(obj.ExptDesign, '\s+'));
                fprintf(['%sExperiment Design: A %d word summary is available. '... 
                            'Use the ExptDesign property.\n'], colPad, numWords);
            end
 
            %== Information in cell array
            infoArray = {obj.Arrays, obj.Samples, obj.Hybridization,...
                         obj.QualityControl, obj.PreProcessing};
            infoFlags = arrayfun(@isempty, infoArray);
            
            if any(infoFlags)
               infoNames = fieldnames(obj);
               infoNames = infoNames(7:11);
               infoNames = infoNames(infoFlag);
               infoStr = 'Information is available on: ';
               
               for i=1:numel(infoNames)
                   infoStr = [infoStr infoNames{i} ', ']; %#ok<AGROW>
               end
               infoStr = infoStr(1:end-2); % remove extra comma and space from the end
               fprintf([colPad '\n'])
               fprintf('%s%s\n', colPad, infoStr);
            end
           %== Other
           if ~isempty(obj.Other)
               fprintf('%sOther notes: \n', colPad);
               disp(obj.Other)
           end
        end % disp
        
        function a = subsindex(varargin)  %#ok
            error(message('bioinfo:MIAME:subsindex:UndefinedMethod'));
        end
        function a = ctranspose(varargin)  %#ok
            error(message('bioinfo:MIAME:ctranspose:UndefinedMethod'));
        end
        function a = transpose(varargin)  %#ok
             error(message('bioinfo:MIAME:transpose:UndefinedMethod'));
        end
        function a = permute(varargin)  %#ok
            error(message('bioinfo:MIAME:permute:UndefinedMethod'));
        end
        function a = reshape(varargin)  %#ok
            error(message('bioinfo:MIAME:reshape:UndefinedMethod'));
        end
        function a = cat(varargin)  %#ok
            error(message('bioinfo:MIAME:cat:UndefinedMethod'));
        end
        function a = horzcat(varargin)  %#ok
            error(message('bioinfo:MIAME:horzcat:UndefinedMethod'));
        end
        function a = vertcat(varargin)  %#ok
            error(message('bioinfo:MIAME:vertcat:UndefinedMethod'));
        end
    end
end % End of MIAME class

%-- Helper functions
function obj = createMIAME(obj, varargin)
% Parse input PV pairs.

% Check for the right number of inputs
if rem(nargin-1, 2)== 1
    error(message('bioinfo:MIAME:MIAME:IncorrectNumberOfArguments', mfilename))
end

% Allowed inputs
okargs = {'investigator', 'lab', 'contact', 'url',...
          'title','abstract',...
          'arrays', 'samples', 'hybrid',...
          'controls','preprocess','pubmed', 'other', 'design'};

for j=1:2:nargin-1 
    try 
        [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    catch ME
       bioinfoprivate.bioclsrethrow('MIAME','MIAME', ME); 
    end
    
    switch(k)
        case 1 % investigator
            obj.Investigator = pval;
        case 2 % lab
            obj.Laboratory = pval;
        case 3 % contact
            obj.Contact = pval;
        case 4 % url
            obj.URL = pval;
        case 5 % title
            obj.Title = pval;
        case 6 % abstract
            obj.Abstract = pval;
        case 7 % arrays
            obj.Arrays = pval;
        case 8 % samples
            obj.Samples = pval;
        case 9 % hybrid
            obj.Hybridization = pval;
        case 10 % controls
            obj.QualityControl = pval;
        case 11 % preprocess
            obj.PreProcessing = pval;
        case 12 % pubmed
            obj.PubMedID = pval;
        case 13 % other
            obj.Other = pval;
        case 14 % design
            obj.ExptDesign = pval;
    end
end
end %createMIAME
%-----------
function outStr = setString(inStr, propName)
if iscellstr(inStr)
    inStr = char(inStr);
end

if ischar(inStr)
    outStr = strtrim(inStr(:)');
else
    error(message('bioinfo:MIAME:MIAME:InvalidStringInputFormat', propName))
end
end 
%-------------------
function outVal = setCellArray(inVal, propName)
if ischar(inVal)
    inVal = cellstr(inVal);
elseif ~iscell(inVal)
    inVal = {inVal};
end

if iscell(inVal)
    outVal = inVal(:);
else
    error(message('bioinfo:MIAME:MIAME:InvalidCellInputFormat', propName))
end
end
%---------------------
function obj = convertGEOSeries(obj, gsStruct)
% Convert the GEO series data structure to a MIAME object

%== Check it contains Header field
hasHeader = isfield(gsStruct, 'Header');

if hasHeader 
    gsStruct = gsStruct.Header;
end

if isstruct(gsStruct) && isfield(gsStruct, 'Series')
    gsStruct = gsStruct.Series;
end

fieldNames = {'title'
    'pubmed_id'
    'summary'
    'overall_design'
    'contributor'
    'sample_id'
    'contact_name'
    'platform_id'
    'supplementary_file'
    'contact_institute'
    'web_link'};

if isstruct(gsStruct) && isfield(gsStruct, 'geo_accession')
    if isfield(gsStruct, fieldNames{1})
        obj.Title = gsStruct.(fieldNames{1});
    end
    
    if isfield(gsStruct, fieldNames{2})
        obj.PubMedID = gsStruct.(fieldNames{2});
    end
    
    if isfield(gsStruct, fieldNames{3})
        obj.Abstract = gsStruct.(fieldNames{3});
    end
    
    if isfield(gsStruct, fieldNames{4})
        obj.ExptDesign = gsStruct.(fieldNames{4});
    end
    
    if isfield(gsStruct, fieldNames{5})
        obj.Investigator = gsStruct.(fieldNames{5});
    end
    
    if isfield(gsStruct, fieldNames{6})
        obj.Samples = gsStruct.(fieldNames{6});
    end
    
    if isfield(gsStruct, fieldNames{7})
        obj.Contact = gsStruct.(fieldNames{7});
    end
    
    if isfield(gsStruct, fieldNames{8})
        obj.Arrays = gsStruct.(fieldNames{8});
    end
    
    if isfield(gsStruct, fieldNames{9})
        obj.Other = gsStruct.(fieldNames{9});
    end
    
    if isfield(gsStruct, fieldNames{10})
        obj.Laboratory = gsStruct.(fieldNames{10});
    end
    
    if isfield(gsStruct, fieldNames{11})
        obj.URL = gsStruct.(fieldNames{11});
    end
else
   error(message('bioinfo:MIAME:MIAME:InvalidStruct'))
end

end
