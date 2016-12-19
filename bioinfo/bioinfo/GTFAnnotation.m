classdef GTFAnnotation <  SequenceAnnotation
%GTFANNOTATION Class representing a collection of GTF annotations.
%
%  A GTFAnnotation is a class representing a collection of annotations to
%  one or more reference sequences conforming with the data available from
%  a GTF formatted file. Annotations may be filtered by feature, by
%  reference or by sequence position. Then data pertaining to one
%  annotation or a subset of annotations can be pulled out into an struct 
%  array.
%
%  GTFAnnotation properties:
%  FieldNames  - available data fields.
%  NumEntries  - number of annotations in the object.
%
%  GTFAnnotation methods:
%  GTFAnnotation      - create a GTFAnnotation object.
%  getData            - retrieve data from the object as an struct array.
%  getIndex           - get an index to specific annotations. 
%  getSubset          - create another object containing a subset of the annotations.
%  getRange           - get the extent of the annotations in the sequence.
%  getReferenceNames  - get the references available in the object.
%  getFeatureNames    - get the types of annotations available in the object.
%  getGeneNames       - get the gene names available in the object.
%  getTranscriptNames - get the transcript names available in the object.
%
%  See also SEQUENCEANNOTATION, GFFANNOTATION.

%   Copyright 2010-2012 The MathWorks, Inc.

%  GTFAnnotation Hidden methods:
%  getBrowserData                  - retrieve data for the NGSBrowser.
%  getStructArray                  - retrieve data as an struct array.
%  getRowInCompactAlignment        - calculate a compact alignment.
%  getRowInGroupedCompactAlignment - calculate a compact alignment.
    
    properties (GetAccess='public', SetAccess='private')
        % FIELDNAMES available data fields.
        FieldNames = {'Reference','Start','Stop','Feature','Gene','Transcript','Source','Score','Strand','Frame','Attributes'};
    end
    properties (GetAccess='protected', SetAccess='private')
        EfficientAccessFields = {'Reference','Feature','Gene','Transcript','Start','Stop'};
        StringSearchableFields = {'Reference','Feature','Gene','Transcript'};
    end
    properties (Dependent)
        Reference
        Start
        Stop
        Feature
        Gene
        Transcript
        Source
        Score
        Strand
        Frame
        Attributes
    end    
    methods
        
%==================================================================
% Constructor
%==================================================================
        function obj = GTFAnnotation(in)
%GTFANNOTATION creates a GTFAnnotation object.
%
%  GTFANNOTATION(FILENAME) creates a GTFAnnotation object from a GTF
%  formatted file FILENAME. 
%
%  Example:
%
%  % Create a GTFAnnotation object from a GTF formatted file.
%  a = GTFAnnotation('hum37_2_1M.gtf')
%
%  See also GFFANNOTATION.
            
            % Check if input is a structure with the proper fields
            if isstruct(in)
                % Check that the input structure has all the required
                % fields for this object (it may have more, that's fine,
                % the StructAnnotationDataAdapter will take care of return
                % the appropriate fields and in the correct order)
                if all(isfield(in,obj.FieldNames))
                    % Construct the StructAnnotationDataAdapter
                    obj.Data =  bioinfoprivate.StructAnnotationDataAdapter(in,obj.FieldNames,obj.EfficientAccessFields,obj.StringSearchableFields);
                else
                    error(message('bioinfo:GTFAnnotation:GTFAnnotation:InvalidStruct'))
                end
                
                % Check if input is a dataset with the proper variables
            elseif isa(in,'dataset')
                if all(ismember(obj.FieldNames,get(in,'VarNames')))
                    obj.Data =  bioinfoprivate.DatasetAnnotationDataAdapter(in,obj.FieldNames,obj.EfficientAccessFields,obj.StringSearchableFields);
                else
                    error(message('bioinfo:GTFAnnotation:GTFAnnotation:InvalidDataset'))
                end
                
                % Check if input is a dataset with sufficient number of fields
                % in the entries
            elseif isa(in,'BioIndexedFile')
                if strcmp(in.FileFormat,'table')
                    fieldBinding = [1 4 5 3 10 11 2 6 7 8 9];
                    numericFields = {'Start','Stop'};
                    obj.Data =  bioinfoprivate.BIFAnnotationDataAdapter(in,obj.FieldNames,obj.EfficientAccessFields,obj.StringSearchableFields,fieldBinding,numericFields);
                else
                    error(message('bioinfo:GTFAnnotation:GTFAnnotation:InvalidFileFormatinBIF'))
                end
                
                % Check if input is a file
            elseif (ischar(in) && (isrow(in) || isempty(in))) && exist(in,'file')
                try
                    fid = fopen(in,'rt');
                    if fid==-1
                        error(message('bioinfo:GTFAnnotation:GTFAnnotation:CannotOpenFile', in))
                    end
                    feat = textscan(fid,'%s%s%s%d%d%s%s%s%s%*[^\n]','delimiter','\t','CommentStyle','#');
                    fclose(fid);
                    
                    % gene_id and transcript_id are mandatory in GTF file
                    % format, but they may be empty (e.g. gene_id "")
                    c = regexp(feat{9},'gene_id\s+"([^"]*)"','tokens','once');
                    feat{10} = [c{:}]';
                    c = regexp(feat{9},'transcript_id\s+"([^"]*)"','tokens','once');
                    feat{11} = [c{:}]';
                    
                    
                    if numel(unique(cellfun(@numel,feat)))>1
                        error(message('bioinfo:GTFAnnotation:GTFAnnotation:InvalidFileFormat'))
                    end
                catch ME 
                    error(message('bioinfo:GTFAnnotation:GTFAnnotation:InvalidFileFormat'))
                end
                fieldBinding = [1 4 5 3 10 11 2 6 7 8 9];
                obj.Data =  bioinfoprivate.CellAnnotationDataAdapter(feat,obj.FieldNames,obj.EfficientAccessFields,obj.StringSearchableFields,fieldBinding);
                
            elseif ischar(in) && (isrow(in) || isempty(in))
                error(message('bioinfo:GTFAnnotation:GTFAnnotation:FileNotFound', in))
            else
                error(message('bioinfo:GTFAnnotation:GTFAnnotation:InvalidInput'))
            end
        end

%==========================================================================
%   BASIC INFO QUERY METHODS ADDED FOR GTF
%==========================================================================        
        
        function geneNames = getGeneNames(obj)
%GETGENENAMES retrieve unique gene names.
%
%  GETGENENAMES retrieve all the possible strings the Gene field contains.
%
%  Example:
%
%  % Create a GTFAnnotation object from a GTF formatted file and retrieve
%  % the gene names:
%  a = GTFAnnotation('hum37_2_1M.gtf')
%  getGeneNames(a)
%
%  See also GTFANNOTATION.
            
            geneNames =   getDictionaryInSearchableStringField(obj.Data,'Gene');
        end
        
        function transcriptNames = getTranscriptNames(obj)
%GETTRANSCRIPTNAMES retrieve unique transcript names.
%
%  GETTRANSCRIPTNAMES retrieve all the possible strings the transcript
%  field contains.
%
%  Example:
%
%  % Create a GTFAnnotation object from a GTF formatted file and retrieve
%  % the trascript names:
%  a = GTFAnnotation('hum37_2_1M.gtf')
%  getTranscriptNames(a)
%
%  See also GTFANNOTATION.
            
            transcriptNames =   getDictionaryInSearchableStringField(obj.Data,'Transcript');
        end

%==========================================================================
%   GETDATA METHOD (overloaded)
%==========================================================================
        
        function data = getData(obj,varargin)
%GETDATA  retrieves data from the object as an struct array.
%
%  GETDATA(OBJ) returns an struct array with data for all the annotations
%  in the object.
%
%  GETDATA(OBJ,X1,X2) returns an struct array with data for the annotations
%  that start before or at X2 and end after or at X1. X1 and X2 are two
%  non-negative integers such that X1 <= X2.
%
%  GETDATA(OBJ,...,'reference',R) or GETDATA(OBJ,X1,X2,...,'reference',R)
%  returns annotations with the field Reference set to R. R is a MATLAB
%  string or a cell string containing multiple references.
%
%  GETDATA(OBJ,...,'feature',F) or GETDATA(OBJ,X1,X2,...,'feature',F)
%  returns annotations with the field Feature set to F. F is a MATLAB
%  string or a cell string containing multiple features.
%
%  GETDATA(OBJ,...,'gene',G) or GETDATA(OBJ,X1,X2,...,'gene',G) returns
%  annotations with the field Gene set to G. G is a MATLAB string or a cell
%  string containing multiple genes.
%
%  GETDATA(OBJ,...,'transcript',T) or GETDATA(OBJ,X1,X2,...,'transcript',T)
%  returns annotations with the field Transcript set to T. T is a MATLAB
%  string or a cell string containing multiple transcripts.
%
%  GETDATA(OBJ,X1,X2,...,'overlap', BP) specifies the minimum number of
%  positions that an annotation must overlap the given range [X1:X2] in 
%  order to be included in the output structure. BP is a number equal to or
%  greater than 1. BP can also contain the string 'full' to indicate that
%  the annotations must be fully contained in the range in order to be
%  included in the output indices. BP can also contain the string 'start'
%  to search for annotations whose start positions lie within the specified
%  range. BP defaults to 1.
%
%  GETDATA(OBJ,IDX) returns a struct array with data for the annotations
%  specified by the numeric index IDX. IDX is a unique column vector with
%  positive integers equal or less than the number of entries in the
%  object.
%
%  Example:
%
%  % Create a GTFAnnotation object from a GTF formatted file and extract 
%  % the information within a given range:
%  a = GTFAnnotation('hum37_2_1M.gtf')
%  st = getData(a,668000,680000)
%
%  See also GTFANNOTATION, GTFANNOTATION/GETGENENAMES, 
%           GTFANNOTATION/GETFEATURENAMES, GTFANNOTATION/GETRANGE, 
%           GTFANNOTATION/GETREFERENCENAMES, GTFANNOTATION/GETSUBSET, 
%           GTFANNOTATION/GETTRANSCRIPTNAMES.
            
            data = getData@SequenceAnnotation(obj,varargin{:});
           
        end

%==========================================================================
%   GETSUBSET METHOD (overloaded)
%==========================================================================
        
        function objs = getSubset(obj,varargin)
%GETSUBSET create another object containing a subset of the annotations.
%
%  GETSUBSET(OBJ,X1,X2) returns an object with a subset of the annotations
%  that start before or at X2 and end after or at X1. X1 and X2 are two
%  non-negative integers such that X1 <= X2.
%
%  GETSUBSET(OBJ,...,'reference',R) or GETSUBSET(OBJ,X1,X2,...,'reference',R)
%  returns an object with a subset of the annotations with the field
%  Reference set to R. R is a MATLAB string or a cell string containing
%  multiple references.
%
%  GETSUBSET(OBJ,...,'feature',F) or GETSUBSET(OBJ,X1,X2,...,'feature',F)
%  returns an object with a subset of the annotations with the field
%  Feature set to F. F is a MATLAB string or a cell string containing
%  multiple features.
%
%  GETSUBSET(OBJ,...,'gene',G) or GETSUBSET(OBJ,X1,X2,...,'gene',G)
%  returns annotations with the field Gene set to G. G is a MATLAB
%  string or a cell string containing multiple genes.
%
%  GETSUBSET(OBJ,...,'transcript',T) or GETSUBSET(OBJ,X1,X2,...,'transcript',T)
%  returns annotations with the field Transcript set to T. T is a MATLAB
%  string or a cell string containing multiple transcripts.
%
%  GETSUBSET(OBJ,X1,X2,...,'overlap', BP) specifies the minimum number of
%  positions that an annotation must overlap the given range [X1:X2] in
%  order to be included in the output object. BP is a number equal to or
%  greater than 1. BP can also contain the string 'full' to indicate that
%  the annotations must be fully contained in the range in order to be
%  included in the output. BP can also contain the string 'start' to search
%  for annotations whose start positions lie within the specified range. BP
%  defaults to 1.
%
%  GETSUBSET(OBJ,IDX) returns an object with a subset of the annotations
%  specified by the numeric index IDX. IDX is a unique column vector with
%  positive integers equal or less than the number of entries in the
%  object.
%
%  Example:
%
%  % Create a GTFAnnotation object from a GTF formatted file and create a 
%  % subset with only the CDS features:
%  a = GTFAnnotation('hum37_2_1M.gtf')
%  b = getSubset(a,'Feature','CDS')
%
%  See also GTFANNOTATION, GTFANNOTATION/GETDATA,
%           GTFANNOTATION/GETGENENAMES, GTFANNOTATION/GETFEATURENAMES,
%           GTFANNOTATION/GETRANGE, GTFANNOTATION/GETREFERENCENAMES,
%           GTFANNOTATION/GETTRANSCRIPTNAMES.

            
           objs = getSubset@SequenceAnnotation(obj,varargin{:});
            
        end        

%==========================================================================
%   GETINDEX METHOD (overloaded)
%==========================================================================
        
        function idx = getIndex(obj,varargin)
%GETINDEX search for specific annotations.
%
%  GETINDEX(OBJ,X1,X2) returns in IDX the index of the annotations that
%  start before or at X2 and end after or at X1. X1 and X2 are two
%  non-negative integers such that X1 <= X2.
%
%  GETINDEX(OBJ,...,'reference',R) or GETINDEX(OBJ,X1,X2,...,'reference',R)
%  returns annotations with the field Reference set to R. R is a MATLAB
%  string or a cell string containing multiple references.
%
%  GETINDEX(OBJ,...,'feature',F) or GETINDEX(OBJ,X1,X2,...,'feature',F)
%  returns annotations with the field Feature set to F. F is a MATLAB
%  string or a cell string containing multiple features.
%
%  GETINDEX(OBJ,...,'gene',G) or GETINDEX(OBJ,X1,X2,...,'gene',G)
%  returns annotations with the field Gene set to G. G is a MATLAB
%  string or a cell string containing multiple genes.
%
%  GETINDEX(OBJ,...,'transcript',T) or GETINDEX(OBJ,X1,X2,...,'transcript',T)
%  returns annotations with the field Transcript set to T. T is a MATLAB
%  string or a cell string containing multiple transcripts.
%
%  GETINDEX (OBJ,X1,X2,...,'overlap', BP) specifies the minimum number of
%  positions that an annotation must overlap the given range [X1:X2] in
%  order to be included in the output indices. BP is a number equal to or
%  greater than 1. BP can also contain the string 'full' to indicate that
%  the annotations must be fully contained in the range in order to be
%  included in the output indices. BP can also contain the string 'start'
%  to search for annotations whose start positions lie within the specified
%  range. BP defaults to 1.
%
%  Example:
%
%  % Create a GTFAnnotation object from a GTF formatted file and retrieve 
%  % the dtsrt podition for all the annotations in a given range:
%  a = GTFAnnotation('hum37_2_1M.gtf')
%  idx = getIndex(a,668000,680000)
%  starts = a.Start(idx)
            
            if nargin==1
                idx = (1:obj.NumEntries)';
                return
            end
            
            if nargin==2
                error(message('bioinfo:GTFAnnotation:getIndex:IncorrectNumberOfArguments'));
            end

            %=== Input parsing
            [overlap, fullOverlap, startOnly, reference, feature, gene, transcript, x1, x2, rangeGiven] = parse_inputs(varargin{:});           
            
            if rangeGiven
                %=== Error check
                if ~isnumeric(x1) || ~isscalar(x1) || any(rem(x1,1~=0)) || ...
                        ~isnumeric(x2) || ~isscalar(x2) || any(rem(x2,1~=0)) || ...
                        any(x1(:)<1) || any(x2(:)<x1(:)) || numel(x1)~=numel(x2)
                    error(message('bioinfo:GTFAnnotation:getIndex:InvalidRange'));
                end
                if numel(x1,1) == 0
                    idx = zeros(0,1);
                    return
                end
            end
          
            if rangeGiven
                idx = findAnnotations(obj.Data,'Range',{[x1,x2],overlap,fullOverlap,startOnly},'Reference',reference,'Feature',feature,'Gene', gene,'Transcript',transcript);
            else
                idx = findAnnotations(obj.Data,'Reference',reference,'Feature',feature,'Gene', gene,'Transcript',transcript);
            end
            
            function [overlap, fullOverlap, startOnly, reference, feature, gene, transcript, x1, x2, rangeGiven] = parse_inputs(varargin)
                % Parse input PV pairs.
                
                %=== defaults
                fullOverlap = false;
                startOnly = false;
                overlap = 1;
                reference = {};
                feature = {};
                gene = {};
                transcript = {};
                x1 = nan;
                x2 = nan;
                rangeGiven = false;                
                
                %=== check for the right number of inputs
                if rem(nargin,2) ~= 0
                    error(message('bioinfo:GTFAnnotation:getIndex:IncorrectNumberOfArguments'))
                end
                
                %=== parse parameter value pairs
                if nargin > 1
                    %=== allowed parameters
                    okargs = {'overlap', 'reference', 'feature', 'gene', 'transcript'};
                    for j = 1:2:nargin-1
                        pname = varargin{j};
                        pval = varargin{j+1};
                        if j==1 && isnumeric(pname)
                            x1 = pname;
                            x2 = pval;
                            rangeGiven = true;
                            k = 0;
                        else                        
                            k = bioinfoprivate.pvpair(pname, pval, okargs,[mfilename,':getIndex']);
                        end
                        switch(k)
                            case 1  % overlap
                                if ~rangeGiven
                                    error(message('bioinfo:GTFAnnotation:getIndex:RangeNotGivenInvalidOverlap'))
                                end
                                if ischar(pval) && strncmpi(pval(:),'full',numel(pval))
                                    fullOverlap = true;
                                    startOnly = false;
                                    overlap = 1;
                                elseif ischar(pval) && strncmpi(pval(:),'start',numel(pval))
                                    fullOverlap = false;
                                    startOnly = true;
                                    overlap = 1;
                                elseif isnumeric(pval) && isscalar(pval) && pval>=1
                                    fullOverlap = false;
                                    startOnly = false;
                                    overlap = pval;
                                else
                                    error(message('bioinfo:GTFAnnotation:getIndex:InvalidOverlap'))
                                end
                            case 2  % reference
                                if (ischar(pval) && isrow(pval))
                                    reference = {pval};
                                elseif iscellstr(pval)
                                    reference = pval(:);
                                else
                                    error(message('bioinfo:GTFAnnotation:getIndex:InvalidReference'))
                                end
                            case 3  % feature
                                if (ischar(pval) && isrow(pval))
                                    feature = {pval};
                                elseif iscellstr(pval)
                                    feature = pval(:);
                                else
                                    error(message('bioinfo:GTFAnnotation:getIndex:InvalidFeature'))
                                end
                            case 4  % gene
                                if (ischar(pval) && (isrow(pval)||isempty(pval)))
                                    gene = {pval};
                                elseif iscellstr(pval)
                                    gene = pval(:);
                                else
                                    error(message('bioinfo:GTFAnnotation:getIndex:InvalidGene'))
                                end
                            case 5  % transcript
                                if (ischar(pval) && isrow(pval))
                                    transcript = {pval};
                                elseif iscellstr(pval)
                                    transcript = pval(:);
                                else
                                    error(message('bioinfo:GTFAnnotation:getIndex:InvalidTranscript'))
                                end
                        end
                    end
                end
                
            end  % function parse_inputs
        end  % function getIndex

    end
    
%==========================================================================
%   HIDDEN METHODS
%==========================================================================
    
    methods (Hidden = true)
        
        
        function data = getBrowserData(obj,varargin)
%GETBROWSERDATA retrieves data for the NGSBrowser.
%
%  GETBROWSERDATA(OBJ) returns an struct array with data for all the
%  annotations in the object. The output structure also contains a field
%  'RowInView' indicating the best row for an annotation in a graphic
%  display and a field 'GroupIndex' indicating how the annotations are
%  grouped in a graphic display.
%
%  GETBROWSERDATA(OBJ,X1,X2) returns an struct array with data for the
%  annotations that start before or at X2 and end after or at X1. X1 and X2
%  are two non-negative integers such that X1 <= X2.
%
%  GETBROWSERDATA(OBJ,...,'reference',R) or GETBROWSERDATA(OBJ,X1,X2,...,'reference',R)
%  returns annotations with the field Reference set to R. R is a MATLAB
%  string or a cell string containing multiple references.
%
%  GETBROWSERDATA(OBJ,...,'feature',F) or GETBROWSERDATA(OBJ,X1,X2,...,'feature',F)
%  returns annotations with the field Feature set to F. F is a MATLAB
%  string or a cell string containing multiple features.
%
%  GETBROWSERDATA(OBJ,...,'gene',G) or GETBROWSERDATA(OBJ,X1,X2,...,'gene',G)
%  returns annotations with the field Gene set to G. G is a MATLAB
%  string or a cell string containing multiple genes.
%
%  GETBROWSERDATA(OBJ,...,'transcript',T) or GETBROWSERDATA(OBJ,X1,X2,...,'transcript',T)
%  returns annotations with the field Transcript set to T. T is a MATLAB
%  string or a cell string containing multiple transcripts.
%
%  GETBROWSERDATA(OBJ,X1,X2,...,'overlap', BP) specifies the minimum number
%  of positions that an annotation must overlap the given range [X1:X2] in
%  order to be included in the output structure. BP is a number equal to or
%  greater than 1. BP can also contain the string 'full' to indicate that
%  the annotations must be fully contained in the range in order to be
%  included in the output indices. BP can also contain the string 'start'
%  to search for annotations whose start positions lie within the specified
%  range. BP defaults to 1.
            
            if isempty(varargin)
                data = getStructArray(obj);
                [row,gi] = getRowInGroupedCompactAlignment(obj,'Transcript');
            else
                idx = getIndex(obj,varargin{:});
                if isempty(idx)
                    fn(1:2*numel(obj.FieldNames)+4) = {cell(0,1)};
                    fn(1:2:end) = [obj.FieldNames 'RowInView' 'GroupIndex'];
                    data = struct(fn{:});
                    return
                end
                data = getStructArray(obj,idx);
                [row,gi] = getRowInGroupedCompactAlignment(obj,'Transcript',idx);
            end
            
            if isempty(row)
                fnames =  [fieldnames(data);{'RowInView';'GroupIndex'}];
                n = 2*numel(fnames);
                data = cell(1,n);
                [data{1:2:n}] = fnames{:};
                for i = 2:2:n
                  data{i} = cell(0,1);
                end
                data = struct(data{:});
            else
                for i = 1:numel(row)
                    data(i).RowInView = row(i);
                    data(i).GroupIndex = gi(i);
                end
            end
            
        end
             
        function [row,gi] = getRowInGroupedCompactAlignment(obj,field,idx)
%GETROWINGROUPEDCOMPACTALIGNMENT calculates a compact alignment
%
%  GETROWINGROUPEDCOMPACTALIGNMENT(OBJ,FIELD) returns in which row each
%  annotation in the object should be placed to achieve a compact display
%  grouping annotations with equal value in the field FIELD along the same
%  row. FIELD is a MATLAB string with either 'Gene' or 'Transcript'.
%
%  GETROWINGROUPEDCOMPACTALIGNMENT(OBJ,FIELD,IDX) calculates only the row
%  position for the annotations specified by IDX. IDX is a vector with
%  unique positive integers equal or less than the number of entries in the
%  object.
%
%  [ROW,GRPIDX] = GETROWINGROUPEDCOMPACTALIGNMENT(...) returns a column
%  vector GRPIDX with the grouping index.
            
            if ~any(strcmp(field,{'Gene','Transcript'}))
                error(message('bioinfo:GTFAnnotation:getRowInGroupedCompactAlignment:invalidField'))
            end
            if nargin < 3
                start = getField(obj.Data,'Start');
                stop  = getField(obj.Data,'Stop');
                gi = grp2idx(getField(obj.Data,field));
            else
                if ~isnumeric(idx) || ~isvector(idx) || islogical(idx) || ...
                        any(rem(idx,1)~=0) || any(idx<=0) || any(idx>obj.Data.NumberOfEntries)
                    error(message('bioinfo:GTFAnnotation:getRowInGroupedCompactAlignment:invalidIndex'))
                end
                %=== pull out start and stop vectors from object
                start = getIndexedField(obj.Data,idx,'Start');
                stop  = getIndexedField(obj.Data,idx,'Stop');
                gi = grp2idx(getIndexedField(obj.Data,idx,field));
            end
            
            %gene_id or transcript_id may contain empty strings, if so
            %grp2idx returns NaNs, we'll treat each annotation without a
            %group as an independent group
            h = isnan(gi);
            if any(h)
                gi(h) = max([gi(~h);0]) + (1:sum(h));
            end
            
            gstart = accumarray(gi,start,[],@min);
            gstop = accumarray(gi,stop,[],@max);
            
            
            %=== sort according to starts in ascending order
            [~, s] = sort(gstart);
            gstart = gstart(s);
            gstop = gstop(s);
            
            %=== initialize variables
            N = length(gstart);
            pad = 1; % minimum number of positions between adjacent reads
            grow = zeros(N,1);
            
            rowCtr = 1;
            rowEnds = [];
            for i = 1:N
                r = find(rowEnds<gstart(i),1);
                if isempty(r)
                    r = rowCtr;
                    rowCtr = rowCtr+1;
                end
                grow(i) = r;
                rowEnds(r) = gstop(i) + pad; %#ok<AGROW>
            end
            
            % Reorder as requested by the input
            grow(s) = grow;
            
            row = grow(gi);
            
        end
        
    end
end

