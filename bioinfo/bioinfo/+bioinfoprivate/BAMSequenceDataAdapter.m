classdef BAMSequenceDataAdapter <  bioinfoprivate.SequenceDataAdapter

%   Copyright 2011-2012 The MathWorks, Inc.

    properties  (Access = private) 
        % FIELDBINDING fields names as seen by the BAMIndexedFile object. 
        %   This property is set at construction time. The information in
        %   this property is used to match between the fields requested to
        %   the adapter and the features stored in BAMIndexedFile object.        
        FieldBinding
    end
    
    methods
        
        function obj = BAMSequenceDataAdapter(data,fields,fieldTypes,searchableFields,efficientFields,fieldBinding)
            
             % Make sure efficient access fields are available in index in the
             % BAMIndexedFile object
             h = [];
             for i =1:numel(fields)
                 if efficientFields(i)
                     % Not all the fields are accessible for efficient
                     % access in a BAMIndexedFile (only those pre-indexed by
                     % the undocumented getFeature methods)
                     g = find(strcmp(fieldBinding{i},data.FeatureName));
                     if isempty(g)
                         warning(message('bioinfo:SequenceDataAdapter:NoEfficientAcces', fields{i}))
                     end
                     h = [h g]; %#ok<AGROW>
                 end
             end
            fin = data.FeatureInIndex;
            fin(:,h) = {true};
            data.FeatureInIndex = fin;
            
            obj.InMemory = false;
            obj.FieldNames = fields;
            obj.FieldNamesMap = containers.Map(fields,1:numel(fields));            
            obj.Data = data;
            obj.NumberOfEntries = data.NumEntries;
            obj.EfficientAccessFields = efficientFields;
            obj.StringSearchableFields = searchableFields;            
            obj.FieldTypes = fieldTypes;
            obj.FieldBinding = fieldBinding;
                          
        end

        function details = getAdapterDetails(obj)
            details.FileName = obj.Data.InputFile;
            details.FileFormat = obj.Data.FileFormat;
            details.InMemory = obj.InMemory;
            details.SubsetIndex = obj.Data.SubsetIndex;
            
            % A SubsetIndex may be empty or may have offsets when the
            % Reference property is a subset of the actual references in
            % the BAM file.
            bi = baminfo(obj.Data.InputFile);
            fileRefs = {bi.SequenceDictionary.SequenceName};
            objRefs = obj.Data.Reference;
            if numel(fileRefs)>numel(objRefs)
                % repeat baminfo, slower b/c it scans the whole file
                bi = baminfo(obj.Data.InputFile,'ScanDictionary',true);
                % find index to used references
                ur = seqmatch(objRefs,bi.ScannedDictionary);
                % find index to not used references
                nr = setdiff(1:numel(bi.ScannedDictionary),ur);
                % cum count in file
                cc = [0;cumsum(uint32(bi.ScannedDictionaryCount))];
                % If SubsetIndex is empty, then create dummy indices that
                % cover the whole reference:
                if isempty(details.SubsetIndex)
                    details.SubsetIndex = ...
                      (1:uint32(sum(bi.ScannedDictionaryCount(ur))));
                end
                % Add offsets to account for the missing references
                for i = 1:numel(nr)
                    h = details.SubsetIndex>cc(nr(i));
                    details.SubsetIndex(h) = details.SubsetIndex(h)+uint32(bi.ScannedDictionaryCount(nr(i)));
                end
            end
               
                details.IsSubset = ~isempty(details.SubsetIndex);
            
        end
        
        function  obj = setField(obj,varargin)
            error(message('bioinfo:SequenceDataAdapter:IllegalPropertyAccess'))
        end
        
        function  data = getField(obj,field,s)
            % Indices to features
            h = field2idx(obj,field);
            
            % Process subsetting index:
            if nargin<3 % All entries indexed by the object are requested
                s = []; % empty will indicate 'no' subsetting
            elseif (ischar(s)&&isvector(s)) || iscellstr(s)
                error(message('bioinfo:SequenceDataAdapter:getField:invalidAccess'))
            elseif isnumeric(s)&&isvector(s)
                if any(s<1) || any(s>obj.NumberOfEntries) || any(rem(s,1))
                    error(message('bioinfo:SequenceDataAdapter:getField:invalidNumeric'))
                end
            elseif isempty(s)
                if ~iscell(field)
                    if strcmp(obj.FieldTypes{h},'string')
                        data = {};
                    else
                        data = zeros(0,0,obj.FieldTypes{h});
                    end
                else
                    data = cell(size(field));
                    for i = 1:numel(data)
                        if strcmp(obj.FieldTypes{h(i)},'string')
                            data{i} = {};
                        else
                            data{i} = zeros(0,0,obj.FieldTypes{h(i)});
                        end
                    end
                end
                return
            else
                error(message('bioinfo:SequenceDataAdapter:getField:invalidInput'))
            end
            
            if ~iscell(field) % only one field was requested with a char array
                if isempty(s)
                    data = obj.Data.getFeature(obj.FieldBinding{h});
                else
                    data = obj.Data.getFeature(obj.FieldBinding{h},s);
                end
             else % one or more fields are requested with a cell string
                if isempty(s)
                    data = obj.Data.getFeature(obj.FieldBinding(h));
                else
                    data = obj.Data.getFeature(obj.FieldBinding(h),s);    
                end
             end
        end

        function sobj = getSubsetReference(obj,r)
            sobj = obj;
            sobj.Data = getSubsetReference(sobj.Data,r);
            sobj.NumberOfEntries = sobj.Data.NumEntries;
        end
                
        function sobj = getSubset(obj,idx)
            sobj = obj;
            sobj.Data = getSubset(sobj.Data,idx);
            sobj.NumberOfEntries = sobj.Data.NumEntries;
        end
        
    end

end
