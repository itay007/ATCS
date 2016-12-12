classdef structSequenceDataAdapter <  bioinfoprivate.SequenceDataAdapter

%   Copyright 2011-2012 The MathWorks, Inc.

    properties (Access = private) 

    end
    
    methods
        function obj = structSequenceDataAdapter(data,fields,fieldTypes,searchableFields,efficientFields,fieldBinding)
           
            obj.InMemory = true;
            obj.FieldNames = fields;
            obj.FieldNamesMap = containers.Map(fields,1:numel(fields));
            obj.NumberOfEntries = numel(data);
            obj.EfficientAccessFields = efficientFields;
            obj.StringSearchableFields = searchableFields;
            obj.FieldTypes = fieldTypes;
            obj.Data = struct;
            
            if isempty(data)
                % Create empty container
                for i = 1:numel(fields)
                    if strcmp(fieldTypes{i},'string')
                        obj.Data.(fields{i}) = cell(0,1);
                    else
                        obj.Data.(fields{i}) = zeros(0,1,fieldTypes{i});
                    end
                end
            elseif nargin>5  % fieldBinding is given
                % Copy an array of structures into a scalar structure with
                % vector elements and resolves the field name binding and class
                % type between the input structure and the data container
                for i = 1:numel(fields)
                    if isfield(data,fieldBinding{i})
                        if strcmp(fieldTypes{i},'string')
                            obj.Data.(fields{i}) = {data.(fieldBinding{i})}';
                        elseif isa(data(1).(fieldBinding{i}),fieldTypes{i})
                            obj.Data.(fields{i}) = [data.(fieldBinding{i})]';
                        else
                            switch fieldTypes{i}
                                case 'uint8'
                                    obj.Data.(fields{i}) = uint8([data.(fieldBinding{i})]');
                                case 'uint16'
                                    obj.Data.(fields{i}) = uint16([data.(fieldBinding{i})]');
                                case 'uint32'
                                    obj.Data.(fields{i}) = uint32([data.(fieldBinding{i})]');
                            end
                        end
                    else
                        if strcmp(fieldTypes{i},'string')
                            obj.Data.(fields{i}) = cell(0,1);
                        else
                            obj.Data.(fields{i}) = zeros(0,1,fieldTypes{i});
                        end
                    end
                end
            end
                               
        end
        
        function  obj = setField(obj,field,d)
            % Indices to features
            h = field2idx(obj,field);
            
            if obj.NumberOfEntries>0 && numel(d)>obj.NumberOfEntries
                error(message('bioinfo:SequenceDataAdapter:setField:NotExpandable'))
            end
            
            if strcmp(obj.FieldTypes{h},'string') 
                if ~iscellstr(d)
                    error(message('bioinfo:SequenceDataAdapter:setField:NotCellStr'))
                end
            else
                if iscell(d)
                    d = [d{:}];
                end
                if ~isnumeric(d)
                    error(message('bioinfo:SequenceDataAdapter:setField:NotNumeric'))   
                end
            end
            
            if strcmp(field,'Sequence')
                if obj.NumberOfEntries==0
                    obj.NumberOfEntries = numel(d);
                elseif numel(d)~=obj.NumberOfEntries
                    error(message('bioinfo:SequenceDataAdapter:setField:InvalidSizeSequence',obj.NumberOfEntries,obj.NumberOfEntries))
                end
            elseif all(numel(d)~=[obj.NumberOfEntries,0])
                if strcmp(obj.FieldTypes{h},'string')
                    error(message('bioinfo:SequenceDataAdapter:setField:InvalidSizeCellStr',obj.NumberOfEntries,obj.NumberOfEntries))
                else
                    error(message('bioinfo:SequenceDataAdapter:setField:InvalidSizeNumeric',obj.NumberOfEntries,obj.NumberOfEntries))
                end
            end
            
            if iscellstr(d) || isa(d,obj.FieldTypes{h})
                obj.Data.(field) = d(:);
            else
                    switch obj.FieldTypes{h}
                        case 'uint8'
                            obj.Data.(field) = uint8(d(:));
                        case 'uint16'
                            obj.Data.(field) = uint16(d(:));
                        case 'uint32'
                            obj.Data.(field) = uint32(d(:));
                    end
            end
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
            
            if ~iscell(field)
                if isempty(s)
                    data = obj.Data.(field);
                else
                    data = obj.Data.(field)(s);
                end
            else
                data = cell(size(field));
                for i = 1:numel(data)
                    if isempty(s)
                        data{i} = obj.Data.(field{i});
                    else
                        data{i} = obj.Data.(field{i})(s);
                    end
                end
            end
            
        end
        
        function sobj = getSubsetReference(obj,r)
            idx = getField(obj,'Reference');
            nidx = zeros(obj.NumberOfEntries,1,'uint32');
            for i = 1:numel(r)
                nidx(idx==r(i)) = i;
            end
            obj = setField(obj,'Reference',nidx);
            sobj = getSubset(obj,find(nidx>0));
        end
        
        
        function obj = getSubset(obj,idx)
            for i = 1:numel(obj.FieldNames)
                if ~isempty(obj.Data.(obj.FieldNames{i}))
                    obj.Data.(obj.FieldNames{i}) = obj.Data.(obj.FieldNames{i})(idx);
                end
            end
            obj.NumberOfEntries = numel(idx);
        end
    end
    
end