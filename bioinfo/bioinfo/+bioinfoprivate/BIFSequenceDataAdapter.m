classdef BIFSequenceDataAdapter <  bioinfoprivate.SequenceDataAdapter

%   Copyright 2011-2012 The MathWorks, Inc.

    properties (Access = private) 
        FieldBinding
        % FIELDBINDING fields names as seen by the BioIndexedFile object. 
        %   This property is set at construction time. The information in
        %   this property is used to match between the fields requested to
        %   the adapter and the features stored in BioIndexedFile object.
        
        ColumnPosition
        % COLUMNPOSITION indicates at which column a field is found.
        %   This property is set at construction time. The information in
        %   this property allows to parse directly the text from the source
        %   file when it is formatted as a table preventing a costly call
        %   to the Interpreter stored in the BioIndexedFile object. When a
        %   field name can not be read with a formatted string because the
        %   source is not a table, set the respective ColumnPosition to 0,
        %   this forces to read a field using the Interpreter.
    end
    
    methods
        
        function obj = BIFSequenceDataAdapter(data,fields,fieldTypes,searchableFields,efficientFields,fieldBinding,columnPosition)
            
            % Make sure efficient access fields are loaded in memory in the
            % BioIndexedFile object
            h = [];
            for i =1:numel(fields)
                if efficientFields(i)
                    % Not all the fields are accessible for efficient access
                    % in a BioIndexedFile (only those pre-indexed by the
                    % undocumented getFeature methods)
                    g = find(strcmp(fieldBinding{i},data.FeatureName));
                    if isempty(g)
                        warning(message('bioinfo:SequenceDataAdapter:NoEfficientAcces', fields{i}))
                    end
                    h = [h g]; %#ok<AGROW>
                end
            end
            fin = data.FeatureInMemory;
            fin(h) = {true};
            data.FeatureInMemory = fin;
            
            obj.InMemory = false;
            obj.FieldNames = fields;
            obj.FieldNamesMap = containers.Map(fields,1:numel(fields));
            obj.Data = data;
            obj.NumberOfEntries = data.NumEntries;
            obj.EfficientAccessFields = efficientFields;
            obj.StringSearchableFields = searchableFields;            
            obj.FieldTypes = fieldTypes;
            obj.FieldBinding = fieldBinding;
            obj.ColumnPosition = columnPosition;
                           
        end
        
        function details = getAdapterDetails(obj)
            details.FileName = obj.Data.InputFile;
            details.FileFormat = obj.Data.FileFormat;
            details.InMemory = obj.InMemory;
            details.SubsetIndex = obj.Data.SubsetIndex;
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
            
            % Features that shall be accessed by column
            hcl = ~strcmp(obj.FieldTypes(:),'string')'; % guaranty a row
          
            if ~iscell(field) % only one field was requested with a char array
                if hcl(h) % get whole numeric column
                    if isempty(s)
                        data = obj.Data.getFeature(obj.FieldBinding{h});
                    else
                        data = obj.Data.getFeature(obj.FieldBinding{h},s);
                    end
                    if ~isa(data,obj.FieldTypes{h})
                        switch obj.FieldTypes{h}
                            case 'uint8'
                                data = uint8(data);
                            case 'uint16'
                                data = uint16(data);
                            case 'uint32'
                                data = uint32(data);
                        end
                    end
                else % read string data from source file
                    if isempty(s)
                        s = 1:obj.NumberOfEntries;
                    end
                    if obj.ColumnPosition(h) == 0
                        data = cell(numel(s),1);
                        block_size = 100000;
                        if numel(s)>block_size
                            for i = 1:block_size:numel(s)
                                j = min(i+block_size-1,numel(s));
                                str = read(obj.Data,s(i:j));
                                data(i:j) = {str.(obj.FieldBinding{h})};
                            end
                        else
                            str = read(obj.Data,s);
                            data =  {str.(obj.FieldBinding{h})}';
                        end
                    else
                        % prepare format string to read with textscan:
                        f = repmat(('%*s')',1,max(obj.ColumnPosition));
                        %f(3,nonzeros(obj.ColumnPosition(hcl))) = 'd';
                        f(2,obj.ColumnPosition(h)) = ' ';
                        f = [strrep(f(1:obj.ColumnPosition(h)*3),' ','') '%*[^\n]'];
                        % read data from file by blocks
                        data = cell(numel(s),1);
                        block_size = 100000;
                        if numel(s)>block_size
                            for i = 1:block_size:numel(s)
                                j = min(i+block_size-1,numel(s));
                                str = textscan(obj.Data.getEntryByIndex(s(i:j)),f);
                                data(i:j) = str{1};
                            end
                        else
                            str = textscan(obj.Data.getEntryByIndex(s),f);
                            data = str{1};
                        end
                    end
                end
            else % one or more fields are requested with a cell string
                data = cell(size(field));
                hp = 1:numel(h); % keeps track of desired position in output
                % reorder to the same order as features appear in the source file
                hc = obj.ColumnPosition(h);
                [~,si] = sort(hc);
                h = h(si);
                hp = hp(si);
                
                % Get features that are in memory                
                for i = find(hcl(h))
                    if isempty(s)
                        data{hp(i)} = obj.Data.getFeature(obj.FieldBinding{h(i)});
                    else
                        data{hp(i)} = obj.Data.getFeature(obj.FieldBinding{h(i)},s);
                    end
                    if ~isa(data{hp(i)},obj.FieldTypes{h(i)})
                        switch obj.FieldTypes{h(i)}
                            case 'uint8'
                                data{hp(i)} = uint8(data{hp(i)});
                            case 'uint16'
                                data{hp(i)} = uint16(data{hp(i)});
                            case 'uint32'
                                data{hp(i)} = uint32(data{hp(i)});
                        end
                    end
                end
                % Get features that are in memory
                i = find(~hcl(h));
                if numel(i)
                    if isempty(s)
                        s = 1:obj.NumberOfEntries;
                    end
                    if any(obj.ColumnPosition(h) == 0)
                        data(hp(i)) = {cell(numel(s),1)};
                        block_size = 100000;
                        if numel(s)>block_size
                            for j = 1:block_size:numel(s)
                                k = min(j+block_size-1,numel(s));
                                str = read(obj.Data,s(j:k));
                                for l=1:numel(i)
                                   data{hp(i(l))}(j:k) = {str.(obj.FieldBinding{h(i(l))})};
                                end
                            end
                        else
                            str = read(obj.Data,s);
                            for l=1:numel(i)
                                data{hp(i(l))} = {str.(obj.FieldBinding{h(i(l))})}';
                            end
                        end                        
                    else
                        % prepare format string to read with textscan:
                        f = repmat(('%*s')',1,max(obj.ColumnPosition));
                        %f(3,nonzeros(obj.ColumnPosition(hcl))) = 'd';
                        f(2,obj.ColumnPosition(h(i))) = ' ';
                        f = [strrep(f(1:max(obj.ColumnPosition(h(i)))*3),' ','') '%*[^\n]'];
                        
                        % read data from file by blocks
                        block_size = 100000;
                        if numel(s)>block_size
                            data(hp(i)) = {cell(numel(s),1)};
                            for j = 1:block_size:numel(s)
                                k = min(j+block_size-1,numel(s));
                                str = textscan(obj.Data.getEntryByIndex(s(j:k)),f);
                                for l=1:numel(i)
                                    data{hp(i(l))}(j:k) = str{l};
                                end
                            end
                        else
                            str = textscan(obj.Data.getEntryByIndex(s),f);
                            data(hp(i))=str;
                        end
                    end
                end
            end
        end
        
        function sobj = getSubsetReference(obj,r)
            idx = getField(obj,'Reference');
            for i = 1:numel(r)
                idx(idx==r(i)) = 0;
            end
            sobj = getSubset(obj,find(idx==0));
        end
        
        function sobj = getSubset(obj,idx)
            sobj = obj;
            sobj.Data = getSubset(sobj.Data,idx);
            sobj.NumberOfEntries = sobj.Data.NumEntries;
        end
        
    end
end
