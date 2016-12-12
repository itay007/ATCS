classdef AnnotationDataAdapter

%   Copyright 2010-2011 The MathWorks, Inc.

    properties
        % DATA data container
        Data
        
        % NUMBEROFENTRIES number of entries seen by the client.
        %   The data in the container may have more entries than what the
        %   client sees; for example, if a subsetting has occurred. With
        %   this property the client can validate a further subscripting
        %   operation. 
        NumberOfEntries
        
        % FIELDNAMES field names as specified by the client.
        %   This property mirrors the same property in the client and it
        %   is passed at construction time (instead of at every method
        %   call). The adapter needs it because the data in its container
        %   may have more fields or the fields may be ordered differently.
        %   With this property the adapter can return the appropriate fields
        %   and in the correct order to the client.
        FieldNames
        
        % STRINGSEARCHABLEFIELDS fields in which the client can search a string.
        %   This property is set at construction time. The generic
        %   implementation of AnnotationDataAdapter does not use it, but
        %   subclass adapters can use this information to create Maps
        %   apriori or to ensure that this fields may be effciently
        %   accessed.
        StringSearchableFields

        % EFFICIENTACCESSFIELDS fields in which the client can frequently access. 
        %   This property is set at construction time. The generic
        %   implementation of AnnotationDataAdapter does not use it, but
        %   subclass adapters can use this information to ensure that
        %   indexed or whole access to individual fields is efficient.
        EfficientAccessFields
        
    end
    
   
    methods (Abstract = true)
        
        % Abstract (virtual) methods must be implemented by the subclasses.
        % This is the minimum set of methods required such that all other
        % generic methods work properly.
        
        data = getField(varargin);
        sobj = getSubset(varargin);
        obj = setField(varargin);
        
    end
    
    methods (Hidden = true)
        
        % These methods persuade the client not to program if-conditional
        % or switch statements based on a particular implementation of the
        % adapter. New adapter implementations shall not add additional
        % functionality other than the defined by the Abstract class
        % AnnotationDataAdapter.
        
        function str = class(obj)
            str = 'bioinfoprivate.AnnotationDataAdapter';
        end
        
        function tf = isa(obj,classname)
            tf = strcmp(classname,'bioinfoprivate.AnnotationDataAdapter');
        end
        
        function disp(obj)
            str = evalc('details(obj)');
            str = regexprep(str,'bioinfoprivate.(\w+)AnnotationDataAdapter','bioinfoprivate.AnnotationDataAdapter');
            disp(str);
        end
        
    end
    
    methods
        
        % Generic methods that depend only on the subclass implementation
        % of the abstract methods (getField and getSubset). Subclasses may
        % implement overloaded methods for these generic methods that
        % access the data more efficiently depending of the particular
        % adapter. 
        
        function tf = isValidField(obj,field)
            if ischar(field) && isrow(field)
                tf = any(strcmp(field, obj.FieldNames));
            elseif iscellstr(field)
                tf = false(size(field));
                for i = 1:numel(field)
                    tf(i) = any(strcmp(field{i}, obj.FieldNames));
                end
            end
        end
        
        function varargout = getFields(obj,fields)
            for i = 1:numel(fields)
                varargout{i} = getField(obj,fields{i});
            end
        end
        
        function data = getIndexedField(obj,idx,field)
            data = getField(obj,field);
            data = data(idx(:));
        end
        
        function obj = setIndexedField(obj,idx,field,d)
            data = getField(obj,field);
            if any(idx>numel(data))
                error(message('bioinfo:AnnotationDataAdapter:AnnotationDataAdapter:NotExpandable'))
            end
            data(idx(:)) = d;
            obj = setField(obj,field,data);
        end
        
        function varargout = getIndexedFields(obj,idx,fields)
            for i = 1:numel(fields)
                data = getField(obj,fields{i});
                varargout{i} = data(idx(:));
            end
        end
        
        function str = getStructArray(obj)
            str = struct([]);
            for i = 1:numel(obj.FieldNames)
                field = obj.FieldNames{i};
                tmp = getField(obj,field);
                if isnumeric(tmp)
                    for j = 1:obj.NumberOfEntries
                        str(j).(field) = tmp(j);
                    end
                else
                    for j = 1:obj.NumberOfEntries
                        str(j).(field) = tmp{j};
                    end
                end
            end
        end
        
        function str = getIndexedStructArray(obj,idx)
            str = struct([]);
            for i = 1:numel(obj.FieldNames)
                field = obj.FieldNames{i};
                tmp = getField(obj,field);
                if isnumeric(tmp)
                    for j = 1:numel(idx)
                        str(j).(field) = tmp(idx(j));
                    end
                else
                    for j = 1:numel(idx)
                        str(j).(field) = tmp{idx(j)};
                    end
                end
            end
        end
        
        function str = getSubStructArray(obj,fields)
            str = struct([]);
            if ischar(fields)
                fields = {fields};
            end
            for i = 1:numel(obj.FieldNames)
                field = obj.FieldNames{i};
                if ismember(field,fields)
                   tmp = getField(obj,field);
                   if isnumeric(tmp)
                       for j = 1:obj.NumberOfEntries
                           str(j).(field) = tmp(j);
                       end
                   else
                       for j = 1:obj.NumberOfEntries
                           str(j).(field) = tmp{j};
                       end
                   end
                end
            end
        end
        
        function str = getIndexedSubStructArray(obj,idx,fields)
            str = struct([]);
            if ischar(fields)
                fields = {fields};
            end            
            for i = 1:numel(obj.FieldNames)
                field = obj.FieldNames{i};
                if ismember(field,fields)
                    tmp = getField(obj,field);
                    if isnumeric(tmp)
                        for j = 1:numel(idx)
                            str(j).(field) = tmp(idx(j));
                        end
                    else
                        for j = 1:numel(idx)
                            str(j).(field) = tmp{idx(j)};
                        end
                    end
                end
            end
        end        
        
        function uni = getDictionaryInSearchableStringField(obj,field)
            uni = unique(getField(obj,field));
        end
          
        function idx = findAnnotations(obj,varargin)
            h = find(strcmp('Range',varargin(1:2:end)));
            if isempty(h)
                idx = findAnnotationsByString(obj,varargin{:});
            elseif numel(varargin)<3
                idx = findAnnotationsByRange(obj,varargin{2}{:});
            else
                idx = findAnnotationsByRangeAndString(obj,varargin{h+1}{:},varargin{setdiff(1:numel(varargin),[h h+1])});
            end
        end
        
        function idx = findAnnotationsByString(obj,varargin)
            lidx = true(obj.NumberOfEntries,1);
            for i = 1:2:numel(varargin)
                if ~isempty(varargin{i+1})
                   tmp =  getField(obj,varargin{i});
                   lidx = lidx & ismember(tmp,varargin{i+1});
                end
            end
            idx = find(lidx);
        end
            
        function idx = findAnnotationsByRangeAndString(obj,X,overlap,fullOverlap,startOnly,varargin)
            idx = findAnnotationsByRange(obj,X,overlap,fullOverlap,startOnly);
            lidx = true(numel(idx),1);
            for i = 1:2:numel(varargin)
                if ~isempty(varargin{i+1})
                   tmp =  getField(obj,varargin{i});
                   lidx = lidx & ismember(tmp(idx),varargin{i+1});
                end
            end
            idx = idx(lidx);
        end
        
        function idx = findAnnotationsByRange(obj,X,overlap,fullOverlap,startOnly)
            start = getField(obj,'Start');
            if isempty(start)
                idx = zeros(0,1);
                return
            end
            stop = getField(obj,'Stop');
            if ~isa(start,'uint32')
                start = uint32(start);
            end
            if ~isa(start,'uint32')
                stop = uint32(stop);
            end
            if ~issorted(start)
                [start,sort_idx] = sort(start);
                stop = stop(sort_idx);
            else
                sort_idx = [];
            end
            idx = bioinfoprivate.getIndexByRangemex(start,stop,uint32(X),uint32(overlap),...
                uint32(0),fullOverlap,startOnly,false);
            if ~isempty(sort_idx)
                idx = sort(sort_idx(idx));
            end          
        end
        
    end 
    
end
        