classdef StructAnnotationDataAdapter <  bioinfoprivate.AnnotationDataAdapter
    
%   Copyright 2010-2011 The MathWorks, Inc.
    properties
        NumericFields
    end
    
    methods
        
        function obj = StructAnnotationDataAdapter(data,fields,searchableFields,efficientFields,numericFields)
            infieldnames = fieldnames(data);
            if numel(infieldnames)~=numel(fields)
                data = rmfield(data,setdiff(infieldnames,fields));
                infieldnames = fieldnames(data);
            end
            if all(strcmp(infieldnames(:),fields(:)))
                obj.Data = data;
            else
                obj.Data = orderfields(data,fields);
            end
            obj.NumberOfEntries = numel(data);
            obj.FieldNames = fields;
            obj.EfficientAccessFields = efficientFields;
            obj.StringSearchableFields = searchableFields;
            if nargin==4
                numericFields = false(1,numel(fields));
                for i = 1:numel(fields)
                    field = fields{i};
                    if all(cellfun(@(x) isnumeric(x)&isscalar(x),{obj.Data.(field)}))
                       numericFields(i) = true;
                    end
                end
                obj.NumericFields = fields(numericFields);
            else
                obj.NumericFields = numericFields;
            end
        end
        
        function data = getField(obj,field)
            if ismember(field,obj.NumericFields)
                data = [obj.Data(:).(field)]';
            else
                data = {obj.Data(:).(field)}';
            end
        end
        
        function sobj = getSubset(obj,idx)
            idx = unique(idx);
            sobj = obj;
            sobj.Data = sobj.Data(idx);
            sobj.NumberOfEntries = numel(idx);
        end
        
        function  obj = setField(obj,field,d)
            if isnumeric(obj.Data(1).(field))
                if ~isnumeric(d)
                    error(message('bioinfo:AnnotationDataAdapter:AnnotationDataAdapter:NotNumeric'))   
                end                
                if ~isa(d,class(obj.Data(1).(field)))
                    d = cast(d,class(obj.Data(1).(field)));
                end                
            else
                if ~iscellstr(d)
                    error(message('bioinfo:AnnotationDataAdapter:AnnotationDataAdapter:NotCellStr'))
                end                
            end
            if (numel(d)== obj.NumberOfEntries)
                for i = 1:numel(d)
                    obj.Data(i).(field) = d(i);
                end
            elseif numel(d)==1 % scalar expansion
                for i = 1:numel(d)
                    obj.Data(i).(field) = d;
                end
            else
                error(message('bioinfo:AnnotationDataAdapter:AnnotationDataAdapter:InvalidSize'))
            end
        end
        
        function str = getStructArray(obj)
            str = obj.Data;
        end
        
        function str = getIndexedStructArray(obj,idx)
            str = obj.Data(idx);
        end
        
    end
    
end