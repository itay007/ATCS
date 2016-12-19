classdef BIFAnnotationDataAdapter <  bioinfoprivate.AnnotationDataAdapter

%   Copyright 2010-2012 The MathWorks, Inc.
    properties
        FieldBinding
        NumericFields
    end
    
    methods
        
        function obj = BIFAnnotationDataAdapter(data,fields,searchableFields,efficientFields,fieldBinding,numericFields)
            obj.FieldNames = fields;
            obj.Data = data;
            obj.NumberOfEntries = data.NumEntries;
            obj.EfficientAccessFields = efficientFields;
            obj.StringSearchableFields = searchableFields;            
            obj.FieldBinding = containers.Map(fields,fieldBinding);
            obj.NumericFields = numericFields;
        end
        
        function  data = getField(obj,field)
            if ismember(field,obj.NumericFields)
                formatStr = [repmat('%*s',1,obj.FieldBinding(field)-1) '%d%*[^\n]'];
            else
                formatStr = [repmat('%*s',1,obj.FieldBinding(field)-1) '%s%*[^\n]'];
            end
            try
                tmp = textscan(getEntryByIndex(obj.Data,1:obj.NumberOfEntries),formatStr,'delimiter','\t','ReturnOnError',0);
            catch ME 
                error(message('bioinfo:AnnotationDataAdapter:AnnotationDataAdapter:InvalidTableFormatInBIF', field))
            end
            data = tmp{1};
            if (iscellstr(data) && all(cellfun(@isempty,data))) || (isnumeric(data) && all(isnan(data)))
                error(message('bioinfo:AnnotationDataAdapter:AnnotationDataAdapter:InvalidTableFormatInBIF', field))
            end
        end
        
        function sobj = getSubset(obj,idx)
            idx = unique(idx);
            sobj = obj;
            sobj.Data = getSubset(sobj.Data,idx);
            sobj.NumberOfEntries = sobj.Data.NumEntries;
        end
        
        function setField(varargin)
            error(message('bioinfo:AnnotationDataAdapter:AnnotationDataAdapter:IllegalPropertyAccess'))
        end
    end
end
