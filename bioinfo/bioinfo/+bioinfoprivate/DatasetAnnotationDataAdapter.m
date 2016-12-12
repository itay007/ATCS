classdef DatasetAnnotationDataAdapter <  bioinfoprivate.AnnotationDataAdapter
    
%   Copyright 2010-2011 The MathWorks, Inc.
    methods
        
        function obj = DatasetAnnotationDataAdapter(ds,fields,searchableFields,efficientFields)
           obj.FieldNames = fields;
           obj.Data = ds;
           obj.NumberOfEntries = size(ds,1);
           obj.EfficientAccessFields = efficientFields;
           obj.StringSearchableFields = searchableFields;             
        end
        
        function data = getField(obj,field)
            data = obj.Data.(field);
        end
        
        function sobj = getSubset(obj,idx)
            idx = unique(idx);
            sobj = obj;
            sobj.Data = sobj.Data(idx,:);
            sobj.NumberOfEntries = numel(idx);
        end

        function  obj = setField(obj,field,d)
            if isnumeric(obj.Data.(field))
                if ~isnumeric(d)
                    error(message('bioinfo:AnnotationDataAdapter:AnnotationDataAdapter:NotNumeric'))   
                end                
                if ~isa(d,class(obj.Data.(field)))
                    d = cast(d,class(obj.Data.(field)));
                end
            else
                if ~iscellstr(d)
                    error(message('bioinfo:AnnotationDataAdapter:AnnotationDataAdapter:NotCellStr'))
                end                
            end
            if (numel(d)== obj.NumberOfEntries) 
                obj.Data.(field) = d(:);
            elseif numel(d)==1 % scalar expansion
                obj.Data.(field)(:) = d; 
            else
                error(message('bioinfo:AnnotationDataAdapter:AnnotationDataAdapter:InvalidSize'))
            end
        end                
        
        function uni = getDictionaryInSearchableStringField(obj,field)
            uni = unique(getField(obj,field));
        end        
        
    end
end