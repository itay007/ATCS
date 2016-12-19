classdef SequenceDataAdapter

%   Copyright 2011-2012 The MathWorks, Inc.

    properties (Access = protected)
        % DATA data container
        Data
    end

    properties (SetAccess = protected, GetAccess = public)        
        % INMEMORY indicates if data is loaded in memory.
        %   Used by the client to make proper guess about accessing the 
        %   data or not.
        InMemory
    end
    
    properties (SetAccess = protected, GetAccess = public)        
        % NUMBEROFENTRIES number of entries seen by the client.
        %   The data in the container may have more entries than what the
        %   client sees; for example, if a subsetting has occurred. With
        %   this property the client can validate a further subscripting
        %   operation. 
        NumberOfEntries
    end
        
    properties (SetAccess = protected, GetAccess = protected)        
        % FIELDNAMES field names as specified by the client.
        %   This property mirrors the same property in the client and it
        %   is passed at construction time (instead of at every method
        %   call). The adapter needs it because the data in its container
        %   may have more fields or the fields may be ordered differently.
        %   With this property the adapter can return the appropriate fields
        %   and in the correct order to the client.
        FieldNames
    end
       
    properties (SetAccess = protected, GetAccess = private)        
        % FIELDNAMESMAP map for the field names.
        %   This property the order in which the property names are stored
        %   in the adapter. Other properties seen by the adapter have the
        %   same respective order. 
        FieldNamesMap
    end
    
    properties (Access = protected)        
        % FIELDTYPES field type as expected by the client.
        %   This property is set at construction time. It indicated the
        %   type of data expected by the client for ecery feature. If
        %   casting is necessary the adaptor does it.
        FieldTypes        
    end
        
    properties (Access = protected)        
        % STRINGSEARCHABLEFIELDS fields in which the client can search a string.
        %   This property is set at construction time. The generic
        %   implementation of SequenceDataAdapter does not use it, but
        %   subclass adapters can use this information to create Maps
        %   apriori or to ensure that this fields may be effciently
        %   accessed.
        StringSearchableFields
    end

    properties (Access = protected)        
        % EFFICIENTACCESSFIELDS fields in which the client can frequently access. 
        %   This property is set at construction time. The generic
        %   implementation of SequenceDataAdapter does not use it, but
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
        sobj = getSubsetReference(varargin);
        obj = setField(varargin);
       
    end
    
    methods (Hidden = true)
        
        % These methods persuade the client not to program if-conditional
        % or switch statements based on a particular implementation of the
        % adapter. New adapter implementations shall not add additional
        % functionality other than the defined by the Abstract class
        % SequenceDataAdapter.
        
        function str = class(obj)
            str = 'bioinfoprivate.SequenceDataAdapter';
        end
        
        function tf = isa(obj,classname)
            tf = strcmp(classname,'bioinfoprivate.SequenceDataAdapter');
        end
        
        function disp(obj)
            str = evalc('details(obj)');
            str = regexprep(str,'bioinfoprivate.(\w+)SequenceDataAdapter','bioinfoprivate.SequenceDataAdapter');
            disp(str);
        end
        
        function details = getAdapterDetails(obj)
            details.FileName = [];
            details.FileFormat = [];
            details.InMemory = obj.InMemory;
            details.IsSubset = false;
            details.SubsetIndex = [];
        end
        
    end
    
    methods (Access = protected)
        
        function idx = field2idx(obj,field)
            if ischar(field) && isrow(field)
                try
                    idx = obj.FieldNamesMap(field);
                catch ME
                    if strcmp(ME.identifier,'MATLAB:Containers:Map:NoKey')
                        error(message('bioinfo:SequenceDataAdapter:field2idx:UnavailableField', field))
                    else
                        rethrow(ME)
                    end
                end
            elseif iscellstr(field)
                idx = zeros(size(field));
                for i = 1:numel(field)
                    try
                        idx(i) = obj.FieldNamesMap(field{i});
                    catch ME
                        if strcmp(ME.identifier,'MATLAB:Containers:Map:NoKey')
                            error(message('bioinfo:SequenceDataAdapter:field2idx:UnavailableField', field{ i }))
                        else
                            rethrow(ME)
                        end
                    end
                end
            else
                error(message('bioinfo:SequenceDataAdapter:field2idx:InvalidField'))
            end
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
    end 
    
end
