classdef SequenceAnnotation
%SEQUENCEANNOTATION Abstract class representing a collection of annotations to a sequence.
%
%  A SequenceAnnotation is an abstract class representing a collection of
%  annotations to a sequence. You cannot create instances of this class
%  directly.  You must create a derived class such as GFFAnnotation or
%  GTFAnnotation.
%
%  SequenceAnnotation properties:
%  FieldNames - available data fields.
%  NumEntries - number of annotations in the object.
%
%  SequenceAnnotation methods:
%  getData           - retrieve data from the object as an struct array.
%  getIndex          - get an index to specific annotations. 
%  getSubset         - create another object containing a subset of the annotations.
%  getRange          - get the extent of the annotations in the sequence.
%  getReferenceNames - get the references available in the object.
%  getFeatureNames   - get the types of annotations available in the object.
%
%  See also GFFANNOTATION, GTFANNOTATION.

%   Copyright 2010-2012 The MathWorks, Inc.

%  SequenceAnnotation Hidden methods:
%  getBrowserData - retrieve data for the NGSBrowser.
%  getStructArray - retrieve data as an struct array.
%  getRowInCompactAlignment - calculates a compact alignment.

    properties (Abstract = true, GetAccess='public', SetAccess='private')
        % FIELDNAMES available data fields.
        FieldNames
    end
    properties (Dependent = true, GetAccess='public', SetAccess='private')
        % NUMENTRIES number of annotations in the object.
        NumEntries
    end
    properties (Abstract = true, GetAccess='protected', SetAccess='private')
        StringSearchableFields
        EfficientAccessFields
    end    
    properties (GetAccess='protected', SetAccess='protected')
        Data
    end
 
    methods 
        
%==========================================================================
%   BASIC INFO QUERY METHODS
%==========================================================================        
                  
        function N = get.NumEntries(obj)
            N = obj.Data.NumberOfEntries;
        end
        
        function range = getRange(obj)
%GETRANGE returns the overall range of the annotations.
%
%  GETRANGE(OBJ) returns a 1-by-2 numeric array with the minimum and
%  maximum covered position by all the sequence annotations.
%
%  Example:
%
%  % Create a GFFAnnotation object from a GFF formatted file and query for 
%  % the first and last coordinate of the available features:
%  a = GFFAnnotation('tair8_1.gff')
%  getRange(a)
%
%  See also SEQUENCEANNOTATION.

            range(1,1) = min(getField(obj.Data,'Start'));
            range(1,2) = max(getField(obj.Data,'Stop'));
        end

        function refNames = getReferenceNames(obj)
%GETREFERENCENAMES retrieve unique reference names.
%
%  GETREFERENCENAMES(OBJ) retrieve all the possible strings the Reference
%  field contains.
%
%  Example:
%
%  % Create a GFFAnnotation object from a GFF formatted file and retrieve
%  % the reference names:
%  a = GFFAnnotation('tair8_1.gff')
%  getReferenceNames(a)
%
%  See also SEQUENCEANNOTATION.

            refNames =   getDictionaryInSearchableStringField(obj.Data,'Reference');
        end
        
        function feaNames = getFeatureNames(obj)
%GETFEATURENAMES retrieve unique feature names.
%
%  GETFEATURENAMES(OBJ) retrieve all the possible strings the feature field 
%  contains.   
%
%  Example:
%
%  % Create a GFFAnnotation object from a GFF formatted file and retrieve
%  % the feature names:
%  a = GFFAnnotation('tair8_1.gff')
%  getFeatureNames(a)
%
%  See also SEQUENCEANNOTATION.

            feaNames =   getDictionaryInSearchableStringField(obj.Data,'Feature');
        end
        
%==========================================================================
%   GETDATA METHOD       
%==========================================================================
        function data = getData(obj,varargin)
%GETDATA retrieves data from the object as an struct array.
%
%  GETDATA(OBJ) returns a struct array with data for all the annotations
%  in the object. 
%
%  GETDATA(OBJ,X1,X2) returns a struct array with data for the annotations
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
%  % Create a GFFAnnotation object from a GFF formatted file and extract 
%  % the information within a given range:
%  a = GFFAnnotation('tair8_1.gff')
%  st = getData(a,10000,20000)
%
%  See also SEQUENCEANNOTATION, SEQUENCEANNOTATION/GETFEATURENAMES,
%           SEQUENCEANNOTATION/GETSUBSET, SEQUENCEANNOTATION/GETRANGE, 
%           SEQUENCEANNOTATION/GETREFERENCENAMES. 

            if isempty(varargin)
                data = getStructArray(obj);
            else
                data = getData(getSubset(obj,varargin{:}));
            end
            
        end
        

%==========================================================================
%   GETSUBSET METHOD       
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
%  % Create a GFFAnnotation object from a GFF formatted file and create a 
%  % subset with only the protein features:
%  a = GFFAnnotation('tair8_1.gff')
%  b = getSubset(a,'Feature','protein')
%
%  See also SEQUENCEANNOTATION, SEQUENCEANNOTATION/GETDATA,
%           SEQUENCEANNOTATION/GETFEATURENAMES, SEQUENCEANNOTATION/GETRANGE, 
%           SEQUENCEANNOTATION/GETREFERENCENAMES. 


            if nargin==2 && isnumeric(varargin{1})
                idx = varargin{1};
                if ~isvector(idx) || islogical(idx)|| any(rem(idx,1)~=0) || any(idx<=0) || ...
                        any(idx>obj.Data.NumberOfEntries)
                    error(message('bioinfo:SequenceAnnotation:getSubset:invalidIndex'))
                end
            else
                idx = getIndex(obj,varargin{:});
            end
            objs = obj;
            objs.Data = getSubset(objs.Data,idx);
        end

%==========================================================================
%   GETINDEX METHOD       
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
%  % Create a GFFAnnotation object from a GFF formatted file and get the
%  % start position fro all the protein features:
%  a = GFFAnnotation('tair8_1.gff')
%  idx = getIndex(a,'Feature','protein');
%  prot_starts = a.Start(idx)
%
%  See also SEQUENCEANNOTATION, SEQUENCEANNOTATION/GETDATA,
%           SEQUENCEANNOTATION/GETFEATURENAMES, SEQUENCEANNOTATION/GETRANGE, 
%           SEQUENCEANNOTATION/GETREFERENCENAMES, SEQUENCEANNOTATION/GETSUBSET.

              if nargin==1
                  idx = (1:obj.NumEntries)';
                  return
              end
              
              if nargin==2
                  error(message('bioinfo:SequenceAnnotation:getIndex:IncorrectNumberOfArguments'));
              end

              %=== Input parsing
              [overlap, fullOverlap, startOnly, reference, feature, x1, x2, rangeGiven] = parse_inputs(varargin{:});
                          
              if rangeGiven
                  %=== Error check
                  if ~isnumeric(x1) || ~isscalar(x1) || any(rem(x1,1~=0)) || ...
                     ~isnumeric(x2) || ~isscalar(x2) || any(rem(x2,1~=0)) || ...
                      any(x1(:)<1) || any(x2(:)<x1(:)) || numel(x1)~=numel(x2)
                      error(message('bioinfo:SequenceAnnotation:getIndex:InvalidRange'));
                  end
                  if numel(x1,1) == 0
                     idx = zeros(0,1);
                     return
                  end
              end
              
              if rangeGiven
                  idx = findAnnotations(obj.Data,'Range',{[x1,x2],overlap,fullOverlap,startOnly},'Reference',reference,'Feature',feature);
              else
                  idx = findAnnotations(obj.Data,'Reference',reference,'Feature',feature);
              end
              
        function [overlap, fullOverlap, startOnly, reference, feature, x1, x2, rangeGiven] = parse_inputs(varargin)
            % Parse input PV pairs.
            
            %=== defaults
            fullOverlap = false;
            startOnly = false;
            overlap = 1;
            reference = {};
            feature = {};
            x1 = nan;
            x2 = nan;
            rangeGiven = false;
                        
            %=== check for the right number of inputs
            if rem(nargin,2) ~= 0
                error(message('bioinfo:SequenceAnnotation:getIndex:IncorrectNumberOfArguments'))
            end
            
            %=== parse parameter value pairs
            if nargin > 1
                %=== allowed parameters
                okargs = {'overlap', 'reference', 'feature'};
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
                                error(message('bioinfo:SequenceAnnotation:getIndex:RangeNotGivenInvalidOverlap'))
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
                                error(message('bioinfo:SequenceAnnotation:getIndex:InvalidOverlap'))
                            end
                        case 2  % reference
                            if (ischar(pval) && isrow(pval))
                                reference = {pval};
                            elseif iscellstr(pval)
                                reference = pval(:);
                            else
                                error(message('bioinfo:SequenceAnnotation:getIndex:InvalidReference'))
                            end
                        case 3  % feature
                            if (ischar(pval) && isrow(pval))
                                feature = {pval};
                            elseif iscellstr(pval)
                                feature = pval(:);
                            else
                                error(message('bioinfo:SequenceAnnotation:getIndex:InvalidFeature'))
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

%==================================================================
% SUBSASGN METHOD
%==================================================================
        function obj = subsasgn(obj,s,in)
%SUBSASGN subscripted assignment.
%   OBJ = SUBSASGN(OBJ,S,D) Subscripted assignment stores data into the
%   OBJ.Data container using the string S(1).subs in the same manner data
%   can be set into a MATLAB structure's field when using 'dot' notation.
%   S(1).subs must be a valid string in the FieldNames property. All other
%   types of subscripted assignment are passed to the built-in subsasgn.

            %=== Only intercept dot-subsasgn 
            if s(1).type(1) ~= '.'
                obj = builtin('subsasgn',obj,s,in);
                return
            end

            fieldName = s(1).subs;
            
            % === Is the fieldName stored in the data container?, if not
            %     verify that fieldName is an accessible member for this
            %     class and call then built-in:
            if ~any(strcmp(fieldName, obj.FieldNames))
                mc = metaclass(obj);
                prop = findobj(mc.PropertyList,'Name',fieldName);
                % Is it a non-public property?
                if ~isempty(prop) && ~strcmp(prop.SetAccess,'public')
                    error(message('MATLAB:class:SetProhibited',prop.Name,class(obj)))
                end
                % Is it a non-public method?
                meth = findobj(mc.MethodList,'Name',fieldName);
                if ~isempty(meth) && ~strcmp(meth.Access,'public')
                    error(message('MATLAB:class:MethodRestricted',meth.Name,class(obj)))
                end
                obj = builtin('subsasgn',obj,s,in);
                return
            end
            
            if numel(s) == 1 % set all the property's entries (e.g. obj.Property = d)
                obj.Data = obj.Data.setField(fieldName,in);
            elseif (numel(s)==2) && (s(2).type(1) == '(') && (numel(s(2).subs)==1) %e.g. obj.Property(x) = d
                if isempty(s(2).subs{1})
                    obj.Data = obj.Data.setIndexedField([],fieldName,in);
                elseif islogical(s(2).subs{1})
                    X = find(s(2).subs{1});
                    obj.Data = obj.Data.setIndexedField(X,fieldName,in); %#ok<FNDSB>
                elseif iscell(s(2).subs{1})
                    error(message('bioinfo:SequenceAnnotation:SequenceAnnotation:badCellRef'))
                else
                    if ischar(s(2).subs{1})
                        if strcmp(s(2).subs{1},':')
                            obj.Data = obj.Data.setField(fieldName,in);
                        else
                            %replicates the behavior of subsindex even if it is a valid header name
                            obj.Data = obj.Data.setIndexedField(0+s(2).subs{1},fieldName,in);
                        end
                    else
                        obj.Data = obj.Data.setIndexedField(s(2).subs{1},fieldName,in);
                    end
                end
            elseif (numel(s)>1) && (s(2).type(1) == '{') && (numel(s(2).subs)==1) && numel(s(2).subs{1})==1 %e.g. obj.Property{x} = d (only valid if x is scalar)
                if (numel(s)>2)
                    v = subsref(obj,s(1:2));
                    in = subsasgn(v,s(3:end),in);
                end
                obj.Data = obj.Data.setIndexedField(s(2).subs{1},fieldName,{in});
            else
                error(message('bioinfo:SequenceAnnotation:SequenceAnnotation:IllegalSubscripting', fieldName))
            end
            
        end
        
%==================================================================
% SUBSREF METHOD
%==================================================================
        function varargout = subsref(obj,s)
%SUBSREF subscripted reference.
%   B = SUBSREF(OBJ,S) Subscripted reference access the data stored in the
%   OBJ.Data container using a string S(1).subs in the same manner data can
%   be accessed from a MATLAB structure's field when using 'dot' notation.
%   S(1).subs must be a valid string in the FieldNames property. All other
%   types of subscripted reference are passed to the built-in subsref.
            
            % === Only intercept dot-subsref
            if s(1).type(1) ~= '.' 
                [varargout{1:nargout}] = builtin('subsref',obj,s);
                return
            end
            
            fieldName = s(1).subs;
            
            % === Is the fieldName stored in the data container?, if not
            %     verify that fieldName is an accessible member for this
            %     class and call then built-in:
            if ~any(strcmp(fieldName, obj.FieldNames))
                mc = metaclass(obj);
                prop = findobj(mc.PropertyList,'Name',fieldName);
                % Is it a non-public property?
                if ~isempty(prop) && ~strcmp(prop.GetAccess,'public')
                    error(message('MATLAB:class:GetProhibited',prop.Name,class(obj)))
                end
                % Is it a non-public method?
                meth = findobj(mc.MethodList,'Name',fieldName);
                if ~isempty(meth) && ~strcmp(meth.Access,'public')
                    error(message('MATLAB:class:MethodRestricted',meth.Name,class(obj)))
                end
                [varargout{1:nargout}] = builtin('subsref',obj,s);
                return
            end
            
            if numel(s) == 1 % get all the property's entries (e.g. obj.Property)
                varargout{1} = obj.Data.getField(fieldName);
            elseif (numel(s)==2) && (s(2).type(1) == '(') && (numel(s(2).subs)==1) %e.g. obj.Property(x)
                if isempty(s(2).subs{1})
                    varargout{1} = obj.Data.getIndexedField([],fieldName);
                elseif islogical(s(2).subs{1})
                    X = find(s(2).subs{1});
                    varargout{1} = obj.Data.getIndexedField(X,fieldName); %#ok<FNDSB>
                elseif iscell(s(2).subs{1})
                    error(message('bioinfo:SequenceAnnotation:SequenceAnnotation:badCellRef'))
                else
                    if ischar(s(2).subs{1})
                        if strcmp(s(2).subs{1},':')
                            varargout{1} =  obj.Data.getField(fieldName);
                        else
                            %replicates the behavior of subsindex even if it is a valid header name
                            varargout{1} = obj.Data.getIndexedField(0+s(2).subs{1},fieldName);
                        end
                    else
                        varargout{1} = obj.Data.getIndexedField(s(2).subs{1},fieldName);
                    end
                end
            elseif (numel(s)>1) && (s(2).type(1) == '{') && (numel(s(2).subs)==1) %e.g. obj.Property{x}
                s(2).type(1) = '(';
                varargout = subsref(obj,s(1:2));
                if ~iscell(varargout)
                    error('bioinfo:SequenceAnnotation:SequenceAnnotation:cellRefFromNonCell','Cell contents reference from a non-cell array data field.')
                end
                if numel(varargout)<nargout
                    error(message('bioinfo:SequenceAnnotation:SequenceAnnotation:needMoreRhsOutputs'))
                end
                if (numel(s)>2)  %e.g. obj.Property{x}(y)
                    if numel(varargout)==1
                        varargout{1} = subsref(varargout{1},s(3:end));
                    else
                        error(message('bioinfo:SequenceAnnotation:SequenceAnnotation:badCellRef'))
                    end
                end
            else
                error(message('bioinfo:SequenceAnnotation:SequenceAnnotation:IllegalSubscripting', fieldName))
            end
        end
        
%==========================================================================
%   GETBROWSERDATA METHOD       
%==========================================================================
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
                row = getRowInCompactAlignment(obj);
            else
                idx = getIndex(obj,varargin{:});
                if isempty(idx)
                    fn(1:2*numel(obj.FieldNames)+4) = {cell(0,1)};
                    fn(1:2:end) = [obj.FieldNames 'RowInView' 'GroupIndex'];
                    data = struct(fn{:});
                    return
                end
                data = getStructArray(obj,idx);
                row = getRowInCompactAlignment(obj,idx);
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
                    data(i).GroupIndex = i;
                end
            end
            
        end
        
%==========================================================================
%   GETROWINCOMPACTALIGNMENT METHOD       
%========================================================================== 
        
        function row = getRowInCompactAlignment(obj,idx)
%GETROWINCOMPACTALIGNMENT calculates a compact alignment            
%
%  GETROWINCOMPACTALIGNMENT(OBJ) returns in which row each annotation 
%  in the object should be placed to achieve a compact display.
%
%  GETROWINCOMPACTALIGNMENT(OBJ,IDX) calculates only the row position for
%  the annotations specified by IDX. IDX is a vector with unique positive
%  integers equal or less than the number of entries in the object.

            if nargin < 2
                 start = getField(obj.Data,'Start');
                 stop  = getField(obj.Data,'Stop');
            else
                if ~isnumeric(idx) || ~isvector(idx) || islogical(idx) || ... 
                     any(rem(idx,1)~=0) || any(idx<=0) || any(idx>obj.Data.NumberOfEntries)
                       error(message('bioinfo:SequenceAnnotation:getRowInCompactAlignment:invalidIndex'))
                end
                %=== pull out start and stop vectors from object
                start = getIndexedField(obj.Data,idx,'Start');
                stop  = getIndexedField(obj.Data,idx,'Stop');    
            end  

            %=== sort according to starts in ascending order
			[~, s] = sort(start);
			start = start(s);
			stop = stop(s);
			
			%=== initialize variables
			N = length(start);
            pad = 1; % minimum number of positions between adjacent reads
			row = zeros(N,1);
            
            rowCtr = 1;
            rowEnds = [];
            for i = 1:N
                r = find(rowEnds<start(i),1);
                if isempty(r)
                    r = rowCtr;
                    rowCtr = rowCtr+1;
                end
                row(i) = r;
                rowEnds(r) = stop(i) + pad; %#ok<AGROW>
            end
            
            % Reorder as requested by the input
            row(s) = row;
            
        end
        
%==========================================================================
%   GETSTRUCTARRAY METHOD       
%========================================================================== 

        function str = getStructArray(obj,idxorfields,fields)
%GETSTRUCTARRAY retrieve data as an struct array.
%
%  GETSTRUCTARRAY(OBJ) returns an struct array with the data contained in
%  the SequenceAnnotation object (or derived) for all the available
%  fields and for all the available entries. 
%
%  GETSTRUCTARRAY(OBJ,IDX) returns only the entries specified by the
%  numeric index IDX. IDX is a vector with positive integers equal or less
%  than the number of entries in the object.
%
%  GETSTRUCTARRAY(OBJ,FIELDS) or GETSTRUCTARRAY(OBJ,IDX,FIELDS) returns
%  only the data for the fields specified by FIELDS. FIELDS is a MATLAB
%  string with the field name or a cell string containing multiple field
%  names. 

            if nargin==1
               str = getStructArray(obj.Data);
            elseif nargin==2
               if iscellstr(idxorfields) || (ischar(idxorfields)&&isrow(idxorfields))
                   if all(ismember(idxorfields,obj.FieldNames))
                       str = getSubStructArray(obj.Data,idxorfields);
                   else
                       error(message('bioinfo:SequenceAnnotation:getStructArray:nonExistentField'))
                   end
               elseif isnumeric(idxorfields) && isvector(idxorfields)
                   if all(rem(idxorfields,1)==0) && all(idxorfields>0) && all(idxorfields<=obj.Data.NumberOfEntries)
                       str = getIndexedStructArray(obj.Data,idxorfields);
                   else
                       error(message('bioinfo:SequenceAnnotation:getStructArray:invalidIndex'))
                   end
               else
                   error(message('bioinfo:SequenceAnnotation:getStructArray:invalidSecondInputArgument'))
               end
            else % (nargin == 3)
               if any(rem(idxorfields,1)~=0) || any(idxorfields<=0) || any(idxorfields>obj.Data.NumberOfEntries)
                   error(message('bioinfo:SequenceAnnotation:getStructArray:invalidIndex'))
               end    
               if any(~ismember(fields,obj.FieldNames))
                   error(message('bioinfo:SequenceAnnotation:getStructArray:nonExistentField'))
               end
               str = getIndexedSubStructArray(obj.Data,idxorfields,fields);    
            end
        end    
%==========================================================================
%   PLOT METHOD       
%==========================================================================
        
        function plot(obj)
%PLOT draw the annotations in a figure.
%
%  PLOT(OBJ) opens a new figure and draws the annotations contained in the
%  object OBJ as colored patches. Zooming and panning actions refresh
%  interactively the figure contents. Mouse clicks over annotations display
%  information about the respective annotation.
            
            calX = @(s)  [s.Start;s.Stop;s.Stop;s.Start;s.Start];
            calY = @(s)  [[s.RowInView]-1;[s.RowInView]-1;s.RowInView;s.RowInView;[s.RowInView]-1];
            calC = @(s)  [s.GroupIndex]';
            
            % Calculate initial display range (approx 1000 annotations)
            r = obj.getRange;
            r0 = r(1) + [0 ceil(diff(r)/max(1, obj.NumEntries/1000))];
            
            % Prepare figure
            hf = figure;
            ha = gca;
            hp = patch([1 2 2 1 1]',[0 0 1 1 0]',1); % create dummy patch
            axis([r0 0 50])
            
            s = struct; %init persistent variable
            updatePatchData
            
            % Listeneres for changes of the axes limits
            addlistener(ha,'YLim','PostSet',@(h,e)constrainYlim);
            XLimListener = addlistener(ha,'XLim','PostSet',@(h,e)updatePatchData);
            
            % Prepare data tip (and its callbacks) that moves while the mouse is pressed
            hTip = text(0,0,'','BackgroundColor', [1 1 0.933333],'Color', [0 0 0],...
                'EdgeColor',[0.8 0.8 0.8],'Visible','off','interpreter','none',...
                'fontsize',8,'fontname','courier','horizontal','left');
            set(hf,'WindowButtonUpFcn',@(h,e)deleteDataTip);
            set(hp,'ButtonDownFcn',@(h,e)updateDataTip)
            
            % Avoids updating the patch data while panning
            set(pan,'ActionPreCallback',@(h,e)set(XLimListener,'Enable','off'));
            set(pan,'ActionPostCallback',@(h,e)set(XLimListener,'Enable','on'));
            
            % =========================================================================
            function updatePatchData
                xl = max(1,round(get(ha,'Xlim')));
                s = getBrowserData(obj,xl(1),xl(2));
                set(hp,'XData',calX(s),'YData',calY(s),'FaceVertexCData',calC(s))
            end  % updatePatchData nested function
            % =========================================================================
            function constrainYlim
                set(ha,'YLim',ylim(ha)-min([0,ylim(ha)]))
            end  % constrainYlim nested function
            % =========================================================================
            function updateDataTip
                w = [1 0]*get(ha,'CurrentPoint')*[1 0 0]';
                h = ceil([1 0]*get(ha,'CurrentPoint')*[0 1 0]');
                i = find([s.RowInView]==h);
                i = i(([s(i).Stop]>=w));
                i = i(([s(i).Start]<=w));
                str = cell(numel(i),1);
                for j = 1:numel(i)
                    str{j} = evalc('disp(s(i(j)))');
                end
                set(hTip,'String',sprintf('%s',str{:}),'visible','on')
                set(hf,'WindowButtonMotionFcn',@(h,e) moveDataTip)
                moveDataTip
            end % updateDataTip nested function
            % =========================================================================
            function moveDataTip
                set(hTip,'Position',[1 1]*get(ha,'CurrentPoint')/2)
            end % moveDataTip nested function
            % =========================================================================
            function deleteDataTip
                set(hf,'WindowButtonMotionFcn',[])
                set(hTip,'String','','visible','off')
            end % deleteDataTip nested function
            % =========================================================================
            
        end % plot method
        
    end
end
        
        
    
    
