classdef BioSeq < internal.matlab.variableeditor.VariableEditorPropertyProvider
	%BIOSEQ class representing a collection of biological sequences.
	%
	%  A BioSeq object is a generic representation of biological sequence
	%  data. It consists of a set of sequences and headers.
    %
    %  NOTE: The BioSeq is the base class for BioMap and BioRead. Usually,
    %  it is preferred to create intances of any of its subclasses BioMap
    %  or BioRead instead of using the BioSeq constructor directly. The
    %  interface for the BioSeq constructor is undocumented and may change
    %  in a future release of the Bioinformatics Toolbox.
	%
	%  BioSeq properties:
	%  Sequence - sequence data.
	%  Header   - header information for each sequence.
	%  NSeqs    - number of sequences represented in the object.
	%  Name     - name of the object.
	%
	%  BioSeq specialized property access:
	%  get            - get multiple properties with a single statment.
	%  getHeader      - subscripted retrieval of the 'Header' property.
	%  getSequence    - subscripted retrieval of the 'Sequence' property.
	%  getSubsequence - subscripted retrieval of partial values of the 'Sequence' property.
	%  set            - set multiple properties with a single statment.
	%  setHeader      - subscripted assignment of the 'Header' property.
	%  setSequence    - subscripted assignment of the 'Sequence' property.
	%  setSubsequence - subscripted assignment of partial values of the 'Sequence' property.
    %    
    %  BioSeq methods:
	%  BioSeq         - create a BioSeq object.
	%  combine        - combine two objects of BioSeq class.
	%  getSubset      - retrieve a subset of elements.
	%  setSubset      - set the values of a subset of elements.
	%
	%  See also BIOMAP, BIOREAD.
	
	%   Copyright 2009-2012 The MathWorks, Inc.
	
	properties (Dependent)
		%SEQUENCE sequence data represented in the object.
		%     The 'Sequence' property is a cell array of strings
		%     containing the sequences represented in the object.
		%
		%     See also BIOSEQ.
        Sequence
    end

	properties (Dependent)		
		%HEADER header information for the sequences in the object.
		%     The 'Header' property is a cell array of strings
		%     containing the header information for each sequence
		%     represented in the object.
		%
		%     See also BIOSEQ.
		Header
    end
    
    properties (Dependent, SetAccess = private)		
		%NSEQS number of sequences in the object.
		%     The 'NSeqs' property is an integer representing the number of
		%     sequences in the object.
		%
		%     See also BIOSEQ.
		NSeqs
	end
	
	properties
		%NAME name of the object.
		%     The 'Name' property is a string specifying the name of a
		%     BioSeq (or derived) object. Default is a empty string.
		%
		%     See also BIOSEQ.
		Name = '';
    end
	
	properties (Access = protected, Hidden = true)
		Index = [];
	end
	
	%======================================================================
	% PUBLIC METHODS
	%======================================================================
	
	methods
		
		%==================================================================
		% Constructor
		%==================================================================
		
		function obj = BioSeq(in, varargin)
			%BIOSEQ create a BioSeq object.
			%
            %  B = BIOSEQ(FILENAME) creates a BioSeq object B from the content of a 
            %  FASTA/FASTQ/SAM formatted file FILENAME. Data is kept in the file and 
            %  the BioSeq object accesses it using an auxiliary index file with the 
            %  extension .IDX. If the index file is not present BioSeq will create one
            %  in the same directory where the source data is.            
            %
			%  B = BIOSEQ(STRUCT) creates a BioSeq object B from the content of a
            %  MATLAB structure STRUCT, such as the structure returned by the 
            %  FASTAREAD, FASTQREAD and SAMREAD functions.
            %            
			%  B = BIOSEQ(SEQS) creates a BioSeq object B from a cell array of strings 
			%  SEQS representing the values of the BioSeq property 'Sequence'.
			%
			%  B = BIOSEQ(SEQS, HEADS) creates a BioSeq object B and sets the BioSeq 
            %  property 'Header' to HEADS. HEADS is a cell array of strings.
            %
            %  B = BIOSEQ(..., 'InMemory', true) forces data to be loaded into memory 
            %  when using an auxiliary index file. InMemory defaults to false. InMemory 
            %  is ignored when the first input is not a file name.  
            %
            %  B = BIOSEQ(..., 'IndexDir', INDEXDIR) sets the path to the directory 
            %  where the index files can be found or created. 
            %
            %  B = BIOSEQ(..., PROPERTY_NAME , VALUE) specifies optional parameter 
            %  name/value pairs to further modify the object properties when allowed. 
            %  PROPERTY_NAME is a string.
			%
			%
			%  Examples:
			%
			%  % Create a BioSeq object from a FASTA file.
			%  obj = BioSeq('primatesaligned.fa')
			%
			%  % Create a BioSeq object from variables in the workspace.
			%  seqs = {randseq(10); randseq(15); randseq(20)};
			%  heads = {'H1'; 'H2'; 'H3'};
			%  obj = BioSeq(seqs, heads)
			%
			%  % Create a BioSeq object from a valid structure.
			%  str = struct('Header', heads, 'Sequence', seqs);
			%  obj = BioSeq(str)
			%
			%  % Create a BioSeq object from a FASTA file and reset
			%  % the 'Header' property.
			%  obj = BioSeq('primatesaligned.fa', 'InMemory', true, 'Header', {})
			%
			%
			%  See also BIOMAP/BIOMAP, BIOREAD/BIOREAD, BIOSEQ.
			
            if (nargin == 0)
                in = [];
            end

            if iscellstr(in)
                if numel(varargin)>0 && iscellstr(varargin{1})
                    %  BIOSEQ(SEQS,HEADS,...)
                    varargin = {'Sequence', in, 'Header', varargin{1:end}};
                    in = [];
                else
                   %  BIOSEQ(SEQS,...)
                   varargin = {'Sequence', in, varargin{1:end}};
                   in = [];
                end
            end

            %=== Check if the (new) first input argument may be a PVP
            %    (building from only PV pairs with properties only is
            %    undocummented but required from some object methods)
            if  ~isempty(in) && ~isstruct(in) &&  ~((ischar(in) && (isrow(in) || isempty(in)))&&any(in=='.'))
                varargin = {in, varargin{:}}; %#ok<CCAT>
                in = []; % Indicate later to build just from PVPs
            end           
            
            %=== Check if within the PVPs, inMemory is set
            %    from PVPs
            removeVarargin = [];
            inMemory = false;
            for i = 1:2:numel(varargin)-1
                pname = varargin{i};
                if ischar(pname)
                    if strncmpi('inmemory',pname,max(1,numel(pname)))
                        removeVarargin = [removeVarargin i i+1]; %#ok<AGROW>
                        inMemory = bioinfoprivate.opttf(varargin{i+1}, 'InMemory', ['BioSeq:' mfilename]); 
                    end
                end
            end            
            
            %=== Check if within the PVPs, a directory for the index is given:
            indexDir = [];
            for i = 1:2:numel(varargin)-1
                pname = varargin{i};
                if ischar(pname)
                    if strncmpi('indexdir',pname,max(1,numel(pname)))
                        removeVarargin = [removeVarargin i i+1]; %#ok<AGROW>
                        indexDir = varargin{i+1}; 
                        if ~ischar(indexDir) || ~isrow(indexDir)
                             error(message('bioinfo:BioSeq:BioSeq:InvalidIndexDir'))
                        end
                        if ~isdir(indexDir)
                            error(message('bioinfo:BioSeq:BioSeq:UnknownIndexDir', indexDir))
                        end
                    end
                end
            end 
            
            %=== Remove PVPs already processed from varargin, all remaining
            %    PVPs shall be only settable object properties
            if numel(removeVarargin)
                varargin(removeVarargin) = [];
            end              
            
            %=== Build BioSeq from FileName
            if ischar(in) && (isrow(in) || isempty(in))
                if ~exist(in,'file')
                    error(message('bioinfo:BioSeq:BioSeq:FileNotFound', in))
                end
                % if file is not in the current 
                % detect type of file
                [pa,~,fileExt] = fileparts(in);
                if isempty(pa) && ~strcmpi(fileparts(which(in)),pwd)
                    % it may be possible that the user is expecting to find
                    % a file in the Matlab path, if so, the path is
                    % appended to the filename when it is different from the
                    % current directory
                    in = which(in);
                end
                if strncmpi('.sam',fileExt,4)
                    if inMemory
                        in = samread(in,'TAGS',false);
                        obj.Index = bioinfoprivate.structSequenceDataAdapter( ...
                           in,...  % data container
                           {'Sequence','Header' },...
                           {'string'  ,'string' },...
                           [false     ,false    ],...% searchable fields
                           [true      ,true     ],...% efficient fields
                           {'Sequence','QueryName' });  % field binding
                    else
                        if isempty(indexDir)
                            in = BioIndexedFile('sam',in,'indexedbykeys',false,'verbose',false);
                        else
                            in = BioIndexedFile('sam',in,indexDir,'indexedbykeys',false,'verbose',false);
                        end
                            
                        % Validate that what was indexed is a SAM formatted file
                        try
                           samread(getEntryByIndex(in,1));
                        catch ME 
                            try
                                delete(in.IndexFile)
                            catch ME2 %#ok<NASGU>
                                %could not delete, just continue with error
                            end                            
                            error(message('bioinfo:BioSeq:BioSeq:InvalidSAMFormat', in.InputFile))
                        end
                        
                        obj.Index = bioinfoprivate.BIFSequenceDataAdapter( ...
                            in,...  % data container
                            {'Sequence','Header'},...
                            {'string','string'},...
                            [false,false],...  % searchable fields
                            [false,false],...    % efficient fields
                            {'Sequence','QueryName'},...  % field binding
                            [10,1]); % columnPosition
                    end                    
                elseif (strncmpi('.fasta',fileExt,6)||strncmpi('.fa',fileExt,3))
                    if inMemory
                        in = fastaread(in, 'trimheaders', true);
                        obj.Index = bioinfoprivate.structSequenceDataAdapter( ...
                           in,...  % data container
                           {'Sequence','Header' },...
                           {'string'  ,'string' },...
                           [false     ,false    ],...% searchable fields
                           [true      ,true     ],...% efficient fields
                           {'Sequence','Header' });  % field binding
                    else
                        if isempty(indexDir)
                            in = BioIndexedFile('fasta',in,'indexedbykeys',false,'verbose',false);
                        else
                            in = BioIndexedFile('fasta',in,indexDir,'indexedbykeys',false,'verbose',false);
                        end
                        
                        % Validate that what was indexed is a FASTA formatted file
                        try
                           fastaread(getEntryByIndex(in,1));
                        catch ME 
                            try
                                delete(in.IndexFile)
                            catch ME2 %#ok<NASGU>
                                %could not delete, just continue with error
                            end                                      
                            error(message('bioinfo:BioSeq:BioSeq:InvalidFASTAFormat', in.InputFile))
                        end
                        
                        obj.Index = bioinfoprivate.BIFSequenceDataAdapter( ...
                            in,...  % data container
                            {'Sequence','Header'},...
                            {'string','string'},...
                            [false,false],...  % searchable fields
                            [false,false],...    % efficient fields
                            {'Sequence','Header'},...  % field binding
                            [0,0]); % columnPosition
                    end
                elseif strncmpi('.fastq',fileExt,6)
                    if inMemory
                        in = fastqread(in, 'trimheaders', true);
                        obj.Index = bioinfoprivate.structSequenceDataAdapter( ...
                           in,...  % data container
                           {'Sequence','Header' },...
                           {'string'  ,'string' },...
                           [false     ,false    ],...% searchable fields
                           [true      ,true     ],...% efficient fields
                           {'Sequence','Header' });  % field binding
                    else
                        if isempty(indexDir)
                            in = BioIndexedFile('fastq',in,'indexedbykeys',false,'verbose',false);
                        else
                            in = BioIndexedFile('fastq',in,indexDir,'indexedbykeys',false,'verbose',false);
                        end
                        
                        % Validate that what was indexed is a FASTAQ formatted file
                        try
                           fastqread(getEntryByIndex(in,1));
                        catch ME 
                            try
                                delete(in.IndexFile)
                            catch ME2 %#ok<NASGU>
                                %could not delete, just continue with error
                            end                                      
                            error(message('bioinfo:BioSeq:BioSeq:InvalidFASTQFormat', in.InputFile))
                        end
                        
                        obj.Index = bioinfoprivate.BIFSequenceDataAdapter( ...
                            in,...  % data container
                            {'Sequence','Header'},...
                            {'string','string'},...
                            [false,false],...  % searchable fields
                            [false,false],...    % efficient fields
                            {'Sequence','Header'},...  % field binding
                            [0,0]); % columnPosition
                    end
                else
                      error(message('bioinfo:BioSeq:BioSeq:InvalidFile'))
                end                  

            %=== Build BioSeq from MATLAB structure (such as returned by
            %    SAMREAD,BAMREAD,FSTAREAD, or FASTQREAD)
            elseif isstruct(in)
                if ~iscellstr(in(1).Sequence) % structure is an array of structures such as returned by file readers
                    if sum(isfield(in,{'Sequence', 'Header'})) == 2 % structure from fastqread or fastaread
                        fieldBinding = {'Sequence','Header'};
                    elseif sum(isfield(in,{'Sequence', 'QueryName'})) == 2 % structure from bamread or samread
                        fieldBinding = {'Sequence','QueryName'};
                    elseif isfield(in,'Sequence') % structure without header
                        fieldBinding = {'Sequence','Header'};
                    else
                        error(message('bioinfo:BioSeq:BioSeq:InvalidInputStructure'))
                    end
                elseif sum(isfield(in, properties(obj))) == numel(properties(obj)) % structure from get(obj)
                    varargin = {'Sequence',in.Sequence,'Header',in.Header,varargin{1:end}};
                    fieldBinding = {'Sequence','Header'};
                    in = [];
                else
                    error(message('bioinfo:BioSeq:BioSeq:InvalidInputStructure'))
                end
                obj.Index = bioinfoprivate.structSequenceDataAdapter( ...
                    in,...  % data container
                    {'Sequence','Header' },...
                    {'string'  ,'string' },...
                    [false     ,false    ],...% searchable fields
                    [true      ,true     ],...% efficient fields
                    fieldBinding); % field binding
            elseif isempty(in)
                obj.Index = bioinfoprivate.structSequenceDataAdapter( ...
                    in,...  % data container
                    {'Sequence','Header' },...
                    {'string'  ,'string' },...
                    [false     ,false    ],...% searchable fields
                    [true      ,true     ],...% efficient fields
                    {'Sequence','Header'}... % field binding
                    );
            else
                error(message('bioinfo:BioSeq:BioSeq:InvalidInput'))
            end
            
            %=== Populate the object with PV pairs
            if rem(numel(varargin),2) ~= 0
                error(message('bioinfo:BioSeq:BioSeq:InvalidInputNumber'))
            end      
            for k = 1:2:numel(varargin)
                %=== make sure Sequence PVP is handled first, so size of BioSeq is predetermined.
                seqPar = find(~cellfun(@isempty, strfind(varargin(1:2:end), 'Sequence')));
                if ~isempty(seqPar)
                    varargin = {varargin{seqPar(1)*2-1} varargin{seqPar(1)*2} varargin{:}}; %#ok<CCAT> % add Sequence PVP at beginning
                    varargin([seqPar(1)*2+1 seqPar(1)*2+2]) = ''; % remove Sequence PVP from original position
                end
                obj.(varargin{k})= varargin{k+1};
            end
            
        end % end constructor             

        %==================================================================
        % get.NSeqs
        %==================================================================
        function n = get.NSeqs(obj)
            n = obj.Index.NumberOfEntries;
        end
		
		%==================================================================
		% set.Name
		%==================================================================
		function obj = set.Name(obj, n)
			if ~ischar(n) || size(n,1) > 1
                error(message('bioinfo:BioSeq:BioSeq:InvalidName'))
			end
			obj.Name = n;
        end
		
        %==================================================================
		% get.Header
		%==================================================================
        function h = get.Header(obj)
            h = getField(obj.Index,'Header');
        end
		%==================================================================
		% set.Header
		%==================================================================
		function obj = set.Header(obj,h)
            obj.Index = setField(obj.Index,'Header',h);
        end
        %==================================================================
		% get.Sequence
		%==================================================================
        function h = get.Sequence(obj)
            h = getField(obj.Index,'Sequence');
        end
		%==================================================================
		% set.Sequence
		%==================================================================
		function obj = set.Sequence(obj,h)
            obj.Index = setField(obj.Index,'Sequence',h);
        end                

function out = get(obj, pname)
%GET retrieve a BioSeq (or derived) object property.
%
%   GET(OBJ) returns all properties of BioSeq (or derived) object OBJ in a
%   scalar structure, where each field name is a property of OBJ, and each
%   field contains the value of that property.
%
%   GET(OBJ, PROP) returns the value of one or more properties of BioSeq (or
%   derived) object OBJ. PROP is a MATLAB string with the property name or
%   a 1-by-N or N-by-1 cell array of strings containing multiple property
%   names. When the input PROP is a cell, the output is also a cell
%   containing the values for each requested property.
%
%   Examples:
%   % Create a BioSeq object from variables in the workspace.
%   seqs = {randseq(10); randseq(15); randseq(20)};
%   heads = {'H1'; 'H2'; 'H3'};
%   obj = BioSeq(seqs, heads)
%
%   % Retrieve the values of the 'Header' property.
%   get(obj, 'Header')
%
%   % Retrieve the values of the 'Header' and 'Sequence' properties.
%   v = get(obj, {'Header', 'Sequence'})
%
%   % Transform the BioSeq object OBJ into a structure.
%   str = get(obj)
%
%   See also BIOMAP, BIOMAP/GETSUBSET, BIOREAD.

%   Copyright 2012 The MathWorks, Inc.

checkScalarInput(obj);

% If property name/s given then validate them
if nargin > 1
    if ischar(pname) && isvector(pname)
        cellOutput = false;
        pname = {pname};
    elseif iscellstr(pname) && isvector(pname) && numel(unique(pname))==numel(pname)
        cellOutput = true;
    else
        error(message('bioinfo:BioSeq:get:IllegalPropertyName'))
    end
else
    pname = properties(obj);
end

% Find properties that are available in the container:
g = isValidField(obj.Index,pname);

% Initialize out
out = cell(size(pname));

% Get properties from container
if any(g)
    out(g) = obj.Index.getField(pname(g));
end

% Fill in with dot notation the rest of the properties that
% are not in the container:
for i = find(~g(:))'
    out{i} = subsref(obj,substruct('.',pname{i}));
end

if nargin == 1
    % Change output to struct if all props requested i.e. get(obj)
    out = cell2struct(out(:), pname(:), 1);
elseif ~cellOutput
    % Change output to not cell if only one property
    % requested, e.g. get(obj,'Sequence')
    out = out{1};
end
end   

	end % methods
	
	%======================================================================
	% HIDDEN METHODS
	%======================================================================
	
	methods (Hidden = true)

		%==================================================================
		% disp
		%==================================================================
		function disp(obj)
			%DISP display the properties of object OBJ.
			%
			%   DISP(OBJ) displays the properties of obj OBJ and their
			%   values. If OBJ is constructed from an index file, then DISP
			%   indicates so, otherwise DISP calls the built-in MATLAB
			%   disp.
			
			%=== if obj is built on an indexed file, overwrite so that
			%    properties are not accessed in the source file
			if numel(obj)==1 && ~obj.Index.InMemory
				str = evalc(['builtin(''disp'',' class(obj) '(''Name'',obj.Name))']);
                props = properties(obj);
                for i = 1:numel(props)
                    if isValidField(obj.Index,props{i})
                        str = regexprep(str,...
                            sprintf(' %s:([^\n]+)',props{i}),...
                            sprintf(' %s: [%dx1 File indexed property]',props{i},obj.NSeqs));
                    elseif strcmp(props{i},'NSeqs')
                        str = regexprep(str,...
                            sprintf(' %s:([^\n]+)',props{i}),...
                            sprintf(' %s: %d',props{i},obj.NSeqs));
                    elseif strcmp(props{i},'SequenceDictionary')
                        if numel(obj.SequenceDictionary)==1
                            str = regexprep(str,...
                                sprintf(' %s:([^\n]+)',props{i}),...
                                sprintf(' %s: ''%s''',props{i},obj.SequenceDictionary{1}));
                        else
                            str = regexprep(str,...
                                sprintf(' %s:([^\n]+)',props{i}),...
                                sprintf(' %s: {1x%d cell}',props{i},numel(obj.SequenceDictionary)));
                        end
                    end
                end
				if feature('hotlinks')==0
                    str = regexprep(str, '<[^>]+>','');
				end
				disp(str)
			else
				builtin('disp',obj)
			end
        end
        
        %==================================================================
        % subsasgn
        %==================================================================
        function obj = subsasgn(obj,s,in)
            %SUBSASGN subscripted assignment.
            
            %=== Only intercept dot-subsasgn and when the property is
            %    stored in the data container which is file indexed,
            %    anyother subsasgn is passed to the built-in subsasgn.
            %    We prevent MATLAB fetching the whole property only to find
            %    later that the adaptor can not write back to a file.
            if s(1).type(1) == '.' && isValidField(obj.Index,s(1).subs) && ~obj.Index.InMemory
                error(message('bioinfo:BioSeq:BioSeq:SetInvalid',s(1).subs,class(obj)))
            end    
            
            %=== Prevent invalid access to object members:
            if s(1).type(1) == '.' && ~isValidField(obj.Index,s(1).subs)
                mc = metaclass(obj);
                prop = findobj(mc.PropertyList,'Name',s(1).subs);
                % Is it a non-public property?
                if ~isempty(prop) && ~strcmp(prop.SetAccess,'public')
                    error(message('MATLAB:class:SetProhibited',prop.Name,class(obj)))
                end
                % Is it a non-public method?
                meth = findobj(mc.MethodList,'Name',s(1).subs);
                if ~isempty(meth) && ~strcmp(meth.Access,'public')
                    error(message('MATLAB:class:MethodRestricted',meth.Name,class(obj)))
                end
            end
            
            obj = builtin('subsasgn',obj,s,in);
        end
        
        %==================================================================
        % subsref
        %==================================================================
        function varargout = subsref(obj,s)
            %SUBSREF subscripted reference for BioSeq (or derived) object.
            %   B = SUBSREF(A,S) is called for the syntax A.PropertyName when A is a
            %   BioSeq (or derived) object.  S is a structure array with the fields:
            %       type -- string containing '()', '{}', or '.' specifying the
            %               subscript type.
            %       subs -- Cell array or string containing the actual subscripts.
            %
            %   B = A.PROPERTYNAME returns the values of a property.
            %   A.PROPERTYNAME may also be followed by further subscripting as
            %   supported by the property. In particular, A.PROPERTYNAME(X) is a valid
            %   syntax, where X is a numeric array, a logical array or a cell array of
            %   strings corresponding to the values of the 'Header' property.
            %
            %  Examples:
            %  % Create a BioRead object from variables in the workspace.
            %  seqs = {randseq(10); randseq(15); randseq(20)};
            %  heads = {'H1', 'H2', 'H3'};
            %  obj = BioRead(seqs, heads);
            %
            %  % Access properties through subscripted reference.
            %  obj.Header(1:2)
            %  obj.Header([true true false])
            %  obj.Sequence([1 2])
            %  obj.Name
            %
            %  See also BIOMAP, BIOMAP/GET, BIOREAD.
            
            
            %=== Only intercept dot-subsref and when the propertie is
            %    stored in the data container, anyother subsref is passed
            %    to the built-in subsref.
            intercept = false;
            if s(1).type(1) == '.' 
                checkScalarInput(obj);
                if isValidField(obj.Index,s(1).subs)
                    intercept = true;
                end
            end
            
            if intercept
                propName = s(1).subs;
                %=== propName must be a valid property in the object since some
                %    containers may have more indexed data available
                if any(strcmp(propName, properties(obj)))
                    if numel(s) == 1 % get all the property's entries (e.g. obj.Property)
                        varargout{1} = obj.Index.getField(propName);
                    elseif (numel(s)==2) && (s(2).type(1) == '(') && (numel(s(2).subs)==1) %e.g. obj.Property(x)
                        if isempty(s(2).subs{1})
                            varargout{1} = obj.Index.getField(propName,[]);
                        elseif islogical(s(2).subs{1})
                            X = find(s(2).subs{1});
                            varargout{1} = obj.Index.getField(propName,X); %#ok<FNDSB>
                        elseif iscell(s(2).subs{1}) 
                            X = getIndexByHeader(obj,s(2).subs{1});
                            varargout{1} = obj.Index.getField(propName,X);
                        else
                            if ischar(s(2).subs{1})
                                if strcmp(s(2).subs{1},':')
                                    varargout{1} = obj.Index.getField(propName);
                                else
                                    %replicates the behavior of subsindex even if it is a valid header name
                                    varargout{1} = obj.Index.getField(propName,0+s(2).subs{1});
                                end
                            else
                                varargout{1} = obj.Index.getField(propName,s(2).subs{1});
                            end
                        end
                    elseif (numel(s)>1) && (s(2).type(1) == '{') && (numel(s(2).subs)==1) %e.g. obj.Property{x}
                        s(2).type(1) = '(';
                        varargout = subsref(obj,s(1:2));
                        if ~iscell(varargout)
                            error(message('bioinfo:BioSeq:BioSeq:cellRefFromNonCell'))
                        end
                        if numel(varargout)<nargout
                            error(message('bioinfo:BioSeq:BioSeq:needMoreRhsOutputs'))
                        end
                        if (numel(s)>2)  %e.g. obj.Property{x}(y)
                            if numel(varargout)==1
                                varargout{1} = subsref(varargout{1},s(3:end));
                            else
                                error(message('bioinfo:BioSeq:BioSeq:badCellRef'))
                            end
                        end
                    else
                        error(message('bioinfo:BioSeq:BioSeq:IllegalSubscripting', propName))
                    end
                else
                    error(message('bioinfo:BioSeq:BioSeq:noSuchMethodOrField', propName, class( obj )))
                end
            else
                %=== Call builtin subsref for everything else we did not intercept
                [varargout{1:nargout}] = builtin('subsref',obj,s);
            end
        end
        
        %==================================================================
        % getNSeqs (discourage use since 12b, use obj.NSeqs instead)
        %==================================================================
        function out = getNSeqs(obj)
            %GETNSEQS retrieve the number of sequences of a BioSeq (or derived) object.
            %
            %   OUT = GETNSEQS(OBJ) retrieves the value of the 'NSeqs' property of a
            %   BioSeq (or derived) object OBJ, i.e. the number of sequences in the object.
            checkScalarInput(obj);
            out = obj.NSeqs;
        end    
        
        %==================================================================
        % getName (discourage use since 12b, use obj.Name instead)
        %==================================================================
        function out = getName(obj)
            %GETNAME retrieve the name of a BioSeq (or derived) object.
            %
            %   OUT = GETNAME(OBJ) retrieves the 'Name' property of a BioSeq (or
            %   derived) object OBJ.
            checkScalarInput(obj);
            out = obj.Name;
        end
        
        %==================================================================
        % setName (discourage use since 12b, use obj.Name instead)
        %==================================================================
        function obj = setName(obj, n)
            %SETNAME set the name of a BioSeq (or derived) object.
            %
            %   SETNAME(OBJ, Y) sets the 'Name' property of a BioSeq (or derived)
            %   object OBJ to Y. Y must be a string.
            bioinfochecknargin(nargin,2,'BioSeq:setName')
            checkScalarInput(obj);
            obj.Name = n;
        end               
 
        
        % Implementation of VariableEditorPropertyProvider to customize 
        % the display of properties in the Variable Editor
        function isVirtual = isVariableEditorVirtualProp(obj,prop)
            % Return true for the file indexed properties to enable the
            % Variable Editor to derive their display without actually
            % accessing the property.
            isVirtual  = numel(obj)==1 && ~obj.Index.InMemory && ...
                any(strcmp(prop,{'SequenceDictionary','Reference','Signature',...
                  'Start','MappingQuality','Flag','MatePosition','Quality','Sequence','Header'}));        
        end
        function isComplex = isVariableEditorComplexProp(~,~)
            % Diagnostics property should not be complex
            isComplex = false;
        end
        function isSparse = isVariableEditorSparseProp(~,~)
            % Diagnostics property should not be sparse
            isSparse = false;
        end 
        function className = getVariableEditorClassProp(~,~)
            % Diagnostics property in the Variable Editor is dataset object
            className = 'File indexed property';
        end
        function sizeArray = getVariableEditorSize(this,~)
            sizeArray = [this.NSeqs; 1];
        end 
        
    end
    
	%======================================================================
	% PRIVATE METHODS
	%======================================================================
	
	methods (Access = 'protected', Hidden = true)
		
        %==================================================================
		% getProperty
		%==================================================================
        function out = getProperty(obj,pname,X)
            checkScalarInput(obj);
            %=== Input check
            bioinfochecknargin(nargin, 2, ['BioSeq:' mfilename])
            %=== If X is not specified, then all entries are to be retrieved
            if nargin == 2
                out = subsref(obj, substruct('.',pname));
            %=== Empty set
            elseif isempty(X)
                out = subsref(obj, substruct('.',pname, '()', {[]}));
            %=== Numeric index
            elseif (isnumeric(X) && isvector(X) && ~islogical(X))
                out = subsref(obj, substruct('.',pname, '()', {X}));
            %=== Logical index
            elseif (islogical(X) && isvector(X))
                out = subsref(obj, substruct('.',pname, '()', {find(X)}));
            %=== Single string  || Cell array of strings
            elseif iscellstr(X) || (ischar(X) && isvector(X))
                X = getIndexByHeader(obj, X);
                out = subsref(obj, substruct('.',pname,'()',{X}));
            else
                error(message('bioinfo:BioSeq:BioSeq:InvalidIndex'))
            end
        end

        %==================================================================
        % setProperty
        %==================================================================
        function obj = setProperty(obj,pname,isCell,Y,X)
            checkScalarInput(obj);
            %=== If X is not specified, then all elements are to be set
            %using Y we use directly obj.(pname)
            if nargin == 4
                obj.(pname) = Y;
                return
            end
            [X,ex] = uniformizeIndexToNumericArray(obj,X,numel(Y));
            if isCell
                if ~iscellstr(Y)
                    error(message('bioinfo:BioSeq:BioSeq:ValueNotCell',pname));
                end
            else
                if ~isnumeric(Y)
                    error(message('bioinfo:BioSeq:BioSeq:ValueNotNumeric',pname));
                end
            end
            if ~isempty(ex)
                Y = Y(ex); % expand to account for duplicate Headers
            end
            if numel(Y) ~= numel(X)
                error(message('bioinfo:BioSeq:BioSeq:IncorrectSize',pname));
            end
            %=== Set Sequence to the new values Y
            obj.(pname)(X) = Y(:);
        end

        %==================================================================
		% uniformizeIndexToNumericArray
		%==================================================================
        function [X,ex] = uniformizeIndexToNumericArray(obj,X,nY)
            %UNIFORMIZEINDEXTONUMERICARRAY guaranties a numeric and valid index.
            %
            % [X,ex] = UNIFORMIZEINDEXTONUMERICARRAY(OBJ,X,nY)
            % returns numeric index that may be used to access the data in
            % the object properties. When object is accessed by Header,
            % headers may be duplicated; therefore ex has an expanding
            % index so the number of values to set equals the number of
            % matching headers. 
            
            hasDuplicate = false;

            ex = [];
            if islogical(X) && isvector(X)
                X = find(X);
            elseif iscellstr(X) || (ischar(X) && isvector(X))
                [X, hasDuplicate, ex] = getIndexByHeader(obj, X);
            elseif ~isnumeric(X) || ~isvector(X)
                error(message('bioinfo:BioSeq:BioSeq:InvalidIndexType'));
            end
            %=== check for out-of-bounds
            if any(X > obj.getNSeqs)
                error(message('bioinfo:BioSeq:BioSeq:IndexOutOfBounds'));
            end
            
            if (nargin>2)
                if (nY~=numel(X))
                    % When the provided number of values in Y does not match
                    % the number of indices (X), we'll try to expand Y 
                    if hasDuplicate && nY ~= length(unique(ex))
                        error(message('bioinfo:BioSeq:BioSeq:IncorrectHeaderSizeWithDuplicates'));
                    end
                else  % even if hasDuplicate, the provided number of values in Y
                      % matches the number of target locations
                    ex = [];
                end
            end
            
        end
        
		%==================================================================
		% getIndexByHeader
		%==================================================================
		function [idx, hasDuplicate, idxDupl] = getIndexByHeader(obj, X)
			%GETINDEXBYHEADER return indexes corresponding to specified 'Header' values.
			%
			% [IDX,HASDUPLICATE,EXPIDX] = GETINDEXBYHEADER(OBJ,X) returns
			% numeric indices IDX of the 'Header' values specified in input
			% X. X is a cellstr or a string. When OBJ is not constructed
			% from BioIndexedFile, it also returns a flag HASDUPLICATE
			% specifying whether the given headers have multiple entries
			% and expand index EXPIDX that indicates which of the specified
			% headers in X correspond to each of the returned indexes in IDX. 
			
            hasDuplicate = false; % scalar flag, true if there are duplicate headers indexed
            idxDupl = []; % index in X of each of the indexed in IDX (not used with BioIndexedFile)
            idx = [];
            
            headers = getHeader(obj);
            if isempty(headers)
                error(message('bioinfo:BioSeq:BioSeq:UnavailableHeaders'))
            end
            
            
            %=== determine which entry correspond to the given list of headers X
            if ischar(X)
                idx = find(strcmp(X, headers));
                idx = idx(:); % make sure it is a column vector
                if isempty(idx)
                    error(message('bioinfo:BioSeq:BioSeq:HeaderNotFound'))
                end
                if length(idx) > 1
                    hasDuplicate = true;
                    idxDupl = ones(length(idx),1);
                end
            else
                for j = 1:numel(X)
                    f = find(strcmp(X{j}, headers));
                    f = f(:); % make sure it is a column vector
                    if isempty(f)
                        error(message('bioinfo:BioSeq:BioSeq:HeaderNotFound'))
                    end
                    if length(f) > 1
                        hasDuplicate = true;
                        fidx = repmat(j, length(f), 1);
                    else
                        fidx = j;
                    end
                    
                    idx = [idx; f];  %#ok<AGROW>
                    idxDupl = [idxDupl; fidx]; %#ok<AGROW>
                end
            end
            
            
            %=== reset idx to empty if there are no duplicates
            if ~hasDuplicate
                idxDupl = [];
            end
        end

		%==================================================================
		% checkScalarInput
		%==================================================================        
        function checkScalarInput(obj)
            if numel(obj)>1
                error(message('bioinfo:BioSeq:BioSeq:NonScalarObject'))
            end
        end

		%==================================================================
		% getAdapterDetails
		%==================================================================        
        function details = getAdapterDetails(obj)
            details = getAdapterDetails(obj.Index);
        end        
		
	end % methods
	
end % BIOSEQ


