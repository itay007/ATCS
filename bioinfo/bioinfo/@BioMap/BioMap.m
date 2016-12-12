classdef BioMap < BioRead
	%BIOMAP class representing a collection of sequences  and their alignment information.
	%
	%  A BioMap object is a generic representation of mapped sequence data. It
	%  consists of a set of sequences, qualities and information relative to
	%  how each sequence aligns to a given reference.
	%
	%  BioMap properties:
	%  Reference          - name of the reference that aligns to each sequence.
	%  Signature          - representation of how each sequence aligns to the reference.
	%  Start              - start position of each aligned sequence relative to the reference.
	%  MappingQuality     - quality of each alignment.
	%  Flag               - bitwise representation of alignment information.
	%  MatePosition       - position of a mapped mate for each sequence.
	%  Quality            - quality of each sequence.
	%  Sequence           - sequence data.
	%  Header             - header information for each sequence.
    %  SequenceDictionary - reference sequences used in the object.
	%  Nseqs              - number of sequences represented in the object.
	%  Name               - name of the object.
	%
    %  BioMap specialized property access:
	%  get                - get multiple properties with a single statment.
	%  getFlag            - subscripted retrieval of the 'Flag' property.
	%  getHeader          - subscripted retrieval of the 'Header' property.
	%  getMappingQuality  - subscripted retrieval of the 'MappingQuality' property.
	%  getMatePosition    - subscripted retrieval of the 'MatePosition' property.
    %  getQuality         - subscripted retrieval of the 'Quality' property.  
	%  getReference       - subscripted retrieval of the 'Reference' property.
	%  getSequence        - subscripted retrieval of the 'Sequence' property.
    %  getSignature       - subscripted retrieval of the 'Signature' property.
	%  getStart           - subscripted retrieval of the 'Start' property.
	%  getSubsequence     - subscripted retrieval of partial values of the 'Sequence' property.    
	%  set                - set multiple properties with a single statment.
	%  setFlag            - subscripted assignment of the 'Flag' property.
	%  setHeader          - subscripted assignment of the 'Header' property.
	%  setMappingQuality  - subscripted assignment of the 'MappingQuality' property.
	%  setMatePosition    - subscripted assignment of the 'MatePosition' property.
    %  setQuality         - subscripted assignment of the 'Quality' property.     
	%  setReference       - subscripted assignment of the 'Reference' property.
	%  setSequence        - subscripted assignment of the 'Sequence' property.
	%  setSignature       - subscripted assignment of the 'Signature' property.
	%  setStart           - subscripted assignment of the 'Start' property.
	%  setSubsequence     - subscripted assignment of partial values of the 'Sequence' property.
    %
	%  BioMap methods:
	%  BioMap              - create a BioMap object.
	%  combine             - combine two objects.
	%  filterByFlag        - filter the elements according to specific criteria.
	%  getAlignment        - reconstruct the alignment within a given region.
	%  getBaseCoverage     - compute the base-by-base coverage within a given region.    
	%  getCompactAlignment - reconstruct a compact alignment within a given region.
    %  getCounts           - counts the number of mapped reads within a given region.
    %  getIndex            - find short reads mapped within a given region.
	%  getInfo             - retrieve the information relative to a single element.
	%  getStop             - compute the stop position in the alignment.
	%  getSubset           - retrieve a subset of elements.
    %  getSummary          - print summary of a BioMap object.        
    %  plotSummary         - quality summary plots.    
	%  setSubset           - set the values of a subset of elements.
	%
	%  See also BIOMAP, BIOREAD/BIOREAD.
	
	%   Copyright 2009-2012 The MathWorks, Inc.
	
	properties
		%SEQUENCEDICTIONARY names of the reference sequences.
		%     The 'SequenceDictionary' property is a unique cell string
		%     containing the name of all the reference sequences referred
		%     by aligned reads. 
		%
		%  See also BIOMAP, BIOREAD.
		SequenceDictionary = cell(0,1);
    end
	properties (Dependent)
        %REFERENCE  reference name for each mapped sequence.
        %    The 'Reference' propertty is a cell array of strings with the 
        %    reference names for each mapped sequence.
        Reference 
        
		%SIGNATURE representation of how each sequence aligns.
		%     The 'Signature' property is a cell array of strings
		%     containing compact representations of how each sequence
		%     aligns against the reference sequence. These signatures
		%     consist of symbols and numbers: the symbols represent
		%     operations such as 'match', 'insertion', 'deletion',
		%     'padding', 'clipping', etc., whereas the number represents
		%     the positions that are concerned to a specific operation.
		%     CIGAR strings are an example of signature formats.
		%
		%  See also BIOMAP, BIOREAD.
		Signature
	
 		%START position in the reference sequence at which each sequence begins to align.
 		%     The 'Start' property is an array of positive integers holding
 		%     the start position in the reference sequence at which each
 		%     sequence alignment begins.
 		%
 		%  See also BIOMAP, BIOREAD.
 		Start
		
		%MAPPINGQUALITY score representing the mapping score for each sequence.
		%     The 'MappingQuality' property is an array of integers holding
		%     the mapping score of each sequence.
		%
		%  See also BIOMAP, BIOREAD.
		MappingQuality
		
		%FLAG numeric representation of mapping information for each sequence.
		%     The 'Flag' property is an array of integers representing
		%     various details of the mapping of each sequence. Each flag
		%     must be interpreted bitwise (0 = false, 1 = true) to obtain
		%     the following information:
		%     0x001: the read is paired in sequencing
		%     0x002: the read is mapped in a proper pair
		%     0x004: the read is unmapped
		%     0x008: the mate is unmapped
		%     0x010: strand of the read (0 = forward, 1 = reverse)
		%     0x020: strand of the mate (0 = forward, 1 = reverse)
		%     0x040: the read is the first in a pair
		%     0x080: the read is the second in a pair
		%     0x100: the alignment is not primary
		%     0x200: the read fails platform/vendor quality checks
		%     0x400: the read is a PCR or optical duplicate
		%
		%  See also BIOMAP, BIOREAD.
		Flag
		
		%MATEPOSITION position of the mate read.
		%     The 'MatePosition' property is an array of integers
		%     corresponding to the positions of read mates for each read. 
		%
		%  See also BIOMAP, BIOREAD.        
		MatePosition
    end
	
    properties (Access = private)
        DictionaryMapping = [];
    end
	%======================================================================
	% PUBLIC METHODS
	%======================================================================
	
	methods
		
		%==================================================================
		% Constructor
		%==================================================================
		
		function obj = BioMap(in, varargin)
			%BIOMAP Create a BioMap object.
			%
            %  B = BIOMAP(FILENAME) creates a BioMap object B from the content of a 
            %  SAM/BAM formatted file FILENAME. Data is kept in the file and the BioMap
            %  object accesses it using auxiliary index files. When FILENAME is a SAM
            %  formatted file BioMap uses an index file with the extension .IDX, for 
            %  BAM formatted files BioMap uses two auxiliary files with the extensions
            %  .BAI and .LINEARINDEX. When the auxiliary files are not present BioMap 
            %  creates them in the same directory where the source data is.
            %
			%  B = BIOMAP(STRUCT) creates a BioMap object B from the content of a 
            %  MATLAB structure STRUCT, such as the structure returned by the SAMREAD 
            %  and BAMREAD functions.
            %
            %  B = BIOMAP(..., 'SelectReference', R) selects one or more references 
            %  when the source data contains sequences mapped to more than one reference. 
            %  R defaults to all the references in the header dictionary of the source file 
            %  or when the header dictionary is not available R defaults to the all the 
            %  reference names found in the source data. R is either a string or a cell 
            %  string. By using this option you can prevent BioMap constructor from creating 
            %  auxiliary index files for references that you will not use in your analysis.
            %
            %  B = BIOMAP(..., 'InMemory', true) forces data to be loaded into memory 
            %  when using an auxiliary index file. InMemory defaults to false. InMemory 
            %  is ignored when the first input is not a file name.  
            %
            %  B = BIOMAP(..., 'IndexDir', INDEXDIR) sets the path to the directory 
            %  where the auxiliary index files can be found or created. 
            %
            %  B = BIOMAP(..., PROPERTY_NAME , VALUE) specifies optional parameter 
            %  name/value pairs to further modify the object properties when allowed. 
            %  PROPERTY_NAME is a string.
			%
			%  Examples:
			%
			%  % Create a BioMap object from a SAM formatted file.
			%  obj = BioMap('ex1.sam')
			%
			%  % Create a BioMap object from a valid structure and set the name.
			%  str = samread('ex1.sam'); 
			%  obj = BioMap(str, 'Name', 'BioMap object from ex1.sam file')
			%
			%  See also BIOMAP, BIOREAD/BIOREAD.
			
            
            if (nargin == 0)
                in = [];
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
                        inMemory = bioinfoprivate.opttf(varargin{i+1}, 'InMemory', ['BioMap:' mfilename]); 
                    end
                end
            end
                            
            %=== Check if within the PVPs, subsetting is requested:
            doSubset = false;
            for i = 1:2:numel(varargin)-1
                pname = varargin{i};
                if ischar(pname)
                    if strncmpi('selectreference',pname,max(2,numel(pname)))
                        removeVarargin = [removeVarargin i i+1]; %#ok<AGROW>
                        doSubset = true;
                        requestedRef = varargin{i+1};
                        if (ischar(requestedRef) && isrow(requestedRef)) || (iscellstr(requestedRef) && (numel(requestedRef)>0))
                            if ~iscellstr(requestedRef)
                                requestedRef = {requestedRef};
                            end                            
                        else
                            error(message('bioinfo:BioMap:BioMap:InvalidSelectReference'))
                        end
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
                             error(message('bioinfo:BioMap:BioMap:InvalidIndexDir'))
                        end
                        if ~isdir(indexDir)
                            error(message('bioinfo:BioMap:BioMap:UnknownIndexDir', indexDir))
                        end
                    end
                end
            end            
            
            %=== Remove PVPs already processed from varargin, all remaining
            %    PVPs shall be only settable object properties
            if numel(removeVarargin)
                varargin(removeVarargin) = [];
            end
            
           
            %=== Build BioMap from FileName
            if ischar(in) && (isrow(in) || isempty(in))
                if ~exist(in,'file')
                    error(message('bioinfo:BioMap:BioMap:FileNotFound', in))
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
                    % Index or open the index file:
                    if isempty(indexDir)
                        bif = BioIndexedFile('sam',in,'indexedbykeys',false,'verbose',false);
                    else
                        bif = BioIndexedFile('sam',in,indexDir,'indexedbykeys',false,'verbose',false);
                    end
                    
                    % Validate that what was indexed is a SAM formatted file
                    try
                        samread(getEntryByIndex(bif,1));
                    catch ME
                        try
                            delete(bif.IndexFile)
                        catch ME2 %#ok<NASGU>
                            %could not delete, just continue with error
                        end
                        error(message('bioinfo:BioMap:BioMap:InvalidSAMFormat', bif.InputFile))
                    end
                    
                    % get the scanned dictionary ordered as
                    % BioIndexedFile.getFeature('Reference') will return
                    % the reference indices.
                    dic = getDictionary(bif);
                    
                    % get the dictionary from the file header
                    info = saminfo(in);
                    if (isfield(info,'SequenceDictionary') && ...
                            numel(info.SequenceDictionary)>0 && ...
                            isfield(info.SequenceDictionary(1),'SequenceName'))
                        obj.SequenceDictionary = {info.SequenceDictionary.SequenceName}';
                    else
                        % Id there is no dictionary in the file header we
                        % default to use the dictionary that was scanned
                        % while indexing the file.
                        obj.SequenceDictionary = dic;
                    end
                    
                    obj.DictionaryMapping = zeros(numel(dic),1,'uint32');
                    for i=1:numel(dic)
                       f = find(strcmp(obj.SequenceDictionary,dic{i}),1);
                       if ~isempty(f)
                           obj.DictionaryMapping(i) = f;
                       end
                    end
                    
                    if ~doSubset
                       requestedRef = obj.SequenceDictionary;
                    end

                    href = false(1,numel(obj.SequenceDictionary));
                    for i=1:numel(requestedRef)
                        k = strcmp(obj.SequenceDictionary,requestedRef{i});
                        if any(k)
                            href(k) = true;
                        else
                            error(message('bioinfo:BioMap:BioMap:UnknownReference'))
                        end
                    end
                    newIdxToSeqDic = cumsum(href(:)).* href(:);
                    g = obj.DictionaryMapping>0;
                    obj.DictionaryMapping(g) = newIdxToSeqDic(obj.DictionaryMapping(g));
                    obj.SequenceDictionary = obj.SequenceDictionary(href);
                    refsub = find(obj.DictionaryMapping>0);
                    refIdx = getFeature(bif,'Reference');
                    logIdx = false(size(refIdx));
                    for i = 1:numel(refsub)
                        logIdx(refIdx==refsub(i)) = true;
                    end
                    if any(logIdx)
                        bif = getSubset(bif,find(logIdx));
                        obj.Index = bioinfoprivate.BIFSequenceDataAdapter( ...
                            bif,...  % data container
                            {'Sequence','Header'   ,'Quality','Signature'  ,'Start'   ,'MappingQuality','Flag'  ,'MatePosition','Stop'       ,'Reference'},...
                            {'string'  ,'string'   ,'string' ,'string'     ,'uint32'  ,'uint8'         ,'uint16','uint32'      ,'uint32'     ,'uint32'   },...
                            [false     ,false      ,false    ,false        ,false     ,false           ,false   ,false         ,false        ,false      ],...  % searchable fields
                            [false     ,false      ,false    ,false        ,true      ,false           ,false   ,false         ,true         ,false      ],...    % efficient fields
                            {'Sequence','QueryName','Quality','CigarString','Position','MappingQuality','Flag'  ,'MatePosition','EndPosition','Reference'},...  % field binding
                            [10,1,11,6,4,5,2,8,0,0]); % columnPosition
                        
                        
                        refIdx = getFeature(bif,'Reference');
                        for i = 1:numel(obj.DictionaryMapping)
                            st = getStart(obj,refIdx==i);
                            if ~issorted(st(st>0))
                                warning(message('bioinfo:BioMap:BioMap:UnsortedReadsInSAMFile'))
                                break; % warn only once
                            end
                        end
                        if inMemory
                            obj = BioMap(get(obj));
                        end
                    else
                        obj = BioMap('SequenceDictionary',obj.SequenceDictionary);
                    end
                    
                    
                elseif strncmpi('.bam',fileExt,4)
                    
                    try
                        if ~doSubset
                            % Starting in R2012b when a reference is not
                            % specified, BioMap defaults to all the
                            % references in the dictionary, [] will intruct
                            % BAMIndexedFile to use all references.
                            requestedRef = [];
                        end
                        if isempty(indexDir)
                            bam = bioinfoprivate.BAMIndexedFile(in,'Reference',requestedRef,'verbose',false);
                        else
                            bam = bioinfoprivate.BAMIndexedFile(in,indexDir,'Reference',requestedRef,'verbose',false);
                        end
                    catch ME
                        if strcmp(ME.identifier,'bioinfo:BAMIndexedFile:BAMIndexedFile:ReferenceNotInBAMFile')
                            error(message('bioinfo:BioMap:BioMap:UnknownReference'))
                        end
                        rethrow(ME)
                    end
                    
                    obj.SequenceDictionary = bam.Reference;
                    
                    
                    obj.Index = bioinfoprivate.BAMSequenceDataAdapter( ...
                        bam,...  % data container
                        {'Sequence','Header','Quality','Signature','Start'   ,'MappingQuality','Flag'  ,'MatePosition','Stop'       ,'Reference'},...
                        {'string'  ,'string','string' ,'string'   ,'uint32'  ,'uint8'         ,'uint16','uint32'      ,'uint32'     ,'uint32'   },...
                        [false     ,false   ,false    ,false      ,false     ,false           ,false   ,false         ,false        ,false      ],...% searchable fields
                        [false     ,false   ,false    ,false      ,true      ,false           ,false   ,false         ,true         ,false      ],...% efficient fields
                        {'Sequence','Header','Quality','Signature','Position','MappingQuality','Flag'  ,'MatePosition','EndPosition','Reference'}... % field binding
                        );

                    if inMemory
                        obj = BioMap(get(obj));
                    end
                    
                else
                    error(message('bioinfo:BioMap:BioMap:InvalidFile'))
                end
                
            %=== Build BioMap from MATLAB structure (such as returned by SAMREAD)
            elseif isstruct(in)
                myNamesSAM = {'QueryName', 'Flag', 'ReferenceName', ...
                    'Position', 'MappingQuality', 'CigarString', ...
                    'MatePosition', 'Sequence', 'Quality'};
                myNamesBAM = {'QueryName', 'Flag', 'ReferenceIndex', ...
                    'Position', 'MappingQuality', 'CigarString', ...
                    'MatePosition', 'Sequence', 'Quality'};
                if sum(isfield(in, myNamesSAM)) == numel(myNamesSAM) %structure from SAM file
                    if doSubset
                        r = {in.ReferenceName};
                        h = false(numel(in),1);
                        for i = 1:numel(requestedRef)
                            g = strcmp(r,requestedRef{i});
                            if all(~g)
                                error(message('bioinfo:BioMap:BioMap:UnknownReference'))
                            end
                            h(g) = true;
                        end
                        in = in(h);
                    end
                    % Create numeric reference
                    [obj.SequenceDictionary,~,refint] = unique({in.ReferenceName}');
                    refint = num2cell(refint);
                    [in(1:numel(in)).Reference] = deal(refint{:});
                    in = rmfield(in,'ReferenceName');
                    
                    fieldBinding = {'Sequence','QueryName','Quality','CigarString','Position','MappingQuality','Flag' ,'MatePosition','Reference'};
                    
                elseif sum(isfield(in, myNamesBAM)) == numel(myNamesBAM)  %structure from BAM file
                    if doSubset
                        [r,~,ri] = unique(cellfun(@(x) sprintf('%d',x),{in.ReferenceIndex},'Uniform',false));
                        h = false(numel(in),1);
                        for i = 1:numel(requestedRef)
                            g = find(strcmp(r,requestedRef{i}),1);
                            if isempty(g)
                                error(message('bioinfo:BioMap:BioMap:UnknownReference'))
                            end
                            h(ri==g) = true;
                        end
                        in = in(h);
                    end
                    % Create numeric reference
                    referenceName = cellfun(@(x) sprintf('%d',x),{in.ReferenceIndex},'Uniform',false);
                    [obj.SequenceDictionary,~,refint] = unique(referenceName(:));
                    refint = num2cell(refint);
                    [in(1:numel(in)).Reference] = deal(refint{:});
                    in = rmfield(in,'ReferenceIndex');
                    
                    fieldBinding = {'Sequence','QueryName','Quality','CigarString','Position','MappingQuality','Flag' ,'MatePosition','Reference'};
                    
                elseif sum(isfield(in, properties(obj))) == numel(properties(obj)) % structure from get(obj)
                    if doSubset
                        error(message('bioinfo:BioMap:BioMap:InputStructureCannotBeSubsetted'))
                    end
                    varargin = {'Sequence', in.Sequence, 'Quality', in.Quality,...
                        'SequenceDictionary', in.SequenceDictionary,...
                        'Header', in.Header, 'Reference', in.Reference, ...
                        'Start', in.Start, 'Signature', in.Signature, ...
                        'MappingQuality', in.MappingQuality, 'Flag', in.Flag, ...
                        'MatePosition', in.MatePosition, varargin{1:end}};
                    in = [];
                    fieldBinding = {'Sequence','QueryName','Quality','Signature','Start','MappingQuality','Flag' ,'MatePosition','Reference'};
                else
                    error(message('bioinfo:BioMap:BioMap:InvalidInputStructure'))
                end
                obj.Index = bioinfoprivate.structSequenceDataAdapter( ...
                    in,...  % data container
                    {'Sequence','Header'   ,'Quality','Signature' ,'Start'   ,'MappingQuality','Flag'  ,'MatePosition' ,'Reference'},...
                    {'string'  ,'string'   ,'string' ,'string'    ,'uint32'  ,'uint8'         ,'uint16','uint32'       ,'uint32'   },...
                    [false     ,false      ,false    ,false       ,false     ,false           ,false   ,false          ,false      ],...% searchable fields
                    [true      ,true       ,true     ,true        ,true      ,true            ,true    ,true           ,true       ],...% efficient fields
                    fieldBinding ... % field binding
                    );
               
            elseif isempty(in)
                obj.Index = bioinfoprivate.structSequenceDataAdapter( ...
                    in,...  % data container
                    {'Sequence','Header'   ,'Quality','Signature' ,'Start'   ,'MappingQuality','Flag'  ,'MatePosition' ,'Reference'},...
                    {'string'  ,'string'   ,'string' ,'string'    ,'uint32'  ,'uint8'         ,'uint16','uint32'       ,'uint32'   },...
                    [false     ,false      ,false    ,false       ,false     ,false           ,false   ,false          ,false      ],...% searchable fields
                    [true      ,true       ,true     ,true        ,true      ,true            ,true    ,true           ,true       ]...% efficient fields
                    );
             else
                error(message('bioinfo:BioMap:BioMap:InvalidInput'))
            end
            
            %=== Populate the object with PV pairs
            if rem(numel(varargin),2) ~= 0
                error(message('bioinfo:BioMap:BioMap:InvalidInputNumber'))
            end
            
            if numel(varargin)
                %=== make sure Sequence and SequenceDictionary PVP are
                %    handled first, so size of BioMap and valid references
                %    are predetermined.
                ho = [find(strcmp(varargin(1:2:end),'SequenceDictionary')) find(strcmp(varargin(1:2:end),'Sequence'))];
                ho = [ho*2-1 setdiff(1:2:numel(varargin),ho*2-1)];
                for k = 1:numel(ho)
                    obj.(varargin{ho(k)})= varargin{ho(k)+1};
                end
            end
            
        end % end constructor

        %==================================================================
        % get.Reference
        %==================================================================   
        function r = get.Reference(obj)
            if isempty(obj.DictionaryMapping)
                r = obj.SequenceDictionary(getField(obj.Index,'Reference'));
            else
                r = obj.SequenceDictionary(obj.DictionaryMapping(getField(obj.Index,'Reference')));
            end
        end       
        %==================================================================
        % set.Reference
        %==================================================================
        function obj = set.Reference(obj, r)
            if (ischar(r)&&isvector(r)&&(numel(obj.SequenceDictionary)<2))
                % Special case: we allow scalar expansion (undocumented) for
                %               backwards compatibility grievance when there is
                %               only one reference in SequenceDictionary
                r = {r};
            else
                if ~iscellstr(r) 
                    error(message('bioinfo:BioMap:BioMap:ReferenceNotCellStr', class( obj )))
                end
                if numel(r) ~= obj.NSeqs              
                    error(message('bioinfo:BioMap:BioMap:InvalidReferenceSize', obj.NSeqs, class(obj), obj.NSeqs))
                end
            end
            [ru,~,rh] = unique(r);
            % If SequenceDictionary has not been assigned, then we build
            % one by default from the input (special case for backwards
            % compatibility grievance)
            if isempty(obj.SequenceDictionary)
                obj.SequenceDictionary = ru(:);
            end
            [i,h] = ismember(ru,obj.SequenceDictionary);
            if ~all(i)
                error(message('bioinfo:BioMap:BioMap:ReferenceNotInDictionary'))
            end
            obj.Index = setField(obj.Index,'Reference',uint32(h(rh))); 
        end
        %==================================================================
		% get.Signature
		%==================================================================
        function h = get.Signature(obj)
            h = getField(obj.Index,'Signature');
        end
		%==================================================================
		% set.Signature
		%==================================================================
		function obj = set.Signature(obj,h)
            obj.Index = setField(obj.Index,'Signature',h);
        end
        %==================================================================
		% get.Start
		%==================================================================
        function h = get.Start(obj)
            h = getField(obj.Index,'Start');
        end
		%==================================================================
		% set.Start
		%==================================================================
		function obj = set.Start(obj,h)
            obj.Index = setField(obj.Index,'Start',h);
        end
        %==================================================================
		% get.MappingQuality
		%==================================================================
        function h = get.MappingQuality(obj)
            h = getField(obj.Index,'MappingQuality');
        end
		%==================================================================
		% set.MappingQuality
		%==================================================================
		function obj = set.MappingQuality(obj,h)
            obj.Index = setField(obj.Index,'MappingQuality',h);
        end
        %==================================================================
		% get.Flag
		%==================================================================
        function h = get.Flag(obj)
            h = getField(obj.Index,'Flag');
        end
		%==================================================================
		% set.Flag
		%==================================================================
		function obj = set.Flag(obj,h)
            obj.Index = setField(obj.Index,'Flag',h);
        end
        %==================================================================
		% get.MatePosition
		%==================================================================
        function h = get.MatePosition(obj)
            h = getField(obj.Index,'MatePosition');
        end
		%==================================================================
		% set.MatePosition
		%==================================================================
		function obj = set.MatePosition(obj,h)
            obj.Index = setField(obj.Index,'MatePosition',h);
        end
		
	end % public methods


	%======================================================================
	% HIDDEN METHODS
	%======================================================================
	
	methods (Hidden = true)

        %==================================================================
        % subsref
        %==================================================================
        function varargout = subsref(obj,s)
            %SUBSREF subscripted reference for BioMap object.
            %   B = SUBSREF(A,S) is called for the syntax A.PropertyName when A is a
            %   BioMap (or derived) object.  S is a structure array with the fields:
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
            %  LIMITATIONS: Subscripting expressions such as A.PropertyName{1:2} do not
            %  return multiple outputs. Subscripting expressions such as
            %  A.PropertyName{X} where X is a cell array of strings are not valid (use
            %  the specific GET methods to retrieve the values of a given property
            %  through a list of 'Header' values, e.g. GETSEQUENCE(A, X)).
            %
            %  Examples:
            %   % Create a BioMap object from a SAM formatted file:
            %   obj = BioMap('ex1.sam')
            %
            %  % Access properties through subscripted reference.
            %  obj.Header(1:2)
            %  obj.Header([true true false])
            %  obj.Sequence([1 2])
            %  obj.NSeqs
            %
            %  See also BIOMAP, BIOMAP/GET, BIOREAD.
            
            %=== First call BioSeq's subsref 
            %varargout = cell(1, max(nargout,1));
            [varargout{1:nargout}] = subsref@BioSeq(obj,s);
            %=== Special case for accessing the Reference using '.'
            if (s(1).type(1) == '.') && strcmp(s(1).subs,'Reference') && isnumeric(varargout{1})
                if isempty(obj.DictionaryMapping)
                    varargout{1} = obj.SequenceDictionary(varargout{1});
                else
                    varargout{1} = obj.SequenceDictionary(obj.DictionaryMapping(varargout{1}));
                end
            end           
        end

        %==================================================================
        % getRange
        %==================================================================
        
        function [minstart,maxend] = getRange(obj)
            checkScalarInput(obj);
            
            if numel(obj.SequenceDictionary)>1
                minstart = zeros(numel(obj.SequenceDictionary),1,'uint32');
                maxend = zeros(numel(obj.SequenceDictionary),1,'uint32');
                for i = 1:numel(obj.SequenceDictionary)
                    [minstart(i),maxend(i)] = getRange(getSubset(obj,'SelectReference',i));
                end
                return
            end
                
            % Gets the range of the mapped sequences. 
            bs = 10000;
            n = obj.NSeqs;
            if n==0
                % There are no reads, just return a safe range
                minstart = 1;
                maxend = 1;
            elseif issorted(getStart(obj,1:min(bs,n)))
                minstart = getStart(obj,1);
                m = max(n-bs+1,1);
                maxend = 0;
                while maxend==0 && m>1
                    maxend = max(getStop(obj,m:n));
                    n = m;
                    m = max(n-bs+1,1);
                end
                maxend = max(maxend,max(getStop(obj,m:n)));
            else %Slow for unsorted BioMap
                st = getStart(obj);
                minstart = min(st(st>0));
                st = getStop(obj);
                maxend= max(st(st>0));
            end    
        end
        
        %==================================================================
        % getRowInCompactAlignment
        %==================================================================
        
		function [row, idx, start] = getRowInCompactAlignment(obj, x1, x2, full)
			% Determine the sequence IDX of reads aligning within reference 
			% position x1 and x2. Also returns in which row each sequence 
			% should be placed to achieve a compact display.
						
			%=== Input check
			bioinfochecknargin(nargin, 4, ['BioMap:' mfilename])
            checkScalarInput(obj);
			
			%=== Error check
			if ~isnumeric(x1) || ~isscalar(x1) || ~isnumeric(x2) || ~isscalar(x2) || ...
					x1 > x2 || x1 < 1 || x2 < 1 || (isa(x1,'float')&&(rem(x1, 1) ~=0)) ...
                    || (isa(x2,'float')&&(rem(x2,1) ~= 0))
                error(message('bioinfo:BioMap:BioMap:InvalidRange'));
			end
			
			%=== determine which entries are mapped within the given range
            if full
                idx = getIndex(obj,x1(:),x2(:),'overlap','full');
            else
                idx = getIndex(obj,x1(:),x2(:));
            end
			if isempty(idx)
				row = [];
				start = [];
				return
			end
			
			%=== determine starts and stops of selected reads
			start = getStart(obj, idx);
			stop = getStop(obj, idx);
			
			%=== filter out entries for which sequence has not been reconstructed (cigar = '*')
			k = start ~= stop;
			start = start(k);
			stop = stop(k);
			idx = idx(k);
			
			%=== sort according to starts in ascending order
			[~, s] = sort(start);
			start = start(s);
			stop = stop(s);
			idx = idx(s);
			
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
            
		end
	end % hidden methods
		
end % BioMap

