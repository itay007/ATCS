classdef BioIndexedFile
	%BIOINDEXEDFILE class allows random read access to text files using an index file.
    %
	%  A BIOINDEXEDFILE object uses an auxiliary index file to store
	%  information that allows the efficient random read access to a text
	%  file with unequal sized entries. The auxiliary index file contains
	%  file offsets to each entry in the source file.
    %
    %  A BIOINDEXEDFILE object lets you access individual entries or a
    %  subset of entries when the source file is too big to fit into
    %  memory. You can access entries using indices or alphanumeric keys.
    %  You can read and parse one or more entries using provided
    %  interpreters or a custom interpreter function.   
    %
    %  BioIndexedFile properties:
    %    FileFormat        - File format of the source file.
    %    InputFile         - Path and filename of the source file.
    %    IndexFile         - Path and filename of the auxiliary index file.
    %    NumEntries        - Number of entries indexed by the object.
    %    IndexedByKeys     - Whether or not access by alphanumeric keys is
    %                        available.
    %    MemoryMappedIndex - Whether or not the indices are memory mapped.
    %    Interpreter       - Function handle to parse entries from the
    %                        source file.  
    %
    %  Interpreter and MemoryMappedIndex are the only properties that can be
    %  changed once the object has been constructed. All other properties
    %  are read-only.
    %        
    %  BioIndexedFile methods:
    %    BioIndexedFile    - Create a BioIndexedFile object.
    %    getEntryByIndex   - Retrieve entries using numeric index.
    %    getEntryByKey     - Retrieve entries using alphanumeric keys.
    %    getIndexByKey     - Retrieve indices using alphanumeric keys.
    %    getKeys           - Retrieve all the alphanumeric keys.
    %    getSubset         - Returns an object that indexes a subset of the 
    %                        entries.
    %    read              - Read entries using a specified interpreter.
    %
    %  Examples:
    %
    %     % Build an index to access a cross reference table between gene names
    %     % and GO terms:
    %     sourcefile = which('yeastgenes.sgd');
    %     gene2go = BioIndexedFile('mrtab',sourcefile,'.','KeyColumn',3,...
    %                              'HeaderPrefix','!')
    %     % Read the entries in the source file associated with the gene YAT2:
    %     getEntryByKey(gene2go,'YAT2')
    %     % Adjust the object interpreter to parse only the GO terms:
    %     gene2go.Interpreter = @(x) regexp(x,'GO:\d+','match');
    %     read(gene2go,'YAT2')
	%
	%  See also TEXTSCAN, BIOMAP/BIOMAP, BIOREAD/BIOREAD.   
    
    %   Copyright 2009-2012 The MathWorks, Inc.
    
    properties (SetAccess = private, GetAccess = public)
        %FILEFORMAT file format of the source file.
		%     The 'FileFormat' property indicates the format of the source
		%     file. It can be any of the following general purpose formats:
        %
        %    'table' - A tab delimited table with multiple columns. Every
        %              text line is a different entry. Alphanumeric keys
        %              are found in a specified column during construction.
        %    'mrtab' - A tab delimited table with multiple columns.
        %              Alphanumeric keys are found in a specified column
        %              during construction. Contiguous text lines with
        %              equal keys are considered a single entry. 
        %    'flat'  - Flat file with concatenated entries separated by a
        %              specified delimiter during construction. 
        %
        %    The 'FileFormat' property can also have any of the following
        %    application specific formats: 'sam', 'fasta', or 'fastq'.
        %
        %  See also BioIndexedFile, BioIndexedFile.BioIndexedFile.
        FileFormat
        %INPUTFILE path and filename of the source file.
        %    The 'InputFile' property contains a char array with the path
        %    and filename of the source file.
        %
        %  See also BioIndexedFile, BioIndexedFile.BioIndexedFile.
        InputFile
        %INDEXFILE path and filename of the auxiliary index file.
        %    The 'IndexFile' property contains a char array with the path
        %    and filename of the auxiliary index file.
        %
        %  See also BioIndexedFile, BioIndexedFile.BioIndexedFile.        
        IndexFile
        %NUMENTRIES number of entries indexed by the object.
        %    The 'NumEntries' property contains a integer indicating the
        %    number of entries that are indexed by the object.
        %
        %  See also BioIndexedFile, BioIndexedFile.BioIndexedFile.        
        NumEntries
        %INDEXEDBYKEYS whether or not access by alphanumeric keys is available.
        %    The 'IndexedByKeys' property is a true/false flag that
        %    indicates if alphanumeric keys exist for each entry in the
        %    source file and have been indexed. If so, the source file can
        %    be accessed using the getEntryByKey method.
        %
        %  See also BioIndexedFile, BioIndexedFile.BioIndexedFile,
        %  BioIndexedFile.getEntryByKey.         
        IndexedByKeys
    end
    properties (Dependent)
        %MEMORYMAPPEDINDEX whether or not the indices are memory mapped.
        %    The 'MemoryMappedIndex' property is a true/false flag that
        %    indicates if the indices are stored in the index file or are
        %    loaded in memory.
        %    
        %    This property can be set to true to unload the indices from
        %    memory or false to load the indices into memory.
        %
        %  See also BioIndexedFile, BioIndexedFile.BioIndexedFile.        
        MemoryMappedIndex
    end     
    properties (SetAccess = public, GetAccess = public)
        %INTERPRETER function handle to parse entries from the source file.
        %    The 'Interpreter' property contains a function handle with the
        %    function that can be used to parse individual or groups of
        %    entries from the source file. The interpreter is used by the
        %    read method.
        %
        %    This property can be set to a new function handle or an empty
        %    array once the object has been constructed.
        %
        %  See also BioIndexedFile, BioIndexedFile.BioIndexedFile,
        %  BioIndexedFile.read.         
        Interpreter
    end
    properties (SetAccess = private, GetAccess = private)
        IndexerType
        OffsetToIndex
        OffsetToKeys
        Keys
        KeyOrder
        Index
        IndexType
        ContiguousEntries
        ContiguousEntriesSourceFile
        EntryDelimiter
        PrivateMemoryMappedIndex
        NumEntriesSourceFile
        PrivateFeatureInMemory
    end
    properties (SetAccess = private, GetAccess = public, Hidden = true)
        SubsetIndex
        FeatureName
        FeatureType
        FeatureOffset
        FeatureData
    end
    properties (Dependent, SetAccess = public, GetAccess = public, Hidden = true)
        FeatureInMemory
    end
    
    methods
		
		%==================================================================
		% Constructor
		%==================================================================
        
        function obj = BioIndexedFile(fileFormat,inputFile,varargin)
        %BIOINDEXEDFILE Create a bioinformatics indexed file object.
        %
        %  BIOINDEXEDFILE(FORMAT,SOURCEFILE) constructs a BioIndexedFile object
        %  that indexes the contents of SOURCEFILE following the parsing rules
        %  defined by FORMAT. SOURCEFILE is a string containing a file name and may
        %  be prefixed by a relative or absolute path. FORMAT is a string
        %  containing one of the following recognized formats:
        %
        %  General purpose formats:
        %
        %    'TABLE' - A tab delimited table with multiple columns. Every
        %              text line is a different entry and keys are found in
        %              the first column.
        %    'MRTAB' - A tab delimited table with multiple columns. Keys are
        %              found in the first column. Contiguous text lines with
        %              equal keys are considered a single entry.
        %    'FLAT'  - Flat file with concatenated entries separated by the
        %              string '//'.
        %
        %  Application specific formats:
        %
        %    'SAM'   - SAM formatted files.
        %    'FASTA' - FASTA formatted files.
        %    'FASTQ' - FASTQ formatted files.
        %
        %  The BIOINDEXEDFILE object uses an auxiliary index file to store
        %  information that allows the efficient direct access to SOURCEFILE. The
        %  name of the auxiliary index file by default is derived from SOURCEFILE by
        %  appending the extension '.idx'. The path to the auxiliary index file
        %  defaults to the same directory as SOURCEFILE. When the auxiliary index
        %  file does not exist, the BIOINDEXEDFILE constructor creates one by
        %  analyzing SOURCEFILE, otherwise, the object is built from the information
        %  stored in the auxiliary index file.
        %
        %  BIOINDEXEDFILE(FORMAT,SOURCEFILE,INDEXDIR) sets a different directory to
        %  write or search for the auxiliary index file.
        %
        %  BIOINDEXEDFILE(FORMAT,SOURCEFILE,INDEXFILE) sets the name of the 
        %  auxiliary index file. INDEXFILE may be prefixed by a relative or
        %  absolute path.
        %
        %  BIOINDEXEDFILE(...,'INDEXEDBYKEYS',TF) indicates whether the object can 
        %  be accessed by keys or not. Defaults to true.
        %
        %  BIOINDEXEDFILE(...,'MEMORYMAPPEDINDEX',TF) indicates whether the object 
        %  leaves the indices in the auxiliary index file and accesses them with 
        %  memory maps (true) or load the indices into memory at construction time 
        % (false). Defaults to true.
        %
        %  BIOINDEXEDFILE(...,'INTERPRETER',FH) sets a function handle that is to be
        %  used when accessing the source file with the READ method. The interpreter
        %  must take one or more concatenated entries of the source file into a
        %  single string and return a single structure or an array of structures
        %  with the interpreted data. When the source file has a general purpose
        %  format (FORMAT is 'TABLE','MRTAB', pr 'FLAT') FH defaults to [] which is
        %  equivalent to FH = @(x) x . For all other application specific formats,
        %  FH defaults to an appropriate function handle to the respective file
        %  reader and usually does not need to be changed.
        %
        %  BIOINDEXEDFILE(...,'VERBOSE',TF) displays information about the object
        %  construction. Defaults to true.
        %
        %  The remaining input arguments are considered only when an auxiliary index
        %  file is to be created and when the source file has a general purpose
        %  format (FORMAT is 'TABLE', 'MRTAB', or 'FLAT'). For all other supported
        %  application specific formats these options have already been predefined
        %  and cannot be changed:
        %
        %  BIOINDEXEDFILE(...,'KEYCOLUMN',KC) when format is 'TABLE' or 'MRTAB', KC
        %  contains an integer indicating the column of the table to be used as the
        %  keys. Defaults to the first column.
        %
        %  BIOINDEXEDFILE(...,'KEYTOKEN',KT) when format is 'FLAT', KT contains a
        %  string that must occur in the entry before key. Defaults to '', which
        %  implies that the key is the first string delimited by blank spaces in
        %  the entry.
        %
        %  BIOINDEXEDFILE(...,'COMMENTPREFIX',CP) specifies a prefix string that
        %  denotes lines that should be ignored by the indexer. Defaults to '',
        %  which implies that the indexer does not check for commented out lines.
        %
        %  BIOINDEXEDFILE(...,'HEADERPREFIX',HP) specifies a prefix string that
        %  denotes the source file header lines and which are ignored by the
        %  indexer.  Defaults to '', which implies that the indexer does not check
        %  for header lines.
        %
        %  BIOINDEXEDFILE(...,'TABLEDELIMITER',TD) when format is 'TABLE' or
        %  'MRTAB', TD sets the delimiter between features in each table line. TD
        %  can be either '\t' (a horizontal tab), ' ' (a blank space), or ','
        %  (a comma). Multiple contiguous blank spaces in the source file are
        %  considered to be single delimiter by the indexer. TD defaults to '\t'.
        %
        %  BIOINDEXEDFILE(...,'ENTRYDELIMITER',ED) when format is 'FLAT', ED sets a
        %  string which delimits multiple entries in a source file. ED defaults to
        %  '//'.
        %
        %  BIOINDEXEDFILE(...,'CONTIGUOUSENTRIES',TF) indicates whether the entries
        %  in the source file are expected to be contiguous or not. Defaults to
        %  false.
        %
        %  Examples:
        %
        %     % Build an index to access a cross reference table between gene names
        %     % and GO terms:
        %     sourcefile = which('yeastgenes.sgd');
        %     gene2go = BioIndexedFile('mrtab',sourcefile,'.','KeyColumn',3,...
        %                              'HeaderPrefix','!')
        %     % Read the entries in the source file associated with the gene YAT2:
        %     getEntryByKey(gene2go,'YAT2')
        %     % Adjust the object interpreter to parse only the GO terms:
        %     gene2go.Interpreter = @(x) regexp(x,'GO:\d+','match');
        %     read(gene2go,'YAT2')
        %
        %  See also BioIndexedFile, BioIndexedFile.getEntryByIndex,
        %           BioIndexedFile.getEntryByKey, BioIndexedFile.getIndexByKey,  
        %           BioIndexedFile.getKeys, BioIndexedFile.read            
        
            fileFormats = {'sam','fastq','fasta','mrtab','table','flat','fixnl'};
            %==== Validate file format
            k = find(strncmpi(fileFormat,fileFormats,numel(fileFormat)));
            if numel(k)==1
                fileFormat = fileFormats{k};
            elseif isempty(k)
                error(message('bioinfo:BioIndexedFile:BioIndexedFile:InvalidFileFormat', fileFormat))
            else
                error(message('bioinfo:BioIndexedFile:BioIndexedFile:AmbiguousFileFormat', fileFormat))
            end
            
            %==== Validate input file name
            if isempty(fileparts(inputFile))
                inputFile = [pwd filesep inputFile];
            end              
            
            if ~exist(inputFile,'file')
                [p,n,e] = fileparts(inputFile);
                error(message('bioinfo:BioIndexedFile:BioIndexedFile:SourceFileNotFound', [ n, e ], p))
            end
            
            %==== Resolve ~ (note that for platforms not accepting a home
            %     directory an error should have already generated)
            if ~isempty(regexp(fileparts(inputFile),'^~','once'))            
                 [~,at] = fileattrib(inputFile);
                 inputFile = at.Name;
            end
            
            %==== Finds out the signature by looking at nvarargin
            if rem(numel(varargin),2)==0
                indexFile = [inputFile,'.idx'];
            elseif ischar(varargin{1})
                indexFile = varargin{1};
                varargin = varargin(2:end);
            end
            
            %==== Parse PV pairs
            switch fileFormat
                case {'sam','fasta','fastq'}
                    [mapKeys,memMapIndex,interpreter,columnWithKeys,keyToken,...
                        contiguousEntries,commentPrefix,headerPrefix,tableDelimiter,...
                        entryDelimiter,verbose] = ...
                        parse_inputs_app_specific_formats(fileFormat,varargin{:});
                case {'table','mrtab','flat'}
                    [mapKeys,memMapIndex,interpreter,columnWithKeys,keyToken,...
                        contiguousEntries,commentPrefix,headerPrefix,tableDelimiter,...
                        entryDelimiter,verbose] = ...
                        parse_inputs_general_purpose_formats(varargin{:});
                otherwise
                    error(message('bioinfo:BioIndexedFile:BioIndexedFile:UnsupportedFileFormat'))
            end
            
            %==== Validate index file name
            if indexFile(end)=='.'
                % the user meant a relative path ending with '.'
                indexFile = [indexFile filesep];
            end
            [p,n,e] = fileparts(indexFile);
            if  isempty(e) && isdir([p,filesep,n])
                % fileparts extracted the last part of a directory as the
                % filename
                p = [p,filesep,n];
                n = '';
            elseif  isempty(e) && isempty(p) && isdir(n)    
                % fileparts extracted the string as a filename, but it is a
                % dir relative to pwd
                p = n;
                n = '';
            elseif isempty(p)
                % the user only passed a file name
                p = fileparts(inputFile); %default dir is taken from source
            end
            if ~isdir(p)
                error(message('bioinfo:BioIndexedFile:BioIndexedFile:InvalidDirectory', p))
            end
            if isempty(n)
                % the user only passed a directory
                [~,n,e] = fileparts(inputFile); %default name derived from source
                indexFile = [p filesep n e '.idx'];
            else
                indexFile = [p filesep n e];
            end
            
            %==== Get the info of the source file
            sourceFileInfo = dir(inputFile);
            
            %==== Define other format specific required variables (these are
            %     not controlled with input arguments)
            switch fileFormat
                case 'sam'
                    indexerType = 'table';
                    featureName     = {'Flag' 'Position' 'MappingQuality' 'MatePosition' 'MateIndex' 'QueryLength' 'EndPosition' 'Reference'};
                    featureType     = {'uint16' 'uint32' 'uint8'          'uint32'       'uint32'    'uint16'      'uint32'      'uint8'};
                    featureOffset   = zeros(1,numel(featureName),'uint64');
                    featureInMemory = {false  false      false            false          false       false         false         false};
                    featureData     = {[]     []         []               []             []          []            []            []};        
                case 'fastq'
                    indexerType = 'fastq';
                    featureName     = {'Length'};
                    featureType     = {'uint16'};
                    featureOffset   = zeros(1,numel(featureName),'uint64');
                    featureInMemory = {false};
                    featureData     = {[]};
                case 'fasta'
                    indexerType     = 'fasta';
                    featureName     = {'Length'};
                    featureType     = {'uint32'};
                    featureOffset   = zeros(1,numel(featureName),'uint64');
                    featureInMemory = {false};
                    featureData     = {[]};
                case 'mrtab'
                    indexerType     = 'mrtab';
                    featureName     = {};
                    featureType     = {};
                    featureOffset   = zeros(1,numel(featureName),'uint64');
                    featureInMemory = {};
                    featureData     = {};
                case 'table'
                    indexerType   = 'table';
                    featureName   = {};
                    featureType   = {};
                    featureOffset = zeros(1,numel(featureName),'uint64');
                    featureInMemory = {};
                    featureData     = {};
                case 'flat'
                    indexerType   = 'flat';
                    featureName   = {};
                    featureType   = {};
                    featureOffset = zeros(1,numel(featureName),'uint64');
                    featureInMemory = {};
                    featureData     = {};
                    contiguousEntries = false;
                otherwise
                    error(message('bioinfo:BioIndexedFile:BioIndexedFile:UnsupportedFileFormat'))
            end
            
            %==== Consolidating some variables that need to be passed to
            %     bioindexermex to avoid complex input checking
            columnWithKeys    = uint8(max(columnWithKeys-1,0));
            tableDelimiter    = sprintf(tableDelimiter);
            entryDelimiter    = sprintf(entryDelimiter);
            mapKeys           = double(logical(mapKeys));
            contiguousEntries = double(logical(contiguousEntries));
            indexType         = uint8(sourceFileInfo.bytes>(2^32-1));
            
            if verbose
                fprintf('Source File: %s\n',sourceFileInfo.name)
                fprintf('   Path: %s\n',fileparts(inputFile))
                fprintf('   Size: %d bytes\n',sourceFileInfo.bytes)
                fprintf('   Date: %s\n',datestr(sourceFileInfo.datenum))
            end
            
            if ~exist(indexFile,'file')
                % Index not found create a new one:
                if verbose
                    fprintf('Creating new index file ...\n')
                end
                
                % Update status in NGS Browser
                if usejava('jvm')                    
                    com.mathworks.toolbox.bioinfo.genome.GenomeDesktop.updateStatusMessage('Creating a new index file...');
                    clearBrowserStatusObj = onCleanup(@()com.mathworks.toolbox.bioinfo.genome.GenomeDesktop.updateStatusMessage('')); 
                end
                
                try
                    [nr,nl] =  BioIndexedFile.bioindexermex(...
                        indexerType,...
                        inputFile,...
                        indexFile,...
                        indexType,...
                        contiguousEntries,...
                        mapKeys,...
                        columnWithKeys,...
                        keyToken,...
                        commentPrefix,...
                        headerPrefix,...
                        tableDelimiter,...
                        entryDelimiter);
                catch ME
                    % When there is an error in the indexer it is likely that
                    % the file left is invalid, then it is deleted.
                    if exist(indexFile,'file')
                        delete(indexFile)
                    end                    
                    
                    % Clear status message
                    clear('clearBrowserStatusObj');
                    
                    bioinfoprivate.bioclsrethrow(mfilename,mfilename,ME)
                end
                if nr==0
                    error(message('bioinfo:BioIndexedFile:BioIndexedFile:indexerError'))
                end
                if verbose
                    fprintf('Indexer found %d entries after parsing %d text lines.\n',nr,nl)
                end
                % Clear status message
                clear('clearBrowserStatusObj');
                
            elseif verbose
                fprintf('Index file exists.\n')
            end
            
            if verbose
                sourceFileInfo = dir(indexFile);
                fprintf('Index File: %s\n',sourceFileInfo.name)
                fprintf('   Path: %s\n',fileparts(indexFile))
                fprintf('   Size: %d bytes\n',sourceFileInfo.bytes)
                fprintf('   Date: %s\n',datestr(sourceFileInfo.datenum))
            end
            
            %==== Object Construction
            
            % Properties that are initialized from the input arguments and
            % other predefined format specific variables at construction time:
            
            obj.IndexerType       = indexerType;
            obj.FileFormat        = fileFormat;
            obj.InputFile         = inputFile;
            obj.IndexFile         = indexFile;
            obj.IndexType         = indexType;
            obj.IndexedByKeys     = mapKeys;
            obj.FeatureName       = featureName;
            obj.FeatureType       = featureType;
            obj.FeatureOffset     = featureOffset;
            obj.FeatureData       = featureData;
            obj.Interpreter       = interpreter;
            
            % Properties that are initialized from the information in the
            % header of a valid index file and updated with readHeaderInfo
            %     obj.NumEntriesSourceFile
            %     obj.OffsetToKeys
            %     obj.OffsetToIndex
            %     obj.ContiguousEntriesSourceFile
            %     obj.FeatureOffset     (already initialized, updated here)
            %     obj.IndexedByKeys     (already initialized, updated here)
            %     obj.EntryDelimiter
            %     obj.IndexType            
            
            obj = readHeaderInfo(obj,true);  % read header of the IDX file
            
            % Properties still needed to be loaded or initialized:
            %     obj.Keys
            %     obj.Index
            %     obj.MemoryMappedIndex
            %     obj.SubsetIndex
            %     obj.NumSubsetEntries

            obj.SubsetIndex       = [];            
            obj.NumEntries  = obj.NumEntriesSourceFile;
            obj.ContiguousEntries = obj.ContiguousEntriesSourceFile;
            
            if verbose
                fprintf('Mapping object to %s ... \n',sourceFileInfo.name)
            end
            
            % Load obj.Keys
            if obj.IndexedByKeys
                fid = fopen(obj.IndexFile,'rb');
                fseek(fid,int64(obj.OffsetToKeys),'bof');
                keys_tmp = textscan(fid,'%s',obj.NumEntriesSourceFile);
                fclose(fid);
                [keys_tmp,sortedIdx] = sort(keys_tmp{1});
                sortedIdx = uint32(sortedIdx);
                % Finds if contiguous strings are not equal
                lv = strdiff(keys_tmp);
                keysAreUnique = all(lv);
                if keysAreUnique
                    obj.KeyOrder = [];
                    obj.Keys = containers.Map(keys_tmp,sortedIdx);
                else % fake a multimap
                    % retain only the unique keys
                    keys_tmp = keys_tmp(lv);
                    % save the original positions
                    obj.KeyOrder = sortedIdx;
                    % Find index to unique keys
                    gg = uint32(find(lv));
                    % Finds how many times each key is repeated
                    rr = uint32(diff([0;gg]));
                    % Check if the index and cardinality can be packed into
                    % a uint32 if not a uint64 is used:
                    if (max(gg)<16777216) && (max(rr)<256)
                       % Codes the number of times the key is repeated into
                       % the most significant byte of a uint32 and the
                       % position to obj.KeyOrder in the three less
                       % significant bytes: 
                       gg = gg + (rr-1).*(16777216);
                    else
                       % Codes the number of times the key is repeated into
                       % the most significant four bytes of a uint64 and
                       % the position to obj.KeyOrder in the four less
                       % significant bytes:  
                       gg = uint64(gg) + uint64(rr-1)*4294967296;
                    end
                    obj.Keys = containers.Map(keys_tmp,gg);
                end
            else
                obj.KeyOrder = [];
                obj.Keys = [];
            end
            
            % The set method of obj.MemoryMappedIndex updates obj.Index if
            % necessary (i.e. load the Index into memory) 
            obj.MemoryMappedIndex = memMapIndex;
 
            % The set method of obj.FeatureInMemory updates obj.FeatureDara if
            % necessary (i.e. load the features into memory)             
            obj.FeatureInMemory   = featureInMemory;

            if verbose
                fprintf('Done. \n')
            end
            
        end % IndexedFile
        %==================================================================
        % set.Interpreter
        %==================================================================
        function obj = set.Interpreter(obj, h)
            if isa(h,'function_handle')
                obj.Interpreter = h;
            elseif isempty(h)
                obj.Interpreter = [];
            else
                error(message('bioinfo:BioIndexedFile:BioIndexedFile:InvalidInterpreter'))
            end
        end
        %==================================================================
        % set.MemoryMappedIndex
        %==================================================================
        function obj = set.MemoryMappedIndex(obj, tf)
            try
                tf = bioinfoprivate.opttf(tf,'MemoryMappedIndex',mfilename);
            catch ME
                if strcmp(ME.identifier,'bioinfo:BioIndexedFile:MemoryMappedIndexOptionNotLogical')
                    error(message('bioinfo:BioIndexedFile:BioIndexedFile:MemoryMappedIndex'))
                else
                    rethrow(ME)
                end
            end
            obj.PrivateMemoryMappedIndex = tf;
            if tf
                 obj.Index = [];
                 obj.ContiguousEntries = obj.ContiguousEntriesSourceFile;
            else
                % Load index with file ofstes into memory:
                if isempty(obj.SubsetIndex) % Object is NOT a subset
                    obj.Index = readWholeIndex(obj);
                    obj.ContiguousEntries = obj.ContiguousEntriesSourceFile;
                else % Object is a subset
                    i = obj.SubsetIndex;
                    % When it is a subset and indices are loaded into
                    % memory, we need to force to ContiguousEntries = false
                    obj.ContiguousEntries = false;
                    if obj.ContiguousEntriesSourceFile
                        obj.Index = readIndex(obj,reshape([i i+1]',2*obj.NumEntries,1));
                    else
                        obj.Index = readIndex(obj,reshape([2*i-1 2*i]',2*obj.NumEntries,1));
                    end
                end
            end
        end 
        %==================================================================
        % get.MemoryMappedIndex
        %==================================================================
        function tf = get.MemoryMappedIndex(obj)
            tf = obj.PrivateMemoryMappedIndex;
        end
        %==================================================================
        % set.FeatureInMemory
        %==================================================================
        function obj = set.FeatureInMemory(obj, tfs)
            % Validate a cell array with true/false flags
            if ~iscell(tfs) || (~isempty(tfs)&&~isvector(tfs))
                  error(message('bioinfo:BioIndexedFile:BioIndexedFile:InvalidFeatureInMemory'))
            end
            if numel(obj.FeatureOffset) ~= numel(tfs)
                error(message('bioinfo:BioIndexedFile:BioIndexedFile:BadSizeFeatureInMemory'))
            end
            try
                for i = 1:numel(tfs)
                   tfs{i} = bioinfoprivate.opttf(tfs{i},'InvalidElementFeatureInMemory','BioIndexedFile:BioIndexedFile'); 
                end
            catch ME
                if strcmp(ME.identifier,'bioinfo:BioIndexedFile:FeatureInMemoryOptionNotLogical')
                    error(message('bioinfo:BioIndexedFile:BioIndexedFile:FeatureInMemory'))

                else
                    rethrow(ME)
                end
            end
            
            if isempty(tfs) % This object does not have any indexable feature
                obj.PrivateFeatureInMemory = {};
                return
            end
            
            % Check what is the current status of features loaded in
            % memory:
            cfs = obj.FeatureInMemory;
            if isempty(cfs) % It is the very first time obj.FeatureInMemory is set
                cfs = repmat({false},1,numel(tfs));
            end
            
            % Load/unload individual features
            for i = 1:numel(tfs)
                % Load/unload only if the requested flag is different to
                % what is already stored:
                if tfs{i} ~= cfs{i}
                    if tfs{i} % Load
                        if isempty(obj.SubsetIndex) % Object is NOT a subset
                            obj.FeatureData{i} = getFeature(obj,obj.FeatureName{i});
                        else % Object is a subset
                            obj.FeatureData{i} = getFeature(obj,obj.FeatureName{i});
                        end
                    else % Unload
                        obj.FeatureData{i} = [];
                    end
                end
            end
            obj.PrivateFeatureInMemory = tfs;
        end
        %==================================================================
        % get.FeatureInMemory
        %==================================================================
        function tf = get.FeatureInMemory(obj)
            tf = obj.PrivateFeatureInMemory;
        end
        
        %==================================================================
        % getEntryByIndex
        %==================================================================
        function str = getEntryByIndex(obj,i)
            %GETENTRYBYINDEX Retrieve entries using a numeric index.
            %
            %  GETENTRYBYINDEX(OBJ,I) Extracts and concatenates the entries
            %  from the source file indexed by I. I is a numeric vector
            %  with positive integers less than or equal to the number of
            %  entries indexed by the object OBJ. There is a one-to-one
            %  relationship between the values in I and the output, order
            %  and quantity are preserved, despite repeated values in I.
            %  However, access to the source file is as efficient as
            %  possible by minimizing the number of file seeks and read
            %  operations. GETENTRYBYINDEX returns a string.
            %
            %  See also BioIndexedFile, BioIndexedFile.BioIndexedFile,
            %           BioIndexedFile.getEntryByKey, BioIndexedFile.getIndexByKey,
            %           BioIndexedFile.getKeys, BioIndexedFile.read
            
            if ~isnumeric(i) || any(i<1) || any(i>obj.NumEntries) || any(rem(i,1))
                error(message('bioinfo:BioIndexedFile:getEntryByIndex:nonNumeric'))
            end
            
            if isscalar(i)
                range = index2range(obj,i);
                fid = fopen(obj.InputFile,'rb');
                fseek(fid,int64(range(1)),'bof');
                str = fread(fid,diff(range),'*char')';
                fclose(fid);
            elseif isempty(i)&&~ischar(i)
                str = '';
            else
                i = i(:);
                ranges = index2range(obj,i);
                
                fid = fopen(obj.InputFile,'rb');
                
                % Check if within i the last entry is queried, and if so
                % check that such entry has a line feed character at the
                % end of the entry, otherwise when the string of the last
                % entry is concatenated with other the result may be
                % invalid for an interpreter:
                if isempty(obj.SubsetIndex)
                    e = find((i==obj.NumEntriesSourceFile),1);
                else
                    e = find((obj.SubsetIndex(i)==obj.NumEntriesSourceFile),1);
                end
                fixLineFeedInLastEntry = false;
                if e
                    fseek(fid,int64(ranges(e,2)-1),'bof');
                    if fread(fid,1,'*char')~=10
                        ranges(e,2) = ranges(e,2)+1;
                        fixLineFeedInLastEntry = true;
                    end
                end
                
                sizes = diff(ranges,[],2);
                strOffsets = [0;uint32(cumsum(double(sizes)))];
                str(1,strOffsets(end)) = ' ';
                
                % Access the file uniquely and orderly for performance
                [~,h,g] = unique(ranges(:,1),'last');
                uranges = ranges(h,:);
               
                % Figure out contiguous entries that can be read at once
                % and merge them ( f(:)==true indicates the limits of the
                % blocks ):
                f = uranges(2:end,1)~=uranges(1:end-1,2);

                blocksize = 2^22;

                % Split blocks in the contiguous ranges that are larger
                % than block size:
                cs = cumsum(double(diff(uranges,[],2))); 
                ft = [false;f];
                m = 0;
                s = find((cs>=blocksize-m) | ft,1);
                while ~isempty(s)
                    ft(s) = false; 
                    f(s-1) = true;
                    m = cs(s-1);
                    s = find((cs>=blocksize+m) | ft,1);
                end
                
                uranges = [uranges([true;f],1),uranges([f;true],2)];
                cui = [find([true;f]),find([f;true])];

                % Read unique ranges and copy them to output str
                
                m = 0;
                for j = 1:size(uranges,1)
                    
                    if j>m
                        o = uranges(j,1);
                        fseek(fid,int64(o),'bof');
                        m = find(uranges(:,2) <= (o+blocksize),1,'last');
                        n = uranges(m,2) - o;
                        if (fixLineFeedInLastEntry && m==size(uranges,1))
                            tmpBlock = [fread(fid,n-1,'*char');char(10)];
                        else
                            tmpBlock = fread(fid,n,'*char');
                        end
                    end
                    
                    si = h(cui(j,1):cui(j,2));
                    if numel(si)==1
                        str(strOffsets(si)+1:strOffsets(si+1)) = tmpBlock(uranges(j,1)-o+1:uranges(j,2)-o)';
                    else
                        strt =  tmpBlock(uranges(j,1)-o+1:uranges(j,2)-o)';
                        p = 1;
                        for k = 1:numel(si)
                            r = [strOffsets(si(k))+1 strOffsets(si(k)+1)];
                            q = p + diff(r);
                            str(r(1):r(2)) = strt(p:q);
                            p = q + 1;
                        end
                    end
                end
                fclose(fid);
                
                % Copy non-unique entries
                g(h) = 0;
                f = find(g); % copy from ...
                g = h(nonzeros(g)); % copy to ...
                for k = 1:numel(g)
                    r = [strOffsets(f(k))+1 strOffsets(f(k)+1)];
                    q = [strOffsets(g(k))+1 strOffsets(g(k)+1)];
                    str(r(1):r(2)) = str(q(1):q(2));
                end
                
                % For FLAT files we need to interleave the entrySeparator
                if strcmp(obj.IndexerType,'flat')
                    delsize = numel(obj.EntryDelimiter)+1;
                    strtmp(1,delsize.*(numel(sizes)-1)+numel(str)) = ' ';
                    k1 = 1; l1 = 1;
                    for j = 1:numel(sizes)-1
                        k2 = k1+sizes(j);
                        l2 = l1+sizes(j)+delsize;
                        strtmp(l1:l2-1) = [str(k1:k2-1) sprintf('%s\n',obj.EntryDelimiter)];
                        k1 = k2;
                        l1 = l2;
                    end
                    j = j+1;
                    k2 = k1+sizes(j);
                    l2 = l1+sizes(j);
                    strtmp(l1:l2-1) = str(k1:k2-1);
                    str = strtmp;
                end
            end
        end %getEntryByIndex
        %==================================================================
        % getKeys
        %==================================================================
        function keys = getKeys(obj)
            %GETKEYS Retrieve alphanumeric keys.
            %
            %  GETKEYS(OBJ) Returns a cell array containing all of the
            %  alphanumeric keys stored, in the same order as they appear
            %  in the source file, even if they are not unique.
            %
            %  See also BioIndexedFile, BioIndexedFile.BioIndexedFile,
            %           BioIndexedFile.getEntryByIndex, BioIndexedFile.getEntryByKey,
            %           BioIndexedFile.getIndexByKey,BioIndexedFile.read
            
            if ~obj.IndexedByKeys
                error(message('bioinfo:BioIndexedFile:getIndexByKey:noKeys'))
            end
            keys = (obj.Keys.keys)';
            if obj.Keys.Count~=obj.NumEntries % Non unique keys
                v = cell2mat(obj.Keys.values);
                if isa(v,'uint32')
                    w = v./16777216;
                    r = v-w.*16777216;
                else
                    w = v./4294967296;
                    r = v-w.*4294967296;
                end
                h = zeros(obj.NumEntries,1,'uint32');
                h(r) = 1;
                h = uint32(cumsum(single(h)))+1;
                h(r) = h(r)-1;
                h(obj.KeyOrder) = h;
                keys = keys(h);
            else
                keys(cell2mat(obj.Keys.values))=keys;
            end
        end %getKeysByIndex
        %==================================================================
        % getIndexByKey
        %==================================================================
        function [idx,otoli] = getIndexByKey(obj,key)
            %GETINDEXBYKEY Indices of entries for alphanumeric keys.
            %
            %  GETINDEXBYKEY(OBJ,KEY) Returns the indices of entries in the
            %  source file that match alphanumeric keys. KEY is a string or
            %  a cell string with multiple keys. If keys are not unique in
            %  the source file, GETINDEXBYKEY returns all indices of
            %  entries that match a given key at its respective
            %  position accordingly with the query. Therefore, when keys
            %  are unique in the source file there is an one-to-one
            %  relationship between the strings in KEY and the output.
            %
            %  [I,J] = GETINDEXBYKEY(OBJ,KEY) Returns a logical vector that
            %  selects only the last match for every key, such that there
            %  is an one-to-one relationship between the strings in KEY and
            %  I(J).
            %
            %  See also BioIndexedFile, BioIndexedFile.BioIndexedFile,
            %           BioIndexedFile.getEntryByIndex, BioIndexedFile.getEntryByKey,
            %           BioIndexedFile.getKeys, BioIndexedFile.read
            if ~obj.IndexedByKeys
                error(message('bioinfo:BioIndexedFile:getIndexByKey:noKeys'))
            end
            if ischar(key)
                key = {key};
            end
            if ~iscellstr(key)
                error(message('bioinfo:BioIndexedFile:getIndexByKey:invalidKey'))
            end
            n = numel(key);
            
            if obj.Keys.Count==obj.NumEntries % Unique keys 
                idx = zeros(n,1,'uint32');
                for i = 1:n
                    idx(i) = obj.Keys(key{i});
                end
                if nargout>1
                    otoli = true(n,1);
                end
            else % Non-unique keys in the source file
                idx = cell(n,1);
                for i = 1:n
                    v = obj.Keys(key{i});
                    if isa(v,'uint32')
                        rl = v./16777216;
                        v = v - rl.*16777216;
                    else
                        rl = v./4294967296;
                        v = v - rl.*4294967296;
                        rl = uint32(rl);
                        v = uint32(v);
                    end
                    idx{i} = obj.KeyOrder(v-rl:v);
                end
                if nargout>1
                    rl = uint32(cumsum(single(cellfun(@numel,idx))));
                    otoli = false(rl(end),1);
                    otoli(rl) = true;
                end
                idx = cell2mat(idx);
            end
            
            
        end
        %==================================================================
        % getEntryByKey
        %==================================================================
        function str = getEntryByKey(obj,key)
            %GETENTRYBYKEY Retrieve entries using an alphanumeric key.
            %
            %  GETENTRYBYKEY(OBJ,KEY) Extracts and concatenates the entries
            %  from the source file that are identified by KEY. KEY is a
            %  string or a cell string with multiple keys. If keys are not
            %  unique in the source file, all entries that match a given
            %  key are concatenated into the output at its respective
            %  position accordingly with the query. Therefore, when keys
            %  are unique in the source file, there is a one-to-one
            %  relationship between the strings in KEY and the output.
            %  Access to the source file is as efficient as possible by
            %  minimizing the number of file seeks and read operations.
            %  GETENTRYBYKEY returns a string.
            %
            %  See also BioIndexedFile, BioIndexedFile.BioIndexedFile,
            %           BioIndexedFile.getEntryByIndex, BioIndexedFile.getIndexByKey,
            %           BioIndexedFile.getKeys, BioIndexedFile.read
            
            idx = getIndexByKey(obj,key);
            str = getEntryByIndex(obj,idx);
            
        end
        %==================================================================
        % read
        %==================================================================
        function data = read(obj,s)
            %READ Read one or more entries from the source file.
            %
            %  READ(OBJ,I) Reads the entries from the source file indexed
            %  by I using the file parser indicated by the object property
            %  Interpreter. The object property Interpreter must contain a
            %  handle to a function that accepts strings containing the
            %  data to be parsed. If the object property Interpreter is
            %  empty, then only the concatenated string with the
            %  corresponding entries is returned. I is a numeric vector
            %  with positive integers less than or equal to the number of
            %  entries indexed by the object OBJ. There is a one-to-one
            %  relationship between the values in I and the output. Order
            %  and quantity are preserved, despite repeated values in I.
            %
            %  READ(OBJ,KEY) Reads the entries from the source file that
            %  are identified by KEY. KEY is a string or a cell string with
            %  multiple keys. If keys are not unique in the source file,
            %  all entries that match a given key are included in the
            %  output at its respective position according to the query.
            %  Therefore, when keys are unique in the source file there is
            %  a one-to-one relationship between the strings in KEY and the
            %  output.
            %
            %  See also BioIndexedFile, BioIndexedFile.BioIndexedFile,
            %           BioIndexedFile.getEntryByIndex, BioIndexedFile.getEntryByKey,
            %           BioIndexedFile.getIndexByKey, BioIndexedFile.getKeys
            if ischar(s)&&isvector(s)
                str = getEntryByKey(obj,s);
            elseif isnumeric(s)&&isvector(s)
                str = getEntryByIndex(obj,s);
            elseif iscellstr(s)
                str = getEntryByKey(obj,s);
            else
                error(message('bioinfo:BioIndexedFile:read:invalidInput'))
            end
            if ~isempty(obj.Interpreter)
                data = obj.Interpreter(str);
            else
                data = str;
            end
        end
        %==================================================================
        % getSubset
        %==================================================================
        function obj = getSubset(obj,s)
            %GETSUBSET Returns an object that access a subset of the entries in the source file.
            %
            %  GETSUBSET(OBJ,I) Returns an object that accesses a subset of
            %  the entries in the source file. Entries for the subset are
            %  selected by I, a vector with unique positive integers which
            %  value is less than or equal to the number of entries indexed
            %  by the object OBJ. There is a one-to-one relationship
            %  between the number and order in I and the order in the
            %  entries indexed by the new object.
            %
            %  GETSUBSET(OBJ,KEY) Returns an object that accesses a subset of
            %  the entries in the source file. Entries for the subset are
            %  selected by KEY, a string or cell array of unique strings
            %  specifying multiple keys. If keys are not unique in the
            %  source file, all entries that match a given key are indexed
            %  in the new object, therefore a one-to-one relationship is
            %  not maintained. When keys are unique in the input object
            %  there is a one-to-one relationship between the number and
            %  order of strings in KEY and the entries indexed by the new
            %  object.
            %
            %  See also BioIndexedFile, BioIndexedFile.BioIndexedFile,
            %           BioIndexedFile.getEntryByIndex, BioIndexedFile.getEntryByKey,
            %           BioIndexedFile.getIndexByKey, BioIndexedFile.getKeys
            %           BioIndexedFile.read


            if ischar(s)&&isvector(s)
                s = {s};
            end

            if isnumeric(s)&&isscalar(s)
                idx = uint32(s);
            elseif isnumeric(s)&&isvector(s)
                if issorted(s)
                    if any(diff(s)==0)
                        error(message('bioinfo:BioIndexedFile:getSubset:invalidNumericInput'))
                    end
                elseif numel(unique(s)) ~= numel(s)
                    error(message('bioinfo:BioIndexedFile:getSubset:invalidNumericInput'))
                end
                idx = uint32(s(:));
            elseif iscellstr(s)
                if numel(unique(s)) ~= numel(s)
                    error(message('bioinfo:BioIndexedFile:getSubset:invalidAlphanumericInput'))
                end
                idx = getIndexByKey(obj,s);
            else
                error(message('bioinfo:BioIndexedFile:getSubset:invalidInput'))
            end
            if ~isempty(idx) && any(idx)==0
                    error(message('bioinfo:BioIndexedFile:getSubset:logicalIndex'))
            end
            if any(idx>obj.NumEntries)
                    error(message('bioinfo:BioIndexedFile:getSubset:outOfBounds'))
            end
            if isempty(obj.SubsetIndex)
                obj.SubsetIndex = idx;
            else
                obj.SubsetIndex = obj.SubsetIndex(idx);
            end
            
            NumEntriesNew = numel(idx);
            
            % Check if obj.Index needs to be updated (is loaded in memory)
            % so it will be updated to contain only the indices of the subset:
            if ~obj.MemoryMappedIndex
                if obj.ContiguousEntries
                    obj.Index = reshape(obj.Index([idx idx+1]'),2*NumEntriesNew,1);
                    obj.ContiguousEntries = false;
                else
                    obj.Index = reshape(obj.Index([(idx-1)*2+1 idx*2]'),2*NumEntriesNew,1);
                end
            end
            
            if obj.IndexedByKeys
                keys_tmp = getKeys(obj);
                keys_tmp = keys_tmp(idx);
                [keys_tmp,sortedIdx] = sort(keys_tmp);
                
                % Finds if contiguous strings are not equal
                lv = strdiff(keys_tmp);
                keysAreUnique = all(lv);
                if keysAreUnique
                    obj.KeyOrder = [];
                    obj.Keys = containers.Map(keys_tmp,sortedIdx);
                else % fake a multimap
                    % retain only the unique keys
                    keys_tmp = keys_tmp(lv);
                    % save the original positions
                    obj.KeyOrder = sortedIdx;
                    % Find index to unique keys
                    gg = uint32(find(lv));
                    % Finds how many times each key is repeated
                    rr = uint32(diff([0;gg]));
                    % Check if the index and cardinality can be packed into
                    % a uint32 if not a uint64 is used:
                    if (max(gg)<16777216) && (max(rr)<256)
                       % Codes the number of times the key is repeated into
                       % the most significant byte of a uint32 and the
                       % position to obj.KeyOrder in the three less
                       % significant bytes: 
                       gg = gg + (rr-1).*(16777216);
                    else
                       % Codes the number of times the key is repeated into
                       % the most significant four bytes of a uint64 and
                       % the position to obj.KeyOrder in the four less
                       % significant bytes:  
                       gg = uint64(gg) + uint64(rr-1)*4294967296;
                    end
                    obj.Keys = containers.Map(keys_tmp,gg);
                end
            end
            
            % Check if obj.FeatureData needs to be update (Features loaded
            % in memory) so it will be updated to contain only the data of
            % the subset: 
            cfs = obj.FeatureInMemory;
            for i = 1:numel(cfs)
                if cfs{i}
                    data = obj.FeatureData{i};
                    obj.FeatureData{i} = data(idx);
                end
            end
            
            obj.NumEntries = NumEntriesNew;
            
            if obj.NumEntries==0
                obj.FeatureInMemory(:) = {1};
            end
                
        end
        %==================================================================
        % getDictionary
        %==================================================================
        function dic = getDictionary(obj)
            %GETDICTIONARY Retrieve the dictionary of a SAM indexed file
            %
            %  DIC = GETDICTIONARY(OBJ) Returns a cell string with the
            %  names of the references used in a SAM indexed file. If OBJ
            %  is not the index of a SAM formatted file, GETDICTIONARY
            %  returns an error.
            if ~strcmp(obj.FileFormat,'sam')
                error(message('bioinfo:BioIndexedFile:getDictionary:invalidFileFormat'))
            end
            [u,h] = unique(obj.getFeature('Reference'));
            % remove unmapped entries
            h(u==0)=[];
            % create output cellstr
            dic = cell(numel(h),1);
            for i = 1:numel(h)
                tmp = samread(getEntryByIndex(obj,h(i)));
                dic{i} = tmp.ReferenceName;
            end
        end
    end % methods
    
    methods  (Hidden = true)
        %==================================================================
        % getFeature
        %==================================================================
        function [data,obj] = getFeature(obj,name,s)
            %GETFEATURE Query a feature from the index file
            %
            %  GETFEATURE(OBJ,NAME) Returns the data corresponding to the
            %  feature NAME. First it tries to query the index file, but if
            %  the feature is not indexed yet GETFEATURE calls INDEXFEATURE
            %  to update the index file from the source file.
            %  When an index file is just created, the first query may be
            %  slow, however, subsequent queries will be more efficient
            %  since the data is brought from contiguous file positions.
            %
            %  GETFEATURE(OBJ,NAME,I) Returns the data corresponding to the
            %  feature NAME and indexed by the numeric vector I. I contains
            %  positive integers less than or equal to the number of
            %  entries indexed by the object. There is a one-to-one
            %  relationship between the values in I and the output. Order
            %  and quantity are preserved, despite repeated values in I.
            %
            %  GETFEATURE(OBJ,NAME,KEY) Returns the data corresponding to the
            %  feature NAME and identified by KEY. KEY is a string or a
            %  cell string with multiple keys. If keys are not unique in
            %  the source file, all entries that match a given key are
            %  included in the output at its respective position according
            %  to the query. Therefore, when keys are unique in the source
            %  file there is a one-to-one relationship between the strings
            %  in KEY and the output.
            %
            %  [DATA,OBJ] = GETFEATURE(...) Returns the updated object in
            %  which the property FeatureOffset reflects the respective
            %  features that are indexed by the index file by having an
            %  offset different from 0.
            
            h = findFeature(obj,name);
            
            if (obj.NumEntries==0)
                data = zeros(0,1,obj.FeatureType{h});
                return
            end
                
            % Shortcut for cases where the features are already loaded in
            % memory:
            if numel(obj.FeatureInMemory)>0 && obj.FeatureInMemory{h}
                if nargin==2 % All entries indexed by the object are requested:
                    data = obj.FeatureData{h}; 
                else % A subset of the entries indexed by the object are requested:
                    if ischar(s)&&isvector(s)
                        i = getIndexByKey(obj,s);
                    elseif isnumeric(s)&&isvector(s)
                        if any(s<1) || any(s>obj.NumEntries) || any(rem(s,1))
                            error(message('bioinfo:BioIndexedFile:getFeature:nonNumeric'))
                        else
                            i = s(:);
                        end
                    elseif iscellstr(s)
                        i = getIndexByKey(obj,s);
                    else
                        error(message('bioinfo:BioIndexedFile:getFeature:invalidInput'))
                    end
                    data = obj.FeatureData{h}(i);
                end
                return
            end
            
            if obj.FeatureOffset(h)==0
                obj = indexFeature(obj,name);
            end
            
            if nargin==2
                % All entries indexed by the object are requested:
                if isempty(obj.SubsetIndex)
                    % If the object is not a subset then we can just read
                    % all the features in the file at once
                    fid = fopen(obj.IndexFile,'r');
                    fseek(fid,int64(obj.FeatureOffset(h)),'bof');
                    data = fread(fid,obj.NumEntriesSourceFile,[obj.FeatureType{h} '=>' obj.FeatureType{h}]);
                    fclose(fid);
                    return
                else
                    % Since the object is a subset we need need to read
                    % only the entries in the source file that are relative
                    % to the current subset:
                    i = obj.SubsetIndex;
                end
            else
                % A subset of the entries indexed by the object are requested:
                if ischar(s)&&isvector(s)
                    i = getIndexByKey(obj,s);
                elseif isnumeric(s)&&isvector(s)
                    if any(s<1) || any(s>obj.NumEntries) || any(rem(s,1))
                        error(message('bioinfo:BioIndexedFile:getFeature:nonNumeric'))
                    else
                        i = s(:);
                    end
                elseif iscellstr(s)
                    i = getIndexByKey(obj,s);
                else
                    error(message('bioinfo:BioIndexedFile:getFeature:invalidInput'))
                end
            
                % At this point the index i is relative to the entries indexed
                % by the object (the object may be only a subset of all the
                % entries in the source file).
                if ~isempty(obj.SubsetIndex)
                    i = obj.SubsetIndex(i);
                end
            end % nargin==2
            
            fType = obj.FeatureType{h};
            if isscalar(i)
                switch fType
                    case 'uint8'
                        fTypeSize = 1;
                    case 'uint16'
                        fTypeSize = 2;
                    case 'uint32'
                        fTypeSize = 4;
                    otherwise
                        error(message('bioinfo:BioIndexedFile:getFeature:unsupportedFeatureType'))
                end
                fid = fopen(obj.IndexFile,'r');
                fseek(fid,int64(obj.FeatureOffset(h)+uint64(i-1)*fTypeSize),'bof');
                data = fread(fid,1,[fType '=>' fType]);
                fclose(fid);
            elseif isempty(i)&&~ischar(i)
                data = zeros(0,1,fType);
            else
                % Figure out the fewest number of file reads
                [u,~,i2u] = unique(i);
                fdu = find(diff(u)>1);
                bs = [0;fdu]; % block starts (zero-indexed)
                bn = [fdu;numel(u)]-bs; % number of elements per block
                if numel(fdu)>100 %if too many file reads, just read once and subindex
                    fid = fopen(obj.IndexFile,'r');
                    fseek(fid,int64(obj.FeatureOffset(h)),'bof');
                    data = fread(fid,obj.NumEntriesSourceFile,[fType '=>' fType]);
                    fclose(fid);
                    data = data(i);
                else
                    data = zeros(numel(u),1,fType);
                    switch fType
                        case 'uint8'
                            fTypeSize = 1;
                        case 'uint16'
                            fTypeSize = 2;
                        case 'uint32'
                            fTypeSize = 4;
                        otherwise
                            error(message('bioinfo:BioIndexedFile:getFeature:unsupportedFeatureType'))
                    end
                    fid = fopen(obj.IndexFile,'r');
                    for j = 1:numel(bs)
                        fseek(fid,int64(obj.FeatureOffset(h)+uint64(u(bs(j)+1)-1)*fTypeSize),'bof');
                        data(bs(j)+(1:bn(j))) = fread(fid,bn(j),[fType '=>' fType]);
                    end
                    fclose(fid);
                    data = data(i2u); % replicate and re-organize as input index
                end
            end
        end
        %==================================================================
        % indexFeature
        %==================================================================
        function obj = indexFeature(obj,name)
            %INDEXFEATURE Updates the index file with feature data.
            %
            %  INDEXFEATURE(OBJ,NAME) Updates an index file with the data
            %  for feature NAME.
            %
            %  OBJ = INDEXFEATURE(...) Returns the updated object in which
            %  the property FeatureOffset reflects the respective
            %  features that are indexed by the index file by having an
            %  offset different from 0.
            
            
            % Check if the feature is already indexed but the object
            % instance is not updated:
            obj = obj.readHeaderInfo(false);
            h = findFeature(obj,name);
            if obj.FeatureOffset(h)>0
                return
            end
            
            % Figure out the offset where the feature will be appended to
            % the index file, which is the file size:
            d = dir(obj.IndexFile);
            offset = uint64(d.bytes);
            
            % Call the feature indexer
            try
                if obj.MemoryMappedIndex || ~isempty(obj.SubsetIndex)
                    BioIndexedFile.biofeatureindexermex(...
                        obj.FileFormat,...
                        obj.FeatureName{h},...
                        obj.FeatureType{h},...
                        obj.InputFile,...
                        obj.IndexFile,...
                        obj.ContiguousEntriesSourceFile,...
                        uint32(obj.NumEntriesSourceFile),...
                        obj.IndexType,...
                        readWholeIndex(obj))
                else
                    BioIndexedFile.biofeatureindexermex(...
                        obj.FileFormat,...
                        obj.FeatureName{h},...
                        obj.FeatureType{h},...
                        obj.InputFile,...
                        obj.IndexFile,...
                        obj.ContiguousEntriesSourceFile,...
                        uint32(obj.NumEntriesSourceFile),...
                        obj.IndexType,...
                        obj.Index)
                end
            catch ME
                bioinfoprivate.bioclsrethrow(mfilename,'indexFeature',ME)
            end
            
            d = dir(obj.IndexFile);
            if (offset == uint64(d.bytes))
                error(message('bioinfo:BioIndexedFile:indexFeature:indexerFailure', name, obj.FileFormat))
            else
                % Update header information and object
                obj = updateHeaderInfo(obj,name,offset);
            end
        end
        
    end % methods (Hidden = true)
    
    methods (Access = private)

        %==================================================================
        % readWholeIndex
        %==================================================================
        function idx = readWholeIndex(obj)
            %READWHOLEINDEX Read file offsets for all entries
            fid = fopen(obj.IndexFile,'r');
            if fid==-1
                error(message('bioinfo:BioIndexedFile:readWholeIndex:CannotOpenFile', obj.IndexFile))
            end
            fseek(fid,int64(obj.OffsetToIndex),'bof');
            if obj.IndexType==1
                if obj.ContiguousEntries
                    idx = fread(fid,obj.NumEntriesSourceFile+1,'uint64=>uint64');
                else
                    idx = fread(fid,obj.NumEntriesSourceFile*2,'uint64=>uint64');
                end
            else % ==0 or []
                if obj.ContiguousEntries
                    idx = fread(fid,obj.NumEntriesSourceFile+1,'uint32=>uint32');
                else
                    idx = fread(fid,obj.NumEntriesSourceFile*2,'uint32=>uint32');
                end
            end
            
            fclose(fid);
        end
        %==================================================================
        % readIndex
        %==================================================================
        function idx = readIndex(obj,i)
            %READINDEX Read file offsets for entries in i
            if isempty(i)
                if obj.IndexType==1
                    idx = zeros(0,2,'uint64');
                else % ==0 or []
                    idx = zeros(0,2,'uint32');
                end
                return
            end
            % Figure out the fewest number of file reads
            [u,~,i2u] = unique(i);
            fdu = find(diff(u)>1);
            bs = [0;fdu]; % block starts (zero-indexed)
            bn = [fdu;numel(u)]-bs; % number of elements per block
            if numel(fdu)>100 %if too many file reads, just read once and subindex
                idx = readWholeIndex(obj);
                idx = idx(i);
            else
                fid = fopen(obj.IndexFile,'r');
                if fid==-1
                    error(message('bioinfo:BioIndexedFile:readIndex:CannotOpenFile', obj.IndexFile))
                end
                if obj.IndexType==1
                    idx = zeros(numel(u),1,'uint64');
                    for j = 1:numel(bs)
                        fseek(fid,int64(obj.OffsetToIndex+uint64(u(bs(j)+1)-1)*8),'bof');
                        idx(bs(j)+(1:bn(j))) = fread(fid,bn(j),'uint64=>uint64');
                    end                    
                else % ==0 or []
                    idx = zeros(numel(u),1,'uint32');
                    for j = 1:numel(bs)
                        fseek(fid,int64(obj.OffsetToIndex+uint64(u(bs(j)+1)-1)*4),'bof');
                        idx(bs(j)+(1:bn(j))) = fread(fid,bn(j),'uint32=>uint32');
                    end
                end
                idx = idx(i2u); % replicate and re-organize as input index
                fclose(fid);
            end
            idx = reshape(idx,size(i));
        end
        %==================================================================
        % findFeature
        %==================================================================
        function h = findFeature(obj,name)
            %FINDFEATURE Find the index in property FeatureName to a feature.
            %
            % H = FINDFEATURE(OBJ,NAME) Finds the index in property
            % FeatureName to a feature and throws an error if absent.
            h = find(strcmp(obj.FeatureName,name),1);
            if isempty(h)
                error(message('bioinfo:BioIndexedFile:findFeature:invalidFeature', name, obj.FileFormat))
            end
        end
        %==================================================================
        % readHeaderInfo
        %==================================================================
        function obj = readHeaderInfo(obj,wf)
            %READHEADERINFO Read header info from the index file.
            %
            % OBJ = READHEADERINFO(OBJ) Updates the properties in the
            % object instance OBJ with information from the header in the
            % index file.
            %
            % OBJ = READHEADERINFO(OBJ,true) Warns inconsistencies.
            
            if nargin==1
                wf = false;
            end
            fid = fopen(obj.IndexFile,'rb');
            if fid==-1
                error(message('bioinfo:BioIndexedFile:readHeaderInfo:CannotOpenFile', obj.IndexFile))
            end
            c = onCleanup(@()fclose(fid));
            Magic                 = fread(fid,7,'uint8=>char')';           
            FormatVersion         = sscanf((fread(fid,4,'uint8=>char')'),'%f');
            if ~strncmp(Magic,sprintf('MBTidx'),6) || (FormatVersion~=1.0 && FormatVersion~=2.0)
                error(message('bioinfo:BioIndexedFile:readHeaderInfo:InvalidVersion', obj.IndexFile))
            end
            LockedFile            = ~(fread(fid,1,'uint8=>double') == 0);       %#ok<NASGU>
            obj.NumEntriesSourceFile = fread(fid,1,'uint32=>double');
            if FormatVersion==2.0
                obj.OffsetToKeys      = fread(fid,1,'uint64=>double');
                obj.OffsetToIndex     = fread(fid,1,'uint64=>double');
                if wf && obj.IndexedByKeys && (obj.OffsetToKeys==0)
                    warning(message('bioinfo:BioIndexedFile:BioIndexedFile:noKeysInIndexFile'))
                    obj.IndexedByKeys = false;
                end
                obj.ContiguousEntriesSourceFile = fread(fid,1,'uint8=>double');
                if obj.ContiguousEntriesSourceFile>1
                    error(message('bioinfo:BioIndexedFile:readHeaderInfo:InvalidIndexFormat'))
                else
                    obj.ContiguousEntriesSourceFile = ~(obj.ContiguousEntriesSourceFile == 0);
                end
                obj.IndexType = fread(fid,1,'uint8=>uint8');
                if obj.IndexType>1
                    error(message('bioinfo:BioIndexedFile:readHeaderInfo:InvalidIndexType'))
                end
            else % FormatVersion==1.0
                obj.OffsetToKeys      = fread(fid,1,'uint16=>double');
                obj.OffsetToIndex     = fread(fid,1,'uint32=>double');
                if wf && obj.IndexedByKeys && (double(obj.OffsetToKeys)==double(obj.OffsetToIndex))
                    warning(message('bioinfo:BioIndexedFile:BioIndexedFile:noKeysInIndexFile'))
                    obj.IndexedByKeys = false;
                end
                obj.ContiguousEntriesSourceFile = ~(fread(fid,1,'uint8=>double') == 0);
                obj.IndexType = uint8(0); % 1.0 could only handle 32 bit
            end
            
            obj.EntryDelimiter    = fscanf(fid,'%s',1);
            fseek(fid,1,'cof'); %jumps the tab after the EntryDelimiter
            offsetToNumFeatures   = ftell(fid);
            numFeatures           = fread(fid,1,'uint8=>double');
            offsetToFeatures      = offsetToNumFeatures +1;
            for i=1:numFeatures
                fseek(fid,offsetToFeatures,'bof');
                name   = fscanf(fid,'%s\t',1);
                offsetToFeatures = offsetToFeatures + numel(name) + 1;
                fseek(fid,offsetToFeatures,'bof');
                type   = fscanf(fid,'%s\t',1);
                offsetToFeatures = offsetToFeatures + numel(type) + 1;
                fseek(fid,offsetToFeatures,'bof');
                if FormatVersion==2.0
                    offset = fread(fid,1,'uint64');
                    offsetToFeatures = offsetToFeatures + 8;
                else % FormatVersion==1.0
                    offset = fread(fid,1,'uint32');
                    offsetToFeatures = offsetToFeatures + 4;
                end
                h = find(strcmp(name,obj.FeatureName),1);
                if isempty(h)
                    if wf
                        warning(message('bioinfo:BioIndexedFile:BioIndexedFile:InvalidFeatureName', name, obj.IndexFile, obj.FileFormat))
                    end
                elseif ~strcmp(obj.FeatureType{h},type)
                    if wf
                        warning(message('bioinfo:BioIndexedFile:BioIndexedFile:InvalidFeatureType', name, obj.IndexFile, type, obj.FeatureType{ h }))
                    end
                    obj.FeatureOffset(h) = 0; % zero is used as an indicator that such feature has not been indexed
                else
                    obj.FeatureOffset(h) = offset;
                end
            end
        end
        %==================================================================
        % updateHeaderInfo
        %==================================================================
        function obj = updateHeaderInfo(obj,name,offset)
            %UPDATEHEADERINFO Write feature info into header of index file.
            %
            % UPDATEHEADERINFO(OBJ,NAME,OFFSET) Appends the feature NAME,
            % its precision, and file position in the index file OFFSET to
            % the header info of the index file.
            %
            % OBJ = UPDATEHEADERINFO(...) Updates also the properties in the
            % object instance OBJ.
            
            fid = fopen(obj.IndexFile,'r+b');
            if fid==-1
                error(message('bioinfo:BioIndexedFile:updateHeaderInfo:CannotOpenFile', obj.IndexFile))
            end
            c = onCleanup(@()fclose(fid));
            Magic                 = fread(fid,7,'uint8=>char')';             
            FormatVersion         = sscanf((fread(fid,4,'uint8=>char')'),'%f');
            if ~strncmp(Magic,sprintf('MBTidx'),6) || (FormatVersion~=1.0 && FormatVersion~=2.0)
                error(message('bioinfo:BioIndexedFile:updateHeaderInfo:InvalidVersion', obj.IndexFile))
            end
            LockedFile            = ~(fread(fid,1,'uint8=>double') == 0);       %#ok<NASGU>
            obj.NumEntriesSourceFile        = fread(fid,1,'uint32=>double');
            
            if FormatVersion==2.0
                obj.OffsetToKeys      = fread(fid,1,'uint64=>double');
                obj.OffsetToIndex     = fread(fid,1,'uint64=>double');
                obj.ContiguousEntriesSourceFile = ~(fread(fid,1,'uint8=>double') == 0);
                obj.IndexType         = fread(fid,1,'uint8=>uint8');
            else % FormatVersion==1.0
                obj.OffsetToKeys      = fread(fid,1,'uint16=>double');
                obj.OffsetToIndex     = fread(fid,1,'uint32=>double');
                obj.ContiguousEntriesSourceFile = ~(fread(fid,1,'uint8=>double') == 0);
                obj.IndexType         = uint8(0); % 1.0 could only handle 32 bit  
            end
            obj.EntryDelimiter    = fscanf(fid,'%s',1);
            fseek(fid,1,'cof'); %jumps the tab after the EntryDelimiter
            offsetToNumFeatures   = ftell(fid);
            numFeatures           = fread(fid,1,'uint8=>double');
            offsetToNewFeature = offsetToNumFeatures + 1;
            for i=1:numFeatures % jump space of feature info already in header
                fseek(fid,offsetToNewFeature,'bof');
                name_tmp = fscanf(fid,'%s\t',1);
                offsetToNewFeature = offsetToNewFeature + numel(name_tmp) + 1;
                fseek(fid,offsetToNewFeature,'bof');
                type_tmp = fscanf(fid,'%s\t',1);
                offsetToNewFeature = offsetToNewFeature + numel(type_tmp) + 1;
                fseek(fid,offsetToNewFeature,'bof');
                if FormatVersion==2.0
                    fread(fid,1,'uint64');
                    offsetToNewFeature = offsetToNewFeature + 8;
                else % FormatVersion==1.0                    
                    fread(fid,1,'uint32');
                    offsetToNewFeature = offsetToNewFeature + 4;
                end
            end
            offsetToNewFeature    = ftell(fid);
            h = find(strcmp(name,obj.FeatureName),1); % update obj, when FeatureOffset~=0 we
            obj.FeatureOffset(h)  = uint64(offset);   % know that the Feature is indexed
            fseek(fid,offsetToNumFeatures,'bof');
            fwrite(fid,numFeatures+1,'uint8');
            fseek(fid,offsetToNewFeature,'bof');
            fprintf(fid,'%s\t%s\t',name,obj.FeatureType{h});
            if FormatVersion==2.0
                fwrite(fid,uint64(offset),'uint64');
            else % FormatVersion==1.0
                % This allows continue indexing new feature columns into an
                % index file which is version 1.0
                fwrite(fid,uint32(offset),'uint32');
            end
        end
        %==================================================================
        % index2range
        %==================================================================
        function range = index2range(obj,i)
            %INDEX2RANGE Range of an entry in the source file.
            %
            %  INDEX2RANGE(OBJ,I) Returns a [Nx2] numeric vector with the
            %  file positions (start and end) of the entries indexed by I.
            %  I is numeric vector with positive integers equal or smaller
            %  than the number of entries indexed by th object OBJ.
            
            if obj.ContiguousEntries
                if obj.MemoryMappedIndex
                    if ~isempty(obj.SubsetIndex)
                        i = obj.SubsetIndex(i);
                    end
                    range = readIndex(obj,[i i+1]);
                else
                    if ~isempty(obj.SubsetIndex)
                        error(message('bioinfo:BioIndexedFile:index2range:ContiguousAndSubset'))
                    end
                    range = obj.Index([i i+1]);
                end
            else
                if obj.MemoryMappedIndex
                    if ~isempty(obj.SubsetIndex)
                        i = obj.SubsetIndex(i);
                    end
                    range = readIndex(obj,[2*i-1 2*i]);
                else
                    range = obj.Index([2*i-1 2*i]);
                end
            end
        end
        
    end % methods (Access = private)
    
    methods (Static, Access = private)
        biofeatureindexermex(ia1,ia2,ia3,ia4,ia5,ia6,ia7,ia8,ia9)
        %BIOFEATUREINDEXERMEX Creates an index file from a source file
        %
        % Note: this is an undocumented private and static method that may
        % not be supported in future releases to the Bioinformatics Toolbox.
        %
        % INPUT ARGUMENTS:
        %
        %   1. File type: {'sam','fasta','fastq'}
        %   2. Feature name: {'Flag','Position','MappingQuality','Mateposition',
        %                     'MateIndex','Length','QueryLength'}
        %   3. Feature Type: {'uint8','uint16','uint32','char'}
        %   4. Source filename
        %   5. Index filename
        %   6. A logical that indicates if the entries are contiguous in the
        %   source file
        %   7. An uint32 with num of entries
        %   8. Index type (0 for uint32 and 1 for uint64)
        %   9. An array of uint32 with the offsets to the source file
        %
        
        [oa1,oa2] = bioindexermex(ia1,ia2,ia3,ia4,ia5,ia6,ia7,ia8,ia9,ia10,ia11,ia12)
        %BIOINDEXERMEX Creates an index file from a source file
        %
        % Note: this is an undocumented private and static method that may
        % not be supported in future releases to the Bioinformatics Toolbox.
        %
        % INPUT ARGUMENTS:
        %
        %   1. File type: {'sam','fasta','fastq','mrtab','table','flat','fixnl'}
        %   2. Source filename
        %   3. Index filename
        %   4. Index type (0 for uint32 and 1 for uint64)
        %   5. A logical that indicates if the entries are contiguous in the
        %   source file
        %   6. A logical that indicates if the source file can be indexed by
        %   keys
        %   7. A scalar that indicates which column contains the key in
        %   table formats
        %   8. keyToken
        %   9. A prefix that is used to denote lines that are comments in
        %   the source file
        %   10. A prefix that is used to denote lines that are headers in the
        %   source file
        %   11. The type of delimiter expected in tables, can be ' ' or \t
        %   12. The type of delimiter expected in flat files
        %
        % OUTPUT ARGUMENTS:
        %
        %   1. The number of entries found in the source file
        %   2. The number of text lines indexed from the source file
        
    end  %methods (Static, Access = private)
    
end % classdef

function lv = strdiff(strs)
    % Finds if contiguous strings are not equal in a cell array of strings
    n = numel(strs);
    lv = true(n,1);
    for j = 1:10000:n-1
        u = j:min(j+9999,n-1);
        lv(u) = ~strcmp(strs(u),strs(u+1));
    end
end

function [mapKeys,memMapIndex,interpreter,columnWithKeys,keyToken,...
    contiguousEntries,commentPrefix,headerPrefix,tableDelimiter,...
    entryDelimiter, verbose] ...
    = parse_inputs_app_specific_formats(fileFormat,varargin)
% Parse input PV pairs for application specific formats.

%=== Check for the right number of inputs
if rem(numel(varargin),2) == 1
    error(message('bioinfo:BioIndexedFile:BioIndexedFile:IncorrectNumberOfArguments'))
end

%=== Defaults
mapKeys           = true;
memMapIndex       = true;
interpreter       = [];
columnWithKeys    = 1;
keyToken          = '';
contiguousEntries = false;
commentPrefix     = '';
headerPrefix      = '';
tableDelimiter    = '\t';
entryDelimiter    = '//';
verbose           = true;

%=== fileFormat specific defaults
switch fileFormat
    case 'sam'
        interpreter  = @samread;
        headerPrefix = '@';
        contiguousEntries = true;
    case 'fastq'
        interpreter  = @(x) fastqread(x,'TRIMHEADERS',true);
        contiguousEntries = true;
    case 'fasta'
        interpreter  = @(x) fastaread(x,'TRIMHEADERS',true);
        contiguousEntries = true;
end

%=== Allowed inputs
okargs = {'mapkeys','memmapindex','interpreter','keycolumn',...
    'keytoken','contiguousentries','commentprefix','headerprefix',...
    'tabledelimiter','entrydelimiter','verbose','memorymappedindex',...
    'indexedbykeys'};

% R2011a: 'memmapindex' input argument changed to have the same name as the
% property 'memorymappedindex', however it will be silently still valid.
% 
% R2011a: 'mapkeys' input argument changed to have the same name as the
% property 'indexedbykeys', however it will be silently still valid.

for j = 1:2:numel(varargin)-1
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case {1,13}   % mapKeys or 'indexedbykeys'
            mapKeys = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case {2,12}   % memMapIndex or 'memorymappedindex'
            memMapIndex = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 3   % interpreter
            if isa(pval,'function_handle')
                interpreter = pval;
            elseif isempty(pval)
                interpreter = [];
            else
                error(message('bioinfo:BioIndexedFile:BioIndexedFile:InvalidInterpreterInput'))
            end
        case 11 % verbose
            verbose = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        otherwise
            warning(message('bioinfo:BioIndexedFile:BioIndexedFile:IgnoredInputs', fileFormat, upper( okargs{ k } )))
    end
end

end

function [mapKeys,memMapIndex,interpreter,columnWithKeys,keyToken,...
    contiguousEntries,commentPrefix,headerPrefix,tableDelimiter,...
    entryDelimiter, verbose] ...
    = parse_inputs_general_purpose_formats(varargin)
% Parse input PV pairs for general purpose formats.

%=== Check for the right number of inputs
if rem(nargin,2) == 1
    error(message('bioinfo:BioIndexedFile:BioIndexedFile:IncorrectNumberOfArguments'))
end

%=== Defaults
mapKeys           = true;
memMapIndex       = true;
interpreter       = [];
columnWithKeys    = 1;
keyToken          = '';
contiguousEntries = false;
commentPrefix     = '';
headerPrefix      = '';
tableDelimiter    = '\t';
entryDelimiter    = '//';
verbose           = true;

%=== Allowed inputs
okargs = {'mapkeys','memmapindex','interpreter','keycolumn',...
    'keytoken','contiguousentries','commentprefix','headerprefix',...
    'tabledelimiter','entrydelimiter','verbose','memorymappedindex',...
    'indexedbykeys'};

% R2011a: 'memmapindex' input argument changed to have the same name as the
% property 'memorymappedindex', however it will be silently still valid.
% 
% R2011a: 'mapkeys' input argument changed to have the same name as the
% property 'indexedbykeys', however it will be silently still valid.

for j = 1:2:nargin-1
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case {1,13}   % mapKeys or 'indexedbykeys'
            mapKeys = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case {2,12}   % memMapIndex or 'memorymappedindex'
            memMapIndex = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 3   % interpreter
            if isa(pval,'function_handle')
                interpreter = pval;
            elseif isempty(pval)
                interpreter = [];
            else
                error(message('bioinfo:BioIndexedFile:BioIndexedFile:InvalidInterpreterInput'))

            end
        case 4   %'keycolumn'
            if isnumeric(pval) && isscalar(pval) && (pval>=0) && (rem(pval,1)==0)
                columnWithKeys = pval;
            else
                error(message('bioinfo:BioIndexedFile:BioIndexedFile:InvalidColumnWithKeys'))

            end
        case 5   % keyToken
            if ischar(pval)
                keyToken = pval;
            else
                error(message('bioinfo:BioIndexedFile:BioIndexedFile:InvalidKeyToken'))

            end
        case 6   % contiguousEntries
            contiguousEntries = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 7   % commentPrefix
            if ischar(pval)
                commentPrefix = pval;
            else
                error(message('bioinfo:BioIndexedFile:BioIndexedFile:InvalidCommentPrefix'))
            end
        case 8   % headerPrefix
            if ischar(pval)
                headerPrefix = pval;
            else
                error(message('bioinfo:BioIndexedFile:BioIndexedFile:InvalidHeaderPrefix'))
            end
        case 9   % tableDelimiter
            if ischar(pval) && (isequal(pval,'\t') || isequal(pval,' ') || isequal(pval,','))
                tableDelimiter = pval;
            else
                error(message('bioinfo:BioIndexedFile:BioIndexedFile:InvalidTableDelimiter'))
            end
        case 10  % entryDelimiter
            if ischar(pval)
                entryDelimiter = pval;
            else
                error(message('bioinfo:BioIndexedFile:BioIndexedFile:InvalidEntryDelimiter'))
            end
        case 11 % verbose
            verbose = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end
end

end

