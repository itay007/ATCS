classdef BioRead < BioSeq
	%BIOREAD class representing a collection of sequences and their quality scores.
	%
	%  A BioRead object represents a collection of biological sequences and
	%  their quality scores.
	%
	%  BioRead properties:
	%  Quality  - quality scores for each position in a given sequence.
	%  Sequence - sequence data.
	%  Header   - header information for each sequence.
	%  NSeqs    - number of sequences represented in the object.
	%  Name     - name of the object.
	%
    %  BioRead specialized property access:
	%  get            - get multiple properties with a single statement.
	%  getHeader      - subscripted retrieval of the 'Header' property.
	%  getQuality     - subscripted retrieval of the 'Quality' property.
	%  getSequence    - subscripted retrieval of the 'Sequence' property.
	%  getSubsequence - subscripted retrieval of partial values of the 'Sequence' property.
	%  set            - set multiple properties with a single statement.
	%  setHeader      - subscripted assignment of the 'Header' property.
	%  setQuality     - subscripted assignment of the 'Quality' property.
	%  setSequence    - subscripted assignment of the 'Sequence' property.
	%  setSubsequence - subscripted assignment of partial values of the 'Sequence' property.
    %        
	%  BioRead methods:
	%  BioRead        - create a BioRead object.
	%  combine        - combine two objects.
	%  getSubset      - retrieve a subset of elements.
    %  plotSummary    - quality summary plots.
	%  setSubset      - set the values of a subset of elements.
	%
	%  See also BIOMAP.
	
	%  Copyright 2009-2012 The MathWorks, Inc.
	
	properties (Dependent)	
		%QUALITY quality scores of each position of the sequences represented in the object.
		%     The 'Quality' property is a cell array of strings containing
		%     the quality scores in ASCII format for each position of each
		%     sequence in the object.
		%
		%     See also BIOREAD.
		Quality
    end
    
    properties (Access = private)

    end
	
	%======================================================================
	% PUBLIC METHODS
	%======================================================================
	
	methods
		
		%==================================================================
		% Constructor
		%==================================================================
		
		
		function obj = BioRead(in, varargin)
			%BIOREAD create a BioRead object.
			%
            %  B = BIOREAD(FILENAME) creates a BioRead object B from the content of a 
            %  FASTQ/SAM formatted file FILENAME. Data is kept in the file and the 
            %  BioRead object accesses it using an auxiliary index file with the 
            %  extension .IDX. If the index file is not present BioRead will create one
            %  in the same directory where the source data is.            
            %
			%  B = BIOREAD(STRUCT) creates a BioRead object B from the content of a 
            %  MATLAB structure STRUCT, such as the structure returned by the FASTQREAD 
            %  and SAMREAD functions.
            %            
			%  B = BIOREAD(SEQS) creates a BioRead object B from a cell array of 
			%  strings SEQS representing the values of the BioRead property 'Sequence'.
			%
			%  B = BIOREAD(SEQS, QUALS) creates a BioRead object B and sets the BioRead 
            %  property 'Quality' to QUALS. QUALS is a cell array of strings.
			%  
			%  B = BIOREAD(SEQS, QUALS, HEADS) creates a BioRead object B and sets the 
            %  BioRead property 'Header' to HEADS. HEADS is a cell array of strings.
            %
            %  B = BIOREAD(..., 'InMemory', true) forces data to be loaded into memory 
            %  when using an auxiliary index file. InMemory defaults to false. InMemory 
            %  is ignored when the first input is not a file name.  
            %
            %  B = BIOREAD(..., 'IndexDir', INDEXDIR) sets the path to the directory 
            %  where the index files can be found or created. 
            %
            %  B = BIOREAD(..., PROPERTY_NAME , VALUE) specifies optional parameter 
            %  name/value pairs to further modify the object properties when allowed. 
            %  PROPERTY_NAME is a string.
			%
			%  Examples:
			%
			%  % Create a BioRead object from a FASTQ file.
			%  obj = BioRead('SRR005164_1_50.fastq')
			%
			%  % Create a BioRead object from variables in the workspace.
			%  seqs = {randseq(10); randseq(15); randseq(20)};
			%  heads = {'H1'; 'H2'; 'H3'};
			%  quals = {repmat('!', 1, 10); repmat('%', 1, 15); repmat('&', 1,20)};
			%  obj = BioRead(seqs, quals, heads)
			%
			%  % Create a BioRead object from a valid structure.
			%  str = struct('Header', heads, 'Sequence', seqs, 'Quality', quals);
			%  obj = BioRead(str)
			%
			%  See also BIOMAP/BIOMAP, BIOREAD.
			
            if (nargin == 0)
                in = [];
            end
			
            if iscellstr(in)
                if numel(varargin)>0 && iscellstr(varargin{1})
                    if numel(varargin)>1  && iscellstr(varargin{2})
                        %  BIOREAD(SEQS,QUALS,HEADS,...)
                        varargin = {'Sequence', in, 'Quality', varargin{1}, 'Header', varargin{2:end}};
                        in = [];
                    else
                        %  BIOREAD(SEQS,QUALS,...)
                        varargin = {'Sequence', in, 'Quality', varargin{1:end}};
                        in = [];
                    end
                else
                   %  BIOREAD(SEQS,...)
                   varargin = {'Sequence', in, varargin{1:end}};
                   in = [];
                end
            end

            %=== Check if the (new) first input argument may be a PVP
            %    (building from only PV pairs with properties only is
            %    undocummented but required from some object methods)
            if  ~isempty(in) && ~isstruct(in) &&   ~((ischar(in) && (isrow(in) || isempty(in)))&&any(in=='.'))
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
                        inMemory = bioinfoprivate.opttf(varargin{i+1}, 'InMemory', ['BioRead:' mfilename]); 
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
                             error(message('bioinfo:BioRead:BioRead:InvalidIndexDir'))
                        end
                        if ~isdir(indexDir)
                            error(message('bioinfo:BioRead:BioRead:UnknownIndexDir', indexDir))
                        end
                    end
                end
            end   
            
            %=== Remove PVPs already processed from varargin, all remaining
            %    PVPs shall be only settable object properties
            if numel(removeVarargin)
                varargin(removeVarargin) = [];
            end       
            
            %=== Build BioRead from FileName
            if ischar(in) && (isrow(in) || isempty(in))
                if ~exist(in,'file')
                    error(message('bioinfo:BioRead:BioRead:FileNotFound', in))
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
                           {'Sequence','Header'   ,'Quality' },...
                           {'string'  ,'string'   ,'string'  },...
                           [false     ,false      ,false     ],...% searchable fields
                           [true      ,true       ,true      ],...% efficient fields
                           {'Sequence','QueryName','Quality' });  % field binding
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
                            error(message('bioinfo:BioRead:BioRead:InvalidSAMFormat', in.InputFile))
                        end
                        
                        obj.Index = bioinfoprivate.BIFSequenceDataAdapter( ...
                            in,...  % data container
                            {'Sequence','Header','Quality'},...
                            {'string','string','string'},...
                            [false,false,false],...  % searchable fields
                            [false,false,false],...    % efficient fields
                            {'Sequence','QueryName','Quality'},...  % field binding
                            [10,1,11]); % columnPosition
                    end
                    
                elseif strncmpi('.fastq',fileExt,6)
                    if inMemory
                        in = fastqread(in, 'trimheaders', true);
                        obj.Index = bioinfoprivate.structSequenceDataAdapter( ...
                           in,...  % data container
                           {'Sequence','Header','Quality' },...
                           {'string'  ,'string','string'  },...
                           [false     ,false   ,false     ],...% searchable fields
                           [true      ,true    ,true      ],...% efficient fields
                           {'Sequence','Header','Quality' });  % field binding
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
                            error(message('bioinfo:BioRead:BioRead:InvalidFASTQFormat', in.InputFile))
                        end
                        
                        obj.Index = bioinfoprivate.BIFSequenceDataAdapter( ...
                            in,...  % data container
                            {'Sequence','Header','Quality'},...
                            {'string','string','string'},...
                            [false,false,false],...  % searchable fields
                            [false,false,false],...    % efficient fields
                            {'Sequence','Header','Quality'},...  % field binding
                            [0,0,0]); % columnPosition
                    end
                    
                else
                    error(message('bioinfo:BioRead:BioRead:InvalidFile'))
                end
                
            %=== Build BioRead from MATLAB structure (such as returned by SAMREAD)
            elseif isstruct(in)
                if ~iscellstr(in(1).Sequence) % structure is an array of structures such as returned by file readers
                    if sum(isfield(in,{'Sequence', 'Quality', 'Header'})) == 3 % structure from fastqread
                        fieldBinding = {'Sequence','Header','Quality'};
                    elseif sum(isfield(in,{'Sequence', 'Quality', 'QueryName'})) == 3 % structure from bamread/samread
                        fieldBinding = {'Sequence','QueryName','Quality'};
                    elseif sum(isfield(in,{'Sequence', 'Header'})) == 2 % structure from fastqread without quality
                        fieldBinding = {'Sequence','Header','Quality'};
                    elseif sum(isfield(in,{'Sequence', 'QueryName'})) == 2 % structure from bamread/samread without quality
                        fieldBinding = {'Sequence','QueryName','Quality'};
                    elseif sum(isfield(in,{'Sequence', 'Quality'})) == 2 % structure without header
                         fieldBinding = {'Sequence','Header','Quality'};
                    elseif isfield(in,'Sequence') % structure without header and quality
                         fieldBinding = {'Sequence','Header','Quality'};                        
                    else
                        error(message('bioinfo:BioRead:BioRead:InvalidInputStructure'))
                    end
                elseif sum(isfield(in, properties(obj))) == numel(properties(obj))  % structure from get(obj)
                    varargin = {'Sequence', in.Sequence, 'Quality', in.Quality,...
                        'Header', in.Header, varargin{1:end}};
                    in = [];
                    fieldBinding = {'Sequence','Header','Quality'};
                else
                    error(message('bioinfo:BioRead:BioRead:InvalidInputStructure'))
                end
                
                obj.Index = bioinfoprivate.structSequenceDataAdapter( ...
                           in,...  % data container
                           {'Sequence','Header','Quality' },...
                           {'string'  ,'string','string'  },...
                           [false     ,false   ,false     ],...% searchable fields
                           [true      ,true    ,true      ],...% efficient fields
                           fieldBinding);  % field binding
            elseif isempty(in)
                obj.Index = bioinfoprivate.structSequenceDataAdapter( ...
                    in,...  % data container
                    {'Sequence','Header','Quality' },...
                    {'string'  ,'string','string'  },...
                    [false     ,false   ,false     ],...% searchable fields
                    [true      ,true    ,true      ],...% efficient fields
                    {'Sequence','Header','Quality' });  % field binding
            else
                error(message('bioinfo:BioRead:BioRead:InvalidInput'))
            end
            
            %=== Populate the object with PV pairs
            if rem(numel(varargin),2) ~= 0
                error(message('bioinfo:BioRead:BioRead:InvalidInputNumber'))
            end                
            for k = 1:2:numel(varargin)
                %=== make sure Sequence PVP is handled first, so size of BioRead is predetermined.
                seqPar = find(~cellfun(@isempty, strfind(varargin(1:2:end), 'Sequence')));
                if ~isempty(seqPar)
                    varargin = {varargin{seqPar(1)*2-1} varargin{seqPar(1)*2} varargin{:}}; %#ok<CCAT> % add Sequence PVP at beginning
                    varargin([seqPar(1)*2+1 seqPar(1)*2+2]) = ''; % remove Sequence PVP from original position
                end
                obj.(varargin{k})= varargin{k+1};
            end
            
        end % end constructor            

        %==================================================================
		% get.Quality
		%==================================================================
        function h = get.Quality(obj)
            h = getField(obj.Index,'Quality');
        end
		%==================================================================
		% set.Quality
		%==================================================================
		function obj = set.Quality(obj,h)
            obj.Index = setField(obj.Index,'Quality',h);
        end
		
	end % methods
	
end %BioRead
