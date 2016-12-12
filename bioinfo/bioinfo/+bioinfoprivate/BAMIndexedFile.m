classdef BAMIndexedFile
	%BAMINDEXEDFILE class allows random read access to BAM files.
    %
	%  A BAMINDEXEDFILE object allows efficient random read access to a BAM
	%  formatted file. BAMINDEXEDFILE requires two auxiliary index files,
	%  BAI and LINEARINDEX, to access the data, if these files do not exist
	%  BAMINDEXEDFILE prepares the BAI and LINEARINDEX files and saves them
	%  for further use.
    %
    %  BAMIndexedFile properties:
    %    FileFormat        - File format of the source file.
    %    InputFile         - Path and filename of the source file.
    %    IndexFile         - Path and filename of the BAI file.
    %    LinearIndexFile   - Path and filename of the LINEARINDEX file.
    %    NumEntries        - Number of entries indexed by the object.
    %    Reference         - Name of the reference sequence.
    %    FeatureName       - Names of accessible features.
    %    FeatureInMemory   - Indicates which features are loaded in memory.
    %    FeatureInIndex    - Indicates which features are saved in the LINEARINDEX file.   
    %    LinearIndexInMemory - Indicates which linear indices are loaded in memory.
    %
    %  FeatureInMemory, FeatureInIndex and LinearIndexInMemory are the only
    %  properties that can be changed once the object has been constructed.
    %  All other properties are read-only. Use FeatureInMemory to load and
    %  unload features into/from memory. Use FeatureInIndex to save
    %  properties in the linear index file for efficient access. Use
    %  LinearIndexInMemory to load and unload linear indices for every
    %  reference in the BAM file. These linear indices must be of the same
    %  size as entries in the source file in order to allow efficient
    %  non-contiguous access, so, contrary to the features loaded in
    %  memory, the size of the linear indices does not decrease when the 
    %  object is subsetted. If you will only work with one of the
    %  references you can unload the linear indices for other reference,
    %  thus saving memory. 
    %        
    %  BAMIndexedFile methods:
    %    BAMIndexedFile    - Create a BAMIndexedFile object.
    %    getFeature        - Retrieve features using numeric index.
    %    getSubset         - Returns an object that indexes a subset of the 
    %                        entries.
    
    %   Copyright 2011-2012 The MathWorks, Inc.
    
    properties (SetAccess = private, GetAccess = public)
        %FILEFORMAT file format of the source file.
		%     The 'FileFormat' property indicates the format of the source
		%     file. 'FileFormat' is always 'bam'.
        %
        %  See also BAMIndexedFile, BAMIndexedFile.BAMIndexedFile.
        FileFormat
        %INPUTFILE path and filename of the source file.
        %    The 'InputFile' property contains a char array with the path
        %    and filename of the source file.
        %
        %  See also BAMIndexedFile, BAMIndexedFile.BAMIndexedFile.
        InputFile
        %INDEXFILE path and filename of the BAI file.
        %    The 'IndexFile' property contains a char array with the path
        %    and filename of the auxiliary index file.
        %
        %  See also BAMIndexedFile, BAMIndexedFile.BAMIndexedFile.        
        IndexFile
        %LINEARINDEXFILE path and filename of the LINEARINDEX file.
        %    The 'LinearIndexFile' property contains a char array with the
        %    path and filename of the auxiliary linear index file.
        %
        %  See also BAMIndexedFile, BAMIndexedFile.BAMIndexedFile.        
        LinearIndexFile
        %NUMENTRIES number of entries indexed by the object.
        %    The 'NumEntries' property contains a integer indicating the
        %    number of entries that are indexed by the object.
        %
        %  See also BAMIndexedFile, BAMIndexedFile.BAMIndexedFile.        
        NumEntries
        %REFERENCE name of the reference sequence.
        %    The 'Reference' property contains a string with the name of
        %    the reference sequence.
        %
        %  See also BAMIndexedFile, BAMIndexedFile.BAMIndexedFile.          
        Reference
        %FEATURENAME names of accessible features.
        %    The 'FeatureName' property contains a cell string with the
        %    names of the features that can be accessed by this object.
        %
        %  See also BAMIndexedFile, BAMIndexedFile.BAMIndexedFile.         
        FeatureName
    end
    
    properties (Dependent, SetAccess = public, GetAccess = public)
        %FEATUREINMEMORY indicates which features are loaded in memory.
        %    The 'FeatureInMemory' property contains a cell array with
        %    true/false values indicating which features are available in
        %    memory. The number of rows is equivalent to the number of
        %    references and the number of columns is equivalent to the
        %    number of features.
        %
        %  See also BAMIndexedFile, BAMIndexedFile.BAMIndexedFile.         
        FeatureInMemory
        %FEATUREININDEX indicates which features are available in the linear index file.
        %    The 'FeatureInIndex' property contains a cell with true/false
        %    values indicating which features are available in the linear
        %    index file. The number of rows is equivalent to the number of
        %    references and the number of columns is equivalent to the
        %    number of features.
        %
        %  See also BAMIndexedFile, BAMIndexedFile.BAMIndexedFile.         
        FeatureInIndex   
        %LINEARINDEXINMEMORY indicates which linear indices are loaded in memory.
        %    The 'LinearIndexInMemory' property contains a cell with
        %    true/false values indicating which linear indices are
        %    available in memory. The number of elements is equivalent to
        %    the number of references.
        LinearIndexInMemory
    end    

    properties (SetAccess = private, GetAccess = private, Hidden)
        IndexedByKeys
    end
   
    properties  (SetAccess = private, GetAccess = public, Hidden)
        SubsetIndex
    end
    properties  (SetAccess = private, GetAccess = private)
        ReferenceIndex
        LinearIndex
        ReverseSubsetIndex
        PrivateLinearIndexInMemory
        PrivateFeatureInMemory
        PrivateFeatureInIndex
        FeatureType
        FeatureData 
        NumEntriesPerRefInSubset
        NumEntriesPerRefInSource
    end
            
    methods
		
		%==================================================================
		% Constructor
		%==================================================================
        
        function obj = BAMIndexedFile(inputFile,varargin)
        %BAMINDEXEDFILE Create a object to access BAM formatted files.
        %
        %  BAMINDEXEDFILE(SOURCEFILE) constructs a BAMIndexedFile object
        %  that provides efficient random access to a BAM formatted file
        %  SOURCEFILE. SOURCEFILE is a string containing a file name and
        %  may be prefixed by a relative or absolute path. 
        %
        %  The BAMINDEXEDFILE object uses two auxiliary files to store 
        %  information that allows the efficient direct access to
        %  SOURCEFILE, the BAI and the LINEARINDEX files. The name of the
        %  auxiliary files by default are derived from SOURCEFILE by 
        %  appending the extensions '.bai' and '.linearindex'. The path to
        %  the auxiliary files defaults to the same directory as
        %  SOURCEFILE. When the auxiliary files do not exist, the
        %  BAMINDEXEDFILE constructor creates them by analyzing SOURCEFILE,
        %  otherwise, the object is built from the information stored in
        %  the auxiliary files. 
        %
        %  BAMINDEXEDFILE(SOURCEFILE,INDEXDIR) sets a different directory
        %  to write or search for the auxiliary files (BAI and LINEARINDEX).
        %
        %  BAMINDEXEDFILE(SOURCEFILE,INDEXFILE) sets the base name of the
        %  auxiliary files, the file extensions '.bai' and '.linearindex'
        %  are appended to the base name. INDEXFILE may be prefixed by a
        %  relative or absolute path.
        %
        %  BAMINDEXEDFILE(...,'REFERENCE',NAME) selects one of the
        %  references in the source file.
        %
        %  BAMINDEXEDFILE(...,'VERBOSE',TF) displays information about the
        %  object construction. Defaults to true.
        %
        %  See also BAMIndexedFile, BAMIndexedFile.getFeature,
        %           BAMIndexedFile.getSubset.
            
            %==== Validate input file name
            if isempty(fileparts(inputFile))
                inputFile = [pwd filesep inputFile];
            end
            if ~exist(inputFile,'file')
                [p,n,e] = fileparts(inputFile);
                error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:SourceFileNotFound', [ n, e ], p))
            end
            
            %==== Finds out the signature by looking at numel(varargin)
            if rem(numel(varargin),2)==0
                baseName = inputFile;
            elseif ischar(varargin{1})
                baseName = varargin{1};
                varargin = varargin(2:end);
            end
                       
            %==== Parse PV pairs
            [reference,verbose] = parse_inputs(varargin{:});
            
            %==== Validate index file name
            if baseName(end)=='.'
                % the user meant a relative path ending with '.'
                baseName = [baseName filesep];
            end
            [p,n,e] = fileparts(baseName);
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
                error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:InvalidDirectory', p))
            end
            if isempty(n)
                % the user only passed a directory get base name for the
                % auxiliary files using the name and extension of the
                % source file (inputFile)
                [~,n,e] = fileparts(inputFile); %default name derived from source
            end
            
            if strcmpi(e,'.bai') || strcmpi(e,'.linearindex')
                % make sure that the baseName does not already have one of
                % the extensions of the auxiliary files
                error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:InvalidIndexFile'))
            end
            
            indexFile = [p filesep n e '.bai'];
            balFile = [p filesep n e '.linearindex'];
            
            indexFileExists = exist(indexFile,'file');
            balFileExists = exist(balFile,'file');
            
            if balFileExists && ~indexFileExists
                 error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:StaleLinearIndexFile'))
            end
            
            if verbose
                d = dir(inputFile);
                fprintf('Source File: %s\n',d.name)
                fprintf('   Path: %s\n',fileparts(inputFile))
                fprintf('   Size: %d bytes\n',d.bytes)
                fprintf('   Date: %s\n',datestr(d.datenum))
            end
            
            %==== Create BAI file when necessary:
            if ~indexFileExists
                % Index not found, create a new one:
                if verbose
                    fprintf('Creating BAI file ...\n')
                end
                
                % Update status in NGS Browser
                if usejava('jvm')                    
                    com.mathworks.toolbox.bioinfo.genome.GenomeDesktop.updateStatusMessage('Creating BAI file ...');
                    clearIdxBrowserStatusObj = onCleanup(@()com.mathworks.toolbox.bioinfo.genome.GenomeDesktop.updateStatusMessage(''));
                end

                if ~bioinfoprivate.bamaccessmex('baminfo',inputFile)
                    % Source file is not ordered, if the file has already
                    % been ordered it should have the following name:
                    orderedBAMfile = regexprep(indexFile,'.bam.bai','.ordered.bam');
                    if ~exist(orderedBAMfile,'file')
                        if exist(['orderedBAMfile' '.bai'],'file') || ...
                           exist(['orderedBAMfile' '.linearindex'],'file')  
                           error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:StaleBAIorLinearIndexFile'))
                        end
                        try
                            bioinfoprivate.bamaccessmex('bamsort', inputFile,regexprep(orderedBAMfile,'.bam',''));
                        catch ME2
                            % When there is an error in the orderee it is likely that
                            % the file left is invalid, then it is deleted.
                            if exist(indexFile,'file')
                                delete(orderedBAMfile)
                            end
                            bioinfoprivate.bioclsrethrow(mfilename,mfilename,ME2)                            
                        end
                    end
                    obj = bioinfoprivate.BAMIndexedFile(orderedBAMfile,varargin{:});
                    return
                    % error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:UnorderedBAMFile'))
                end
                try
                    bioinfoprivate.bamaccessmex('indexbamfile',inputFile,indexFile);
                catch ME
                    % When there is an error in the indexer it is likely that
                    % the file left is invalid, then it is deleted.
                    if exist(indexFile,'file')
                        delete(indexFile)
                    end
                    bioinfoprivate.bioclsrethrow(mfilename,mfilename,ME)
                end
                if verbose
                    fprintf('%s file created.\n',indexFile)
                end
                
                % Clear NGS browser status message
                clear('clearIdxBrowserStatusObj');
                
            elseif verbose
                fprintf('Auxiliary file BAI exists.\n')
            end
            
            %==== Create BAL file when necessary:
            if ~balFileExists
                % Index not found, create a new one:
                if verbose
                    fprintf('Creating LINEARINDEX file ...\n')
                end
                
                % Update status in NGS Browser
                if usejava('jvm')                    
                    com.mathworks.toolbox.bioinfo.genome.GenomeDesktop.updateStatusMessage('Creating LINEARINDEX file ...');
                    clearBrowserStatusObj = onCleanup(@()com.mathworks.toolbox.bioinfo.genome.GenomeDesktop.updateStatusMessage(''));
                end
                
                BAL_version = '1.0';
                try
                    save(balFile,'BAL_version','-v7.3')
                catch ME
                    if exist(balFile,'file')
                        delete(balFile)
                    end
                    bioinfoprivate.bioclsrethrow(mfilename,mfilename,ME)
                end
                if verbose
                    fprintf('%s file created.\n',balFile)
                end
                
            elseif verbose
                fprintf('Auxiliary file LINEARINDEX exists.\n')
            end
                
            if verbose
                d = dir(indexFile);
                fprintf('Index File: %s\n',d.name)
                fprintf('   Path: %s\n',fileparts(indexFile))
                fprintf('   Size: %d bytes\n',d.bytes)
                fprintf('   Date: %s\n',datestr(d.datenum))
            end

            %=== Validate the reference name
            nrefsinfile = bioinfoprivate.bamaccessmex('getnumberofreferences',inputFile);
            if isempty(reference)
                % In 12b when reference is empty (undefined) the object
                % will index all available references in the file.
                referenceIndex = int32(0:nrefsinfile-1);
                nrefs = nrefsinfile;
            elseif ischar(reference) || iscellstr(reference)
                if ischar(reference)
                    reference = {reference};
                end
                nrefs = numel(reference);
                referenceIndex = -ones(1,nrefs,'int32');
                for i = 1:nrefs
                    referenceIndex(i) = bioinfoprivate.bamaccessmex('validatereferencestring',inputFile,reference{i});
                end
                if any(referenceIndex<0)
                    hole = regexprep(sprintf('%s, ',reference{referenceIndex<0}),',\s$','');
                    error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:ReferenceNotInBAMFile', hole))
                end
                referenceIndex = sort(referenceIndex);
            elseif isnumeric(reference) && isvector(reference)    
                referenceIndex = sort(int32(reference(:)));
                nrefs = numel(reference);
                if any(referenceIndex<0) || any(referenceIndex>=nrefsinfile)
                    error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:ReferenceIndexNotInBAMFile', nrefsinfile))
                end
            else
                error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:InvalidReference'))
            end
            
            reference = cell(nrefs,1);
            for i = 1:nrefs
                reference{i} = bioinfoprivate.bamaccessmex('getreferencestring',inputFile,referenceIndex(i));
            end
            
            if verbose
                fprintf('Reference Index: %s\n',num2str(referenceIndex(:)'))
            end
           
            %==== Object Construction
            obj.FileFormat        = 'bam';
            obj.InputFile         = inputFile;
            obj.IndexFile         = indexFile;
            obj.LinearIndexFile   = balFile;
            obj.ReferenceIndex    = referenceIndex;
            obj.SubsetIndex       = [];
            obj.ReverseSubsetIndex = [];
            obj.Reference         = reference;
            obj.IndexedByKeys     = false;
            obj.FeatureName       = {'Flag' 'Position' 'MappingQuality' 'MatePosition' 'MateIndex' 'QueryLength' 'EndPosition' 'Signature' 'Sequence' 'Header' 'Quality'};
            obj.FeatureType       = {'uint16' 'uint32' 'uint8'          'uint32'       'uint32'    'uint16'      'uint32'      'String'    'String'   'String' 'String'};
            nfeats = numel(obj.FeatureName); 
            obj.FeatureData       = repmat({[]},nrefs,nfeats);
            obj.LinearIndex    = repmat({[]},nrefs,1);
            obj.PrivateLinearIndexInMemory = repmat({false},nrefs,1);
            obj.PrivateFeatureInMemory = repmat({false},nrefs,nfeats);
            obj.PrivateFeatureInIndex  = repmat({false},nrefs,nfeats);
            obj.NumEntriesPerRefInSubset = zeros(nrefs,1);
            obj.NumEntriesPerRefInSource = zeros(nrefs,1);
            
            %==== Update obj.PrivateFeatureInIndex according with
            %     information in BAL file
            wo = properties(matfile(obj.LinearIndexFile));
            for j = 1:nrefs
                for i = 1:nfeats
                  strtmp = sprintf('%s_%d',obj.FeatureName{i},obj.ReferenceIndex(j));
                  obj.PrivateFeatureInIndex{j,i} = any(strcmp(strtmp,wo));
                end
            end
            
            %==== Make sure Position, EndPosition and Flag are in the
            %     LINEARINDEX file for each one of the references contained
            %     in this object
            
            % At construction time obj.LinearIndex{j} is empty, only
            % then we use the bamaccessmex('getstartsinref',...) instead of
            % the getFeature(obj,...) method.
            % bamaccessmex('getstartsinref',...) will load the start and
            % end positions and the flags for the whole reference in a BAM
            % file without the need of linear indexing or subsetting
            % support.            
            tfs = repmat({true true false false false false true false false false false},nrefs,1);
            cfs = obj.FeatureInIndex;
            ftl = (cell2mat(cfs)==false) & (cell2mat(tfs)==true);
            mf = matfile(obj.LinearIndexFile,'writable',true);
            for j = 1:numel(obj.ReferenceIndex)
                ftlr = ftl(j,:); % row for this reference
                if any(ftlr([1 2 7]))
                    s1 = sprintf('Position_%d',obj.ReferenceIndex(j));
                    s2 = sprintf('EndPosition_%d',obj.ReferenceIndex(j));
                    s3 = sprintf('Flag_%d',obj.ReferenceIndex(j));
                    [mf.(s1),mf.(s2),mf.(s3)] = bioinfoprivate.bamaccessmex('getstartsinref',...
                        obj.InputFile,obj.IndexFile(1:end-4),obj.ReferenceIndex(j));
                end
                tfs(j,[1 2 7]) = {true}; % indicate that all tree are indexed
            end
            obj.PrivateFeatureInIndex = tfs;
            
           
            %==== Load LinearIndex (should be available in the Index
            %     File by now). The LinearIndex are never subseted 
            %     and will be later used for efficient linear indexing into
            %     the BAM file.
            %obj.LinearIndexInMemory(:) = {true};
            mf = matfile(obj.LinearIndexFile);                
            for j = 1:numel(obj.ReferenceIndex)
                % In 12b we use lazy load for LinearIndex
                %obj.LinearIndex{j} = mf.(sprintf('Position_%d',obj.ReferenceIndex(j)))-1;
                n = numel(mf.(sprintf('Position_%d',obj.ReferenceIndex(j))));
                obj.NumEntriesPerRefInSubset(j) = n;
                obj.NumEntriesPerRefInSource(j) = n; 
            end
                
            obj.NumEntries = sum(obj.NumEntriesPerRefInSubset);
            
            % Clear NGS Browser status message
            clear('clearBrowserStatusObj');
            
            if verbose
                fprintf('Done. \n')
            end
            
            
        end % IndexedFile

        
        %==================================================================
        % set.FeatureInIndex
        %==================================================================        
        function obj = set.FeatureInIndex(obj, tfs)
            % Validate a cell array with true/false flags
            if ~iscell(tfs) 
                  error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:InvalidFeatureInIndex'))
            end
            if (numel(obj.ReferenceIndex) ~= size(tfs,1)) || (numel(obj.FeatureName) ~= size(tfs,2))
                error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:BadSizeFeatureInIndex'))
            end
            try
                for i = 1:numel(tfs)
                   tfs{i} = bioinfoprivate.opttf(tfs{i},'FeatureInIndex',mfilename); 
                end
            catch ME
                if strcmp(ME.identifier,'bioinfo:BAMIndexedFile:FeatureInIndexOptionNotLogical')
                    error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:InvalidElementFeatureInIndex'))
                else
                    rethrow(ME)
                end
            end
            
            % Check what is the current status of features saved in the index file:
            cfs = obj.FeatureInIndex;  
            
            % Features to remove from the index file:
            ftu = (cell2mat(cfs)==true) & (cell2mat(tfs)==false);
            if any(ftu(:))
                % Features already in the index file cannot be removed
                tfs{ftu} = true;
            end

            % Non numeric features
            if any(strcmp('String',obj.FeatureType(any(cell2mat(tfs),1))))
                 error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:NonNumericFeatureInIndex'))
            end
            
            % Features to add to the index file
            ftl = (cell2mat(cfs)==false) & (cell2mat(tfs)==true);
            
            if any(any(ftl(:,[3:6 8:end])))
                 error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:NoSupportFeatureInIndex'))
            end
            
            obj.PrivateFeatureInIndex = tfs;

        end

        %==================================================================
        % set.FeatureInMemory
        %==================================================================
        function obj = set.FeatureInMemory(obj, tfs)
            % Validate a cell array with true/false flags
            if ~iscell(tfs) 
                  error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:InvalidFeatureInMemory'))
            end
            if (numel(obj.ReferenceIndex) ~= size(tfs,1)) || (numel(obj.FeatureName) ~= size(tfs,2))
                error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:BadSizeFeatureInMemory'))
            end
            try
                for i = 1:numel(tfs)
                   tfs{i} = bioinfoprivate.opttf(tfs{i},'FeatureInMemory',mfilename); 
                end
            catch ME
                if strcmp(ME.identifier,'bioinfo:BAMIndexedFile:FeatureInMemoryOptionNotLogical')
                    error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:InvalidElementFeatureInMemory'))
                else
                    rethrow(ME)
                end
            end
           
            % Check what is the current status of features loaded in memory:
            cfs = obj.FeatureInMemory;
            
            % Features to unload:
            ftu = (cell2mat(cfs)==true) & (cell2mat(tfs)==false);
            if any(ftu(:))
              obj.FeatureData(ftu) = {[]};
            end
            
            % Features to load:
            ftl = (cell2mat(cfs)==false) & (cell2mat(tfs)==true);
            
            if any(strcmp('String',obj.FeatureType(any(ftl,1))))
                 error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:NonNumericFeatureInMemory'))
            end
            
            % Load features for each reference: getFeature method will get
            % them from the most efficient place and will subset if required.
            for j = 1:numel(obj.ReferenceIndex)
                ftlr = ftl(j,:);
                if any(ftlr)
                    idx = find(indexToRefIndex(obj,uint32(1):uint32(obj.NumEntries))==j);
                    if ~isempty(idx)
                        obj.FeatureData(j,ftlr) = getFeature(obj,obj.FeatureName(ftlr),idx);
                    end
                end
            end
            
            % Update object with features that are available in memory
            obj.PrivateFeatureInMemory = tfs;
        end
        
        %==================================================================
        % set.LinearIndexInMemory
        %==================================================================
        function obj = set.LinearIndexInMemory(obj, spr)    
            
            spc = obj.LinearIndexInMemory;
            % to unload:
            spu = (cell2mat(spc)==true) & (cell2mat(spr)==false);
            if any(spu)
              obj.LinearIndex{spu} = [];
            end
            % to load:
            spl = (cell2mat(spc)==false) & (cell2mat(spr)==true);
            
            if any(spl)
                mf = matfile(obj.LinearIndexFile);
                rr = find(spl);
                for i = 1:numel(rr)
                   obj.LinearIndex{rr(i)} = mf.(sprintf('Position_%d',obj.ReferenceIndex(rr(i))))-1;
                end
            end
            obj.PrivateLinearIndexInMemory = spr;
        end
        
        %==================================================================
        % get.FeatureInMemory
        %==================================================================
        function tf = get.FeatureInMemory(obj)
            tf = obj.PrivateFeatureInMemory;
        end
        %==================================================================
        % get.FeatureInIndex
        %==================================================================
        function tf = get.FeatureInIndex(obj)
            tf = obj.PrivateFeatureInIndex;
        end    
        %==================================================================
        % get.LinearIndexInMemory
        %==================================================================
        function tf = get.LinearIndexInMemory(obj)
            tf = obj.PrivateLinearIndexInMemory;
        end  
        
        %==================================================================
        % getSubsetReference
        %==================================================================
        function obj = getSubsetReference(obj,r)
            %GETSUBSETREFERENCE Returns an object that access a subset of the entries in the source file.
            %
            %  GETSUBSETREFERENCE(OBJ,R) Returns an object that accesses a
            %  subset of the entries in the source file. Entries for the
            %  subset are selected by R, a cell string with unique
            %  reference names.
            %
            %  See also BAMIndexedFile, BAMIndexedFile.BAMIndexedFile,
            %            BAMIndexedFile.getFeature. 
            
            if ~(iscellstr(r) && numel(unique(r))==numel(r))
                error(message('bioinfo:BAMIndexedFile:getSubsetReference:InvalidInput'))
            end
            
            h = false(1,numel(obj.ReferenceIndex));
            
            for i = 1:numel(r)
                h(strcmp(obj.Reference,r{i})) = true;
            end
                        
            UpdateSubsetIndexFlag = ~isempty(obj.SubsetIndex);
            
            if UpdateSubsetIndexFlag
                sr = indexToRefIndex(obj,1:obj.NumEntries);
                nes = [0;obj.NumEntriesPerRefInSource];
                nes([false h]) = 0;
                nes = uint32(cumsum(nes));
                p = int32(find(h));
            end

            obj.Reference = obj.Reference(h);
            obj.FeatureData = obj.FeatureData(h,:);
            obj.PrivateFeatureInIndex = obj.PrivateFeatureInIndex(h,:);
            obj.PrivateFeatureInMemory = obj.PrivateFeatureInMemory(h,:);
            obj.PrivateLinearIndexInMemory = obj.PrivateLinearIndexInMemory(h);
            obj.LinearIndex = obj.LinearIndex(h);
            obj.ReferenceIndex = obj.ReferenceIndex(h);
            obj.NumEntriesPerRefInSubset = obj.NumEntriesPerRefInSubset(h);
            obj.NumEntriesPerRefInSource = obj.NumEntriesPerRefInSource(h);
            
            if UpdateSubsetIndexFlag
                sl = false(numel(sr),1);
                if isempty(obj.ReverseSubsetIndex)
                    for j = 1:numel(p)
                        w = sr==p(j);
                        sl(w) = true;
                    end
                else
                    obj.ReverseSubsetIndex = zeros(numel(sr),1,'uint32');
                    for j = 1:numel(p)
                        w = sr==p(j);
                        sl(w) = true;
                        sw = uint32(sum(w));
                        obj.ReverseSubsetIndex(w) = 1:sw;
                    end
                    obj.ReverseSubsetIndex = obj.ReverseSubsetIndex(sl);
                end
                obj.SubsetIndex = obj.SubsetIndex(sl)-nes(sr(sl));
            end
            obj.NumEntries = sum(obj.NumEntriesPerRefInSubset);
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
            %  See also BAMIndexedFile, BAMIndexedFile.BAMIndexedFile,
            %            BAMIndexedFile.getFeature. 

            %  Access to BAM files is not indexed by key yet:
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

            
            if ischar(s)&&isvector(s)
                s = {s};
            end

            UpdateReverseIndexFlag = ~isempty(obj.ReverseSubsetIndex);
            if isnumeric(s)&&isscalar(s)
                if s<1 || s>obj.NumEntries || rem(s,1)
                    error(message('bioinfo:BAMIndexedFile:getSubset:invalidNumericInput'))
                end
                idx = uint32(s);
            elseif isnumeric(s)&&isvector(s)
                if any(s<1) || any(s>obj.NumEntries) || any(rem(s,1))
                    error(message('bioinfo:BAMIndexedFile:getSubset:invalidNumericInput'))
                end
                if issorted(s)
                    if any(diff(s)==0)
                        error(message('bioinfo:BAMIndexedFile:getSubset:invalidNumericInput'))
                    end
                elseif numel(unique(s)) ~= numel(s)
                    error(message('bioinfo:BAMIndexedFile:getSubset:invalidNumericInput'))
                else
                    UpdateReverseIndexFlag = true;
                end
                if isa(s,'uint32')
                    idx = s(:);
                else
                    idx = uint32(s(:));
                end
            elseif iscellstr(s)
                error(message('bioinfo:BAMIndexedFile:getSubset:invalidAlphanumericInput'))
            else
                error(message('bioinfo:BAMIndexedFile:getSubset:invalidInput'))
            end
            
            sr = indexToRefIndex(obj,idx);
            
           
            if isempty(obj.SubsetIndex)
                obj.SubsetIndex = idx;
            else
                obj.SubsetIndex = obj.SubsetIndex(idx);
            end
            
            % Change the global index to a reference index
            if isempty(obj.ReverseSubsetIndex)
                % When previous subset indices have been sorted or when no
                % subsetting has occurred:
                cs = uint32(cumsum([0;obj.NumEntriesPerRefInSubset]));
                idx = idx-cs(sr);
            else
                idx = obj.ReverseSubsetIndex(idx);
            end
            
            % Check if obj.FeatureData needs to be updated (Features loaded
            % in memory) so it will be updated to contain only the data of
            % the new subset: 
            cfs = cell2mat(obj.FeatureInMemory); 
            for j = 1:numel(obj.ReferenceIndex)
                cfsr = cfs(j,:);
                if any(cfsr)
                    for i = 1:numel(cfsr)
                        if cfsr(i)
                            data = obj.FeatureData{j,i};        
                            obj.FeatureData{j,i} = data(idx(sr==j));
                        end
                    end
                end
            end

            % Update ReverseSubsetIndex (when needed),
            % NumEntriesPerRefInSubset and NumEntries
            if UpdateReverseIndexFlag
                obj.ReverseSubsetIndex = zeros(numel(sr),1,'uint32');
                for j = 1:numel(obj.ReferenceIndex)
                    w = sr==j;
                    sw = uint32(sum(w));
                    obj.ReverseSubsetIndex(w) = 1:sw;
                    obj.NumEntriesPerRefInSubset(j) = sw;
                end
            else
                for j = 1:numel(obj.ReferenceIndex)
                    sw = uint32(sum(sr==j));
                    obj.NumEntriesPerRefInSubset(j) = sw;
                end
            end
            obj.NumEntries = sum(obj.NumEntriesPerRefInSubset);
            obj.LinearIndexInMemory(obj.NumEntriesPerRefInSubset==0) = {false};
            
        end

        %==================================================================
        % getFeature
        %==================================================================
        function [data,obj] = getFeature(obj,name,s)
            %GETFEATURE Query a feature from the index file
            %
            %  GETFEATURE(OBJ,NAME) Returns the data corresponding to the
            %  feature NAME. GETFEATURE determines the most efficient place
            %  from which the data can be pulled from, since the data may
            %  have already been loaded into memory or may exist stored
            %  linearly in an auxiliar file. Otherwise, the data will be
            %  pulled from the original source file. NAME may be also a
            %  cell string with multiple non-repeated feature names.
            %
            %  GETFEATURE(OBJ,NAME,I) Returns the data corresponding to the
            %  feature NAME and indexed by the numeric vector I. I contains
            %  positive integers less than or equal to the number of
            %  entries indexed by the object. There is a one-to-one
            %  relationship between the values in I and the output. Order
            %  and quantity are preserved, despite repeated values in I.
            %
            %  See also BAMIndexedFile, BAMIndexedFile.BAMIndexedFile,
            %           BAMIndexedFile.getSubset.  
            
            %  (R2012b) Access to BAM files is not indexed by key yet:
            %  GETFEATURE(OBJ,NAME,KEY) Returns the data corresponding to the
            %  feature NAME and identified by KEY. KEY is a string or a
            %  cell string with multiple keys. If keys are not unique in
            %  the source file, all entries that match a given key are
            %  included in the output at its respective position according
            %  to the query. Therefore, when keys are unique in the source
            %  file there is a one-to-one relationship between the strings
            %  in KEY and the output.
            
            % === Validate subsetting index and find out references to be
            %     accessed:
            if nargin<3 % All entries indexed by the object are requested
                r = uint32(1:numel(obj.ReferenceIndex));
                so = obj.NumEntries;
                si =  uint32(1):uint32(so);
                sr = indexToRefIndex(obj,si);
            elseif (ischar(s)&&isvector(s)) || iscellstr(s)
                error(message('bioinfo:BAMIndexedFile:getFeature:invalidAlphanumericInput'))
            elseif  isnumeric(s)&&isempty(s)
                si =  uint32([]);
                sr =  uint32([]);
                r = uint32([]);
                so = 0;
            elseif isnumeric(s)&&isvector(s)
                if any(s<1) || any(s>obj.NumEntries) || any(rem(s,1))
                    error(message('bioinfo:BAMIndexedFile:getFeature:nonNumeric'))
                else
                    if isa(s,'uint32')
                        si = s(:);
                    else
                        si = uint32(s(:));
                    end
                end
                % === Find indices of the references to be accessed
                sr = indexToRefIndex(obj,si);
                r = unique(sr);
                so = numel(si);
            else
                error(message('bioinfo:BAMIndexedFile:getFeature:invalidInput'))
            end
           
            if iscell(name) 
                dataIsCell = true;
            else
                dataIsCell = false;
            end
                
            % === Indices to features
            [h,hr,name] = findFeature(obj,name);

            % === Initialize output according with the indexing
            data = cell(size(name));
            
            for i = 1:numel(h)
                cl = obj.FeatureType{h(i)};
                if strcmp(cl,'String')
                    data{i} = cell(so,1); 
                else
                    data{i} = zeros(so,1,cl);
                end
            end
            
            for j = 1:numel(r)
                rr = r(j); % index to this reference
                w = sr == rr; % where to put the data
                if ~any(w)
                    % if there are non requested entries for this reference
                    continue;
                end
                 
                % Features that are available in memory:
                hm = find(cell2mat(obj.FeatureInMemory(rr,h)));
                % Features that are available in index file but not available in memory:
                hi = find(~cell2mat(obj.FeatureInMemory(rr,h)) & cell2mat(obj.FeatureInIndex(rr,h)));
                % Features that are neither in memory or the index file and
                % will be read from the source file
                hf = find(~cell2mat(obj.FeatureInMemory(rr,h)) & ~cell2mat(obj.FeatureInIndex(rr,h)));
                               
                % Get features that are in memory
                if numel(hm)
                    if isempty(obj.ReverseSubsetIndex)
                        % start offset for data in an ordered subset, if
                        % data is not subsetted yet, NumEntriesPerRefInSubset is
                        % equal to NumEntriesPerRefInSource anyways.
                        offset = uint32(sum(obj.NumEntriesPerRefInSubset(1:rr-1))); 
                        sil = si(w)-offset;
                    else
                        % if object was subsetted changing the order, we
                        % must use ReverseSubsetIndex
                        sil = obj.ReverseSubsetIndex(si(w));
                    end
                    for i = 1:numel(hm)
                        data{hm(i)}(w) = obj.FeatureData{rr,h(hm(i))}(sil);
                    end
                end
                
                % Translate the global indexing to a subsetted index (when
                % required) and then make it relative to the reference in
                % the source file. 
                if numel(hi) || numel(hf)
                    % start offset for data in source and index files
                    offset = uint32(sum(obj.NumEntriesPerRefInSource(1:rr-1))); 
                    if ~isempty(obj.SubsetIndex)
                        sil = obj.SubsetIndex(si(w)) - offset;
                    else
                        sil = si(w)-offset;
                    end
                end

                % Get features from BAL index file
                if numel(hi)
                    matobj = matfile(obj.LinearIndexFile);
                    if issorted(sil) && all(diff(sil))
                        % Requested indices in 'sil' are ordered and unique,
                        % then data can be copied directly in to varargout.
                        for i = 1:numel(hi)
                            data{hi(i)}(w) = readmatfile(matobj,obj.FeatureName{h(hi(i))},obj.ReferenceIndex(rr),sil);
                        end
                    else
                        % Requested indices in 'sil' are NOT ordered or
                        % unique, then first we need to figure out a index
                        % that minimizes the number of I/O operations and
                        % which traverses the file the most orderly
                        % possible. Data is copied to temporal variables
                        % and then reordered and arranged accordongly with
                        % the requested indexing.
                        [usi,~,sil2] = unique(sil);
                        for i = 1:numel(hi)
                            tmpdata = readmatfile(matobj,obj.FeatureName{h(hi(i))},obj.ReferenceIndex(rr),usi);
                            data{hi(i)}(w) = tmpdata(sil2(:));
                        end
                    end
                end % if numel(hi)
                
                % Get features from source file
                if numel(hf)
                    if ~obj.LinearIndexInMemory{rr}
                        % Only one vector of start positions at a time in
                        % memory:
                        obj.LinearIndexInMemory(:) = {false};
                        obj.LinearIndexInMemory{rr} = true;
                    end
                    
                    if issorted(sil) && all(diff(sil))
                        % Requested indices in 'sil' are ordered and unique,
                        % then data can be copied directly into varargout.
                        [tmpcell{1:numel(hf)}] = ...
                            bioinfoprivate.bamaccessmex('getfeaturesinrefbyindex',...
                            obj.InputFile,obj.IndexFile(1:end-4),obj.FeatureName(h(hf)),...
                            int32(obj.ReferenceIndex(rr)),sil,obj.LinearIndex{rr});
                        for i=1:numel(hf)
                            data{hf(i)}(w) = tmpcell{i};
                        end
                    else
                        % Requested indices in 'sil' are NOT ordered or
                        % unique, then first we need to figure out an index
                        % that minimizes the number of I/O operations and
                        % which traverses the file the most orderly
                        % possible. Data is copied to temporal variables
                        % and then reordered and arranged accordongly with
                        % the requested indexing.
                        [usi,~,sil(:)] = unique(sil);
                        tmpcell = cell(1,numel(hf));
                        [tmpcell{1:numel(hf)}] = ...
                            bioinfoprivate.bamaccessmex('getfeaturesinrefbyindex',...
                            obj.InputFile,obj.IndexFile(1:end-4),obj.FeatureName(h(hf)),...
                            int32(obj.ReferenceIndex(rr)),usi,obj.LinearIndex{rr});
                        for i=1:numel(hf)
                            data{hf(i)}(w) = tmpcell{i}(sil);
                        end
                    end
                end % if numel(hf)
                
            end
            
            % Fill in with the reference index when necessary
            if ~isempty(hr)
                data(hr+1:end+1) = data(hr:end);
                data{hr} = sr(:);
            end
            
            % Remove cell in case only one feature was requested
            if ~dataIsCell
                data = data{1};
            end
        end

    end % methods
    
    methods (Access = private)

        %==================================================================
        % indexToRefIndex
        %==================================================================
        function ri = indexToRefIndex(obj,si)
            %INDEXTOREFINDEX Returns the reference index for each global index.
            %
            % RI = INDEXTOREFINDEX(OBJ,IDX) Find the index of the reference
            % for each entry indexed by IDX. 
            
            if ~isempty(obj.SubsetIndex)
                si = obj.SubsetIndex(si);
            end
            v = uint32(obj.NumEntriesPerRefInSource);
            if ~isempty(v)
                s = v(1);
                ri = ones(size(si),'uint32');
                for i = 2:numel(v)
                    ri(si>s) = ri(si>s) + 1;
                    s = s + v(i);
                end
            else
                ri = [];
            end
        end
        %==================================================================
        % findFeature
        %==================================================================
        function [h,hr,name] = findFeature(obj,name)
            %FINDFEATURE Find the index in property FeatureName to a feature.
            %
            % [H, HR] = FINDFEATURE(OBJ,NAME) Finds the index in property
            % FeatureName to a feature and throws an error if absent. NAME
            % may also be a cell string with multiple non-repeated feature
            % names. HR is the index in NAME where the special case for the
            % feature named 'Reference' is found; 'reference is not indexed
            % in the BAM file and needs to be filled in manually.
            
            if ~iscell(name)
                name = {name};
            end
            hr = find(strcmp(name,'Reference'));
            if numel(hr)>1
                error(message('bioinfo:BAMIndexedFile:findFeature:nonUniqueFeature'))
            end
            name(hr) = [];    
            h = zeros(numel(name),1);
            for i = 1:numel(name)
              hh = find(strcmp(obj.FeatureName,name{i}),1);
              if isempty(hh)
                error(message('bioinfo:BAMIndexedFile:findFeature:invalidFeature', name{ i }))
              end  
              h(i) = hh;
            end
            if numel(unique(h))~=numel(h)
              error(message('bioinfo:BAMIndexedFile:findFeature:nonUniqueFeature'))
            end
        end
        
    end % methods (Access = private)
       
end % classdef

function d = readmatfile(mf,var,refi,si)

if nargin<4
    d = mf.(sprintf('%s_%d',var,refi));
elseif all(diff(si)==1)
    d = mf.(sprintf('%s_%d',var,refi))(si(:)',1);
else
    ti = min(si):max(si);
    td = mf.(sprintf('%s_%d',var,refi))(ti,1);
    d = td(si-si(1)+1);
end

end

function [reference,verbose] = parse_inputs(varargin)

%=== Check for the right number of inputs
if rem(numel(varargin),2) == 1
    error(message('bioinfo:BAMIndexedFile:BAMIndexedFile:IncorrectNumberOfArguments'))
end

%=== Defaults
reference = '';
verbose   = true;

%=== Allowed inputs
okargs = {'reference','verbose'};

for j = 1:2:numel(varargin)-1
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1   % reference
            reference = pval;
        case 2 % verbose
            verbose = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end
end
end
