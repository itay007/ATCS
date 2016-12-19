function varargout = bamread(filename, reference, range, varargin)
%BAMREAD reads BAM format file.
%
%   DATA = BAMREAD(FILENAME, REFERENCE, RANGE) reads the alignment records
%   which align to the reference sequence REFERENCE in the range RANGE in a
%   BAM format file FILENAME, returning the data in the file as a
%   structure.  The corresponding BAM index (BAM.BAI) file is required for
%   reading the BAM format file.  REFERENCE may either be a MATLAB string
%   containing the name of the reference sequence, or a positive integer R,
%   indicating the Rth reference sequence in the BAM format file. R is also
%   the index of the reference in the info.Reference structure returned by
%   the BAMINFO function. RANGE is a two element vector specifying the
%   first and last position of the range on the reference sequence. RANGE
%   uses 1-based positions. The output DATA will be an Nx1 array of
%   structures, where N is the number of records which align to the
%   requested range.
%
%   [DATA, HEADER] = BAMREAD(FILENAME, REFERENCE, RANGE) returns the
%   alignment data as well as the header information stored in the file.
%
%   BAMREAD(..., 'FULL', TRUE) only returns those alignment records which
%   are fully contained in the specified range. Default is false.
%
%   BAMREAD(..., 'TAGS', false) reads the first eleven fields for each
%   alignment record, but does not read any of the optional tags.  Default
%   is true.
%
%   BAMREAD(..., 'TOFILE', FILENAME) stores the alignment records in the
%   requested range in a SAM format file. FILENAME is a string specifying a
%   file name or a path and file name to which to save the data. If you
%   specify only a file name, the file is stored in the current folder.  If
%   the data you are trying to read does not fit in memory, use this option
%   with no output arguments. You can then access records in the resulting
%   SAM formatted file using the BioIndexedFile class, or you can access
%   the mapped records in the original BAM formatted file using the BioMap
%   class. 
%
%   BAMREAD(..., 'ZEROBASED', TRUE) uses 0-based positions in the output
%   DATA structure. By default, ZEROBASED is false,  BAMREAD uses 1-based
%   positions. The SAM format file created with the 'TOFILE' option
%   contains 1-based positions, regardless of the value of 'ZEROBASED'.
%
%   BAMREAD(FILENAME, 'Unmapped') or BAMREAD(FILENAME, 0) returns reads
%   that are not mapped to any reference.
%
%   Examples:
%
%       % Read alignment records stored in 'ex1.bam'
%       [data header] = bamread('ex1.bam', 'seq1', [100 200]);
%
%       % Read only alignment records fully contained in the range [100 200]
%       data = bamread('ex1.bam', 'seq1', [100 200], 'full', true);
%
%       % Store output in a SAM format file; do not read optional tags
%       bamread('ex1.bam','seq1', [1 500], 'tags', false, ...
%               'tofile', 'ex1_example.sam');
%
%
%   See also BIOMAP, BAMINFO, BAMINDEXREAD, SAMREAD, SAMINFO, FASTAINFO,
%   FASTAWRITE, FASTQINFO, FASTQREAD, FASTQWRITE.

%   Copyright 2009-2012 The MathWorks, Inc.

% BAM format specified here:
% http://samtools.sourceforge.net/SAM1.pdf

if ~ischar(filename) || size(filename, 1) > 1
    error(message('bioinfo:bamread:InvalidInput'))
elseif ~((exist(filename,'file') || exist(fullfile(pwd,filename),'file')))
    error(message('bioinfo:bamread:FileNotFound', filename));
end

% Check file extension
[pa,~,theExt] = fileparts(filename);
if ~strcmp(theExt,'.bam')
    error(message('bioinfo:bamread:InvalidExtension'));
end

% it may be possible that the user is expecting to find a file in the
% Matlab path, if so, the path is appended to the filename when it is
% different from the current directory
if isempty(pa) && ~strcmpi(fileparts(which(filename)),pwd)
    filename = which(filename);
end

% Check if signature is BAMREAD(FILENAME, 'Unmapped') or BAMREAD(FILENAME, 0)
if (ischar(reference) && strcmpi('Unmapped',reference)) || (isnumeric(reference) && isscalar(reference) && (reference==0))
    if nargin>2
        varargin = [{range} varargin];
    end
    readUnmapped = true;
else
    readUnmapped = false;
end

if ~readUnmapped
    % Check if BAI file exists, if not it indexes the BAM file
    if ~numel(dir([filename '.bai']))
        if ~bioinfoprivate.bamaccessmex('baminfo',filename)
            error(message('bioinfo:bamread:UnorderedBAMFile'))
        end
        try
            bioinfoprivate.bamaccessmex('indexbamfile',filename,[filename '.bai']);
        catch ME %#ok<NASGU>
            % When there is an error in the indexer it is likely that
            % the file left is invalid, then it is deleted.
            if numel(dir([filename '.bai']))
                delete(indexFile)
            end
            error(message('bioinfo:bamread:FileNotIndexable'))
        end
    end
    
    % Check RANGE input (RANGE must be 1-based)
    if ~isnumeric(range) || numel(range)> 2 || isempty(range) || ~all(range>0)
        error(message('bioinfo:bamread:BadRange'));
    end
    
    % Uniformize RANGE
    range = sort(range);
    if numel(range)<2
        range = [range range];
    end
    range = double(range); % Range is 1-based
    
    if range(2)>536870911 % 2^29-1
        error(message('bioinfo:bamread:BadRangeUpperLimit'));
    end
    
    % Obtain the header information:
    header = rmfield(baminfo(filename),{'Filename','FilePath','FileSize','FileModDate'});
    dic = {header.SequenceDictionary.SequenceName};
    
    % Check REFERENCE input
    if ischar(reference)
        i = find(strcmp(reference,dic),1);
        if isempty(i)
            error(message('bioinfo:bamread:ReferenceNotFound',reference,filename))
        end
        reference = i;
    end
    
    if ~(isnumeric(reference) && isscalar(reference) && reference > 0) 
        error(message('bioinfo:bamread:BadReference'));
    end
    
    if  reference > numel(dic)
        error(message('bioinfo:bamread:NumericReferenceNotFound',reference,numel(dic),filename));
    end
        
    % Standarize reference
    reference = dic{reference};
end

[full zeroBasedPos samfile readTags] = parse_inputs(varargin{:});

if isempty(samfile)
    outfile = [tempname,'.sam'];
else 
    outfile = samfile;
end

% Copy header to outfile   
headertext = bioinfoprivate.bamaccessmex('getheader',filename);
fid = fopen(outfile,'w');
if fid<0
    error(message('bioinfo:bamread:CanNotCreateSAMFile'));
elseif isempty(samfile)
    c = onCleanup(@()delete(outfile));
end
fprintf(fid,'%s',headertext);
fclose(fid);

if readUnmapped
    bioinfoprivate.bamaccessmex('writesamunaligned',filename,outfile,readTags);
    if nargout>0
        varargout{1} = samread(outfile,'TAGS',readTags);
        if numel(varargout{1})==0
            if readTags
                varargout{1} = repmat(struct('QueryName',[],'Flag',[],'Position',[],...
                    'MappingQuality',[],'CigarString',[],'MatePosition',[],...
                    'InsertSize',[],'Sequence',[],'Quality',[],'Tags',[],...
                    'ReferenceIndex',[],'MateReferenceIndex',[]),0,1);
            else
                varargout{1} = repmat(struct('QueryName',[],'Flag',[],'Position',[],...
                    'MappingQuality',[],'CigarString',[],'MatePosition',[],...
                    'InsertSize',[],'Sequence',[],'Quality',[],...
                    'ReferenceIndex',[],'MateReferenceIndex',[]),0,1);
            end
        else
            for i = 1:numel(varargout{1})
                varargout{1}(i).ReferenceIndex = int32(0);
                varargout{1}(i).MateReferenceIndex = int32(0);
            end
            varargout{1} = rmfield(varargout{1},{'ReferenceName','MateReferenceName'});
        end
    end
else
    bioinfoprivate.bamaccessmex('writesambyregion',filename,filename,...
        outfile,sprintf('%s:%d-%d',reference,range(1),range(2)),readTags,full)
    if nargout>0
        varargout{1} = samread(outfile,'TAGS',readTags);
        % For backwards compatibility, calculate the reference indices for the
        % read and its mate instead of returning the full name.
        if numel(varargout{1})==0
            if readTags
                varargout{1} = repmat(struct('QueryName',[],'Flag',[],'Position',[],...
                    'MappingQuality',[],'CigarString',[],'MatePosition',[],...
                    'InsertSize',[],'Sequence',[],'Quality',[],'Tags',[],...
                    'ReferenceIndex',[],'MateReferenceIndex',[]),0,1);
            else
                varargout{1} = repmat(struct('QueryName',[],'Flag',[],'Position',[],...
                    'MappingQuality',[],'CigarString',[],'MatePosition',[],...
                    'InsertSize',[],'Sequence',[],'Quality',[],...
                    'ReferenceIndex',[],'MateReferenceIndex',[]),0,1);
            end
        else
            dic = [{'*'} dic];
            ri = int32(seqmatch({varargout{1}.ReferenceName},dic,'exact',true))-1;
            for i = 1:numel(varargout{1})
                varargout{1}(i).ReferenceIndex = ri(i);
            end
            dic = [{'='} dic]; % symbol used to indicate that reference of mate is the same
            ri2 = int32(seqmatch({varargout{1}.MateReferenceName},dic,'exact',true))-2;
            ri2(ri2<0) = ri(ri2<0);
            for i = 1:numel(varargout{1})
                varargout{1}(i).MateReferenceIndex = ri2(i);
            end
            varargout{1} = rmfield(varargout{1},{'ReferenceName','MateReferenceName'});
            
            if zeroBasedPos
                for i = 1:numel(varargout{1})
                    varargout{1}(i).Position = varargout{1}(i).Position-1;
                    varargout{1}(i).MatePosition = varargout{1}(i).MatePosition-1;
                end
            end
        end
    end
end

if nargout>1
    varargout{2} = saminfo(outfile);
    varargout{2} = rmfield(varargout{2},{'Filename','FilePath','FileSize','FileModDate'});
end

%--------------------------------------------------------------------------
function [full zeroBasedPos samfile readTags] = parse_inputs(varargin)

% default values
full = 0;
zeroBasedPos = 0;
readTags = 1;
samfile = [];

% get input arguments
if rem(nargin,2) == 1
    error(message('bioinfo:bamread:IncorrectNumberOfArguments', mfilename));
end
okargs = {'full', 'index', 'zerobased', 'tofile', 'tags'};
for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1 %full
            full = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 2 %bamindex structure
            error(message('bioinfo:bamread:Incompatibility'))
        case 3 %0-based positions
            zeroBasedPos = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 4 %write data to sam file
            if ~ischar(pval) || size(pval, 1) > 1
                error(message('bioinfo:bamread:InvalidFilename'));
            end
            samfile = pval;
            if exist(samfile, 'file') == 2
                error(message('bioinfo:bamread:ExistingFile', samfile));
            end
        case 5 %read tags
            readTags = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end
end
