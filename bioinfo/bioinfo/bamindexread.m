function index = bamindexread(filename)
%BAMINDEXREAD reads BAM.BAI files
%
%   INDEX = BAMINDEXREAD(FILENAME) returns a structure which specifies the
%   offsets into the compressed BAM file and uncompressed data block for
%   each reference and range of positions (bins) on the reference.  This
%   information is used by BAMREAD for indexing a BAM file to extract
%   alignment records in a specified range.  FILENAME is a string
%   containing a file name, or a path and a file name, of a BAM file or the
%   corresponding BAM.BAI file. INDEX is a structure with the following
%   fields:
%
%   Filename - Name of the BAM.BAI file used to generate the index
%   structure
%
%   Index - 1xNRefs array of structures with the following fields
%
%           BinID - array of bin IDs
%
%           BGZFOffsetStart - offset to the start of the first BGZF block
%           in the BAM file where alignment records with the corresponding
%           bin ID are stored
%
%           BGZFOffsetEnd - offset to the start of the last BGZF block in
%           the BAM file where alignment records with the corresponding
%           bin ID are stored
%
%           DataOffsetStart - start offset in decompressed data block where
%           alignment records with the corresponding bin ID are stored
%
%           DataOffsetEnd - end offset in decompressed data block where
%           alignment records with the corresponding bin ID are stored
%
%           LinearBGZFOffset - Offset in BAM file to the first alignment in
%           the corresponding 16384bp interval
%
%           LinearDataOffset - Offset in decompressed data block to the
%           first alignment in the corresponding 16384bp interval
%
%   Example:
%
%       % Read BAM index structure and use it to read multiple alignment
%       % records from a BAM file
%       ind = bamindexread('ex1.bam');
%       data1 = bamread('ex1.bam', 'seq1', [100 200], 'index', ind);
%       data2 = bamread('ex1.bam', 'seq2', [100 200], 'index', ind);
%
%   See also BAMREAD, BAMINFO.

%   Copyright 2009-2010 The MathWorks, Inc.


% BAM and BAM.BAI file formats specified here:
% http://samtools.sourceforge.net/SAM1.pdf

if ~ischar(filename) || size(filename, 1) > 1
    error(message('bioinfo:bamread:InvalidInput'))
end
if strcmp(filename(end-2:end), 'bam')
    filename = [filename '.bai'];
end
if ~((exist(filename,'file') || exist(fullfile(pwd,filename),'file')))
    error(message('bioinfo:bamread:FileNotFound', filename));
end

fid = fopen(filename, 'r', 'l');
c = onCleanup(@()fclose(fid));

magic = fread(fid, 4, '*char');

if ~strcmp(magic', sprintf('BAI\1'))
    error(message('bioinfo:bamindexread:InvalidMagicNumber'));
end

n_ref = fread(fid, 1, '*int32');

ind(1, n_ref) = struct('BinID', [], ...
    'BGZFOffsetStart', [], 'BGZFOffsetEnd', [], ...
    'DataOffsetStart', [], 'DataOffsetEnd', [], ...
    'LinearBGZFOffset', [], 'LinearDataOffset', []);

for i_ref = 1: n_ref
    
    n_bins = fread(fid, 1, 'int32');
    
    if n_bins
        bins(n_bins) = struct('BinID', [], ...
            'BGZFOffsetStart', [], 'BGZFOffsetEnd', [], ...
            'DataOffsetStart', [], 'DataOffsetEnd', []); %#ok<AGROW>
        
        for i_bin = 1:n_bins
            
            binID = fread(fid, 1, '*uint32');
            n_chunk = fread(fid, 1, '*int32');
            
            bins(i_bin).BinID = binID * uint32(ones(1, n_chunk)); %#ok<AGROW>
            bins(i_bin).DataOffsetStart = uint16(zeros(1, n_chunk)); %#ok<AGROW>
            bins(i_bin).BGZFOffsetStart = uint64(zeros(1, n_chunk)); %#ok<AGROW>
            bins(i_bin).DataOffsetEnd = uint16(zeros(1, n_chunk)); %#ok<AGROW>
            bins(i_bin).BGZFOffsetEnd = uint64(zeros(1, n_chunk)); %#ok<AGROW>
            
            for i_chunk = 1:n_chunk
                bins(i_bin).DataOffsetStart(i_chunk) = fread(fid, 1, '*ubit16'); %#ok<AGROW>
                bins(i_bin).BGZFOffsetStart(i_chunk) = fread(fid, 1, '*ubit48'); %#ok<AGROW>
                bins(i_bin).DataOffsetEnd(i_chunk) = fread(fid, 1, '*ubit16'); %#ok<AGROW>
                bins(i_bin).BGZFOffsetEnd(i_chunk) = fread(fid, 1, '*ubit48'); %#ok<AGROW>
            end
        end
        
        ind(i_ref).BinID = [bins(1:n_bins).BinID];
        ind(i_ref).DataOffsetStart = [bins(1:n_bins).DataOffsetStart];
        ind(i_ref).DataOffsetEnd = [bins(1:n_bins).DataOffsetEnd];
        ind(i_ref).BGZFOffsetStart = [bins(1:n_bins).BGZFOffsetStart];
        ind(i_ref).BGZFOffsetEnd = [bins(1:n_bins).BGZFOffsetEnd];
    end
    
    n_intv = fread(fid, 1, '*int32');
    
    if ~isempty(n_intv)
        ind(i_ref).LinearDataOffset = fread(fid, n_intv, '*ubit16', 48)';
        fseek(fid, -8*n_intv+2, 'cof');
        ind(i_ref).LinearBGZFOffset = fread(fid, n_intv, '*ubit48', 16)';
        fseek(fid, -2, 'cof');
    end
end

index.Filename = filename;
index.Index = ind;
