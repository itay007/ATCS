function data = bowtieread(filename, varargin)
%BOWTIEREAD reads BOWTIE format file.
%
%   NOTE:  BOWTIEREAD will be removed in a future release. When using the
%   BOWTIE mapper/aligner make sure that you specify the appropriate option
%   in order to create either a SAM or BAM output file, then use the BioMap
%   class or the SAMREAD/BAMREAD functions to access the mapped short reads.
%
%   BOWTIEREAD(FILENAME) reads a BOWTIE format file FILENAME, returning the
%   data in the file as a MATLAB structure. FILENAME can also be a string
%   that contains the text of a BOWTIE format file. The output is an Nx1
%   structure, where N is the number of alignment records stored in the
%   BOWTIE file.
%
%   BOWTIEREAD(...,'BLOCKREAD', M) allows you to read in a single entry or
%   block of entries from a file containing multiple alignment records. If
%   M is a scalar then the M'th entry in the file is read. If M is a two
%   element vector then the block of entries starting at entry M(1) and
%   ending at entry M(2) will be read.  Use Inf for M(2) to read all
%   entries in the file starting at position M(1).
%
%   BOWTIEREAD(...,'ZEROBASED', true) returns the position of the mapped
%   queries using 0-based indexing. By default (false) the position is
%   1-based indexed.
%
%   BOWTIEREAD(...,'ALIGNDETAILS', false) reads the first seven fields for
%   each record in the BOWTIE file, but does not read the mismatch
%   descriptors. Default is true.
%
%   Examples:
%
%       % Read the alignment records stored in 'sample01.bowtie'
%       data = bowtieread('sample01.bowtie');
%
%       % Read a block of entries from a file
%       data = bowtieread('sample01.bowtie','blockread', [5 10]);
%
%   See also SAMREAD, FASTQREAD, SOAPREAD

%   Copyright 2010-2012 The MathWorks, Inc.

% For a detailed description of the BOWTIE alignment format, see
% http://bowtie-bio.sourceforge.net/manual.shtml

warning(message('bioinfo:bowtieread:incompatibility'));
 
bioinfochecknargin(nargin,1,mfilename);
if ~ischar(filename)
    error(message('bioinfo:bowtieread:InvalidInput'));
end

% Parse inputs
[blockRead, blockSize, skip, offset, aligndetails] = parse_inputs(varargin{:});

if  (exist(filename,'file') || exist(fullfile(pwd,filename),'file'))
    fid = fopen(filename, 'r');
    c = onCleanup(@()fclose(fid));
    if blockRead
        if isinf(skip)
            ftext = cell(0,1); % BOWTIEREAD(...,'BLOCKREAD', Inf)
        else
            if skip>0
                textscan(fid,'%*s',skip,'Delimiter','\n'); % jump records
            end
            if isinf(blockSize)
                ftext = textscan(fid, '%s','Delimiter', '\n'); % read everything
            else
                ftext = textscan(fid,'%s', blockSize,'Delimiter','\n');
            end
            ftext = strtrim(ftext{1});
        end
        if isempty(ftext)
            error(message('bioinfo:bowtieread:StartTooBig', sprintf('%d',skip + 1)));
        end
    else
        try
            ftext = textscan(fid,'%s','Delimiter','\n');
            ftext = strtrim(ftext{1});
            
        catch theErr
            if strcmpi(theErr.identifier,'MATLAB:nomem')
                error(message('bioinfo:bowtieread:FileTooBig'));
            else
                rethrow(theErr);
            end
        end
    end
else  % input must be a string
    %assume that if the string does not contain a delimiter, it is intended
    %to be a filename
    if all(filename~=sprintf('\t'))
        error(message('bioinfo:bowtieread:FileNotFound', filename));
    end
    ftext = textscan(filename,'%s','Delimiter','\n');
    ftext = strtrim(ftext{1});
end

try
    % Initialize output structure
    % AlignDetails contains information on Mismatches, Insertions and
    % Deletions if any available
    if aligndetails
        data(numel(ftext),1) = struct('QueryName', [],'Strand', [],'ReferenceName',[], 'Position',[], ...
            'Sequence',[], 'Quality', [], 'NumHits', [], 'AlignDetails', []);
    else
        data(numel(ftext),1) = struct('QueryName', [],'Strand', [],'ReferenceName',[], 'Position',[], ...
            'Sequence',[], 'Quality', [], 'NumHits', []);
    end
    
    for row = 1:numel(ftext)
        % Read the text and parse each line
        if aligndetails
            rtext = textscan(ftext{row}, '%s%s%s%u32%s%s%u32%s',1, 'Delimiter', '\t');
        else
            rtext = textscan(ftext{row}, '%s%s%s%u32%s%s%u32',1, 'Delimiter', '\t');
        end

        if isempty(rtext{7}) % when textscan fails to read one of the fields they are returned empty
            error(message('bioinfo:bowtieread:IncompleteEntry', row))
        else
            % index into cell and extract data
            rtext([1 2 3 5 6]) = [rtext{[1 2 3 5 6]}];
        end
        if (offset) % 1-based indexing
            rtext{4} = rtext{4}+1;
        end
        
        if aligndetails
            if ~isempty(rtext{8})
                rtext(8)=[rtext{8}];
            else
                rtext{8} = '';
            end
            data(row) = cell2struct(rtext', {'QueryName', 'Strand', 'ReferenceName', 'Position', ...
                'Sequence', 'Quality', 'NumHits', 'AlignDetails'});
        else
            data(row) = cell2struct(rtext', {'QueryName', 'Strand', 'ReferenceName', 'Position', ...
                'Sequence', 'Quality', 'NumHits'});
        end
    end %for row
catch ME
    if (strfind(ME.identifier, 'bioinfo:bowtieread:IncompleteEntry'))
        rethrow(ME)
    elseif (strfind(ME.identifier, 'MATLAB:nomem'))
        error(message('bioinfo:bowtieread:FileTooBig'));
    else
        % Generic Error
        error(message('bioinfo:bowtieread:IncorrectDataFormat'));
    end
end

%--------------------------------------------------------------------------
function [blockRead, blockSize, skip, offset, aligndetails] = parse_inputs(varargin)
% Parse input PV pairs.

% default value
blockRead = false;
offset = 1;
aligndetails = true;
blockSize = inf;
skip = 0;

% get input arguments
if rem(nargin,2) ~= 0
    error(message('bioinfo:bowtieread:IncorrectNumberOfArguments', mfilename));
end
okargs = {'blockread', 'zerobased', 'aligndetails'};
for j=1:2:nargin-1
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % range
            if ~isnumeric(pval) || numel(pval)> 2 || ...
                    isempty(pval) || any(pval<1) || any(rem(pval,1))
                error(message('bioinfo:bowtieread:BadBlockRead'))
            end
            blockRead = true;
            pval = double(pval);
            if numel(pval)==1
                skip = pval-1;
                blockSize = 1;
            else
                skip = min(pval)-1;
                blockSize = max(pval)-skip;
            end
        case 2 % format
            offset = ~(bioinfoprivate.opttf(pval,okargs{k},mfilename));
        case 3 %aligndetails
            aligndetails = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end
end

