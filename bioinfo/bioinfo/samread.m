function [data header] = samread(filename, varargin)
%SAMREAD reads SAM format file.
%
%   DATA = SAMREAD(FILENAME) reads a SAM format file FILENAME, returning
%   the data in the file as a structure. FILENAME can also be a MATLAB
%   string that contains the text of a SAM format file. The output
%   DATA will be an Nx1 structure, where N is the number of alignment
%   records stored in the SAM file.
%
%   [DATA, HEADER] = SAMREAD(FILENAME) returns the alignment data as well
%   as the header information stored in the file.
%
%   SAMREAD(...,'TAGS',false) reads the first eleven fields for each
%   alignment in the SAM file, but does not read any of the optional tags.
%   Default is true.
%
%   SAMREAD(...,'READGROUP',RG) reads only the records belonging to the
%   specified read group. RG is a MATLAB string. By default SAMREAD reads all the
%   groups.
%
%   SAMREAD(...,'BLOCKREAD', M) allows you to read in a single entry or
%   block of entries from a file containing multiple alignment records. If
%   M is a scalar then the M'th entry in the file is read. If M is a two
%   element vector then the block of entries starting at entry M(1) and
%   ending at entry M(2) will be read.  Use Inf for M(2) to read all
%   entries in the file starting at position M(1).
%
%   Examples:
%
%       % Read the alignment records stored in 'ex1.sam'
%       [data header] = samread('ex1.sam');
%
%       % Read a block of entries from a file
%       data = samread('ex1.sam','blockread', [ 5 10], 'tags', false);
%
%       % Read records belonging to the read group 'L1'
%       data = samread('ex1.sam', 'readgroup', 'L1');
%
%   See also SAMINFO, EMBLREAD, FASTAINFO, FASTAWRITE, FASTQINFO, FASTQREAD,
%   FASTQWRITE, GENBANKREAD, GENPEPTREAD, HMMPROFDEMO, MULTIALIGNREAD,
%   SEQALIGNVIEWER, SEQPROFILE, SEQVIEWER, SFFINFO, SFFREAD.

%   Copyright 2009-2012 The MathWorks, Inc.

% SAM format specified here:
% http://samtools.sourceforge.net/SAM1.pdf

if ~ischar(filename) || size(filename, 1)>1
    error(message('bioinfo:samread:InvalidInput'));
end

% Parse input
[blockRead,blockSize,skip,readTagsFlag,filterByReadGroup,readgroup] = parse_inputs(varargin{:});

if  isempty(regexp(filename,'\t','once')) && ((exist(filename,'file') || exist(fullfile(pwd,filename),'file')))
    fid = fopen(filename, 'r');
    c = onCleanup(@()fclose(fid));
    
    %Read header text
    header = textscan(fid, '@%s', 'Delimiter', '\n');
    header = header{1};
    
    if blockRead
        if isinf(skip)
            ftext = cell(0,1); % SAMREAD(...,'BLOCKREAD', Inf)
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
            if any(cellfun(@isempty,ftext))
                ftext = ftext(~cellfun(@isempty,ftext));
            end            
        end
        if isempty(ftext)
            error(message('bioinfo:samread:StartTooBig',sprintf('%d',skip + 1)));
        end
    else
        try
            ftext = textscan(fid,'%s','Delimiter','\n');
            ftext = strtrim(ftext{1});
            if any(cellfun(@isempty,ftext))
                ftext = ftext(~cellfun(@isempty,ftext));
            end
        catch theErr
            if strcmpi(theErr.identifier,'MATLAB:nomem')
                error(message('bioinfo:samread:FileTooBig'));
            else
                rethrow(theErr);
            end
        end
    end
else  % input must be a string
    %assume that if the string does not contain a delimiter, it is intended
    %to be a filename
    if all(filename~=sprintf('\t')) && all(filename~=sprintf('\n'))
        error(message('bioinfo:samread:FileNotFound', filename));
    end
    
    header = textscan(filename,'@%s','Delimiter','\n');
    header = header{1};
    n_headerrows = numel(header);
    
    if blockRead
        if isinf(skip)
            ftext = cell(0,1); % SAMREAD(...,'BLOCKREAD', Inf)
        else
            if isinf(blockSize)
                ftext = textscan(filename, '%s','Delimiter', '\n', 'HeaderLines', n_headerrows+skip); % read everything
            else
                ftext = textscan(filename,'%s', blockSize,'Delimiter','\n', 'HeaderLines', n_headerrows+skip);
            end
            ftext = strtrim(ftext{1});
            if any(cellfun(@isempty,ftext))
                ftext = ftext(~cellfun(@isempty,ftext));
            end
        end
        if isempty(ftext)
            error(message('bioinfo:samread:StartTooBigInString',sprintf('%d',skip + 1)));
        end
    else
        ftext = textscan(filename,'%s','Delimiter','\n', 'HeaderLines', n_headerrows);
        ftext = strtrim(ftext{1});
        if any(cellfun(@isempty,ftext))
            ftext = ftext(~cellfun(@isempty,ftext));
        end        
    end
end %if exist(filename)

if filterByReadGroup
    readgrouptag = sprintf('\tRG:Z:%s', readgroup);
    belongsToGroup = regexp(ftext, [readgrouptag '(?:\t|\>)'], 'once');
    belongsToGroup = cellfun(@(x)~isempty(x), belongsToGroup);
    ftext = ftext(belongsToGroup);
    if isempty(ftext)
        if blockRead
            error(message('bioinfo:samread:InvalidReadGroupInBlock', readgroup));
        else
            error(message('bioinfo:samread:ReadGroupNotFound', readgroup, filename));
        end
    end
end

try
    if nargout>1
        header = saminfo(strcat('@',header));
        if isempty(fieldnames(header))
            header = [];
        end
    end
    
    if readTagsFlag
        s = struct('QueryName', [],'Flag', [],'ReferenceName',[], 'Position',[], 'MappingQuality',[],...
            'CigarString',[], 'MateReferenceName',[], 'MatePosition',[], 'InsertSize',[], 'Sequence',[], 'Quality', [], 'Tags', []);
    else
        s = struct('QueryName', [],'Flag', [],'ReferenceName',[], 'Position',[], 'MappingQuality',[], ...
            'CigarString',[], 'MateReferenceName',[], 'MatePosition',[], 'InsertSize',[], 'Sequence',[], 'Quality', []);
    end
    
    data = repmat(s,numel(ftext),1);
    
    for row = 1:numel(ftext)
        [rtext pos] = textscan(ftext{row}, '%s%u16%s%u32%u8%s%s%u32%d32%s%s',1, 'Delimiter', '\t');
        
        if isempty(rtext{11}) % when textscan fails to read one of the fields they are returned empty
            error(message('bioinfo:samread:IncompleteEntry', row))
        else
            % index into cell and extract data
            rtext([1 3 6 7 10 11])=[rtext{[1 3 6 7 10 11]}];
        end
        
        if readTagsFlag
            
            optional_fields = regexp(ftext{row}(pos+1:end), '([^\s:]{2}):([AifZH]):(\S+)\t?', 'tokens');
            
            if ~isempty(optional_fields)
                optional_fields = [optional_fields{:}]';
                tags = optional_fields(1:3:end);
                type = [optional_fields{2:3:end}];
                value = optional_fields(3:3:end);
                
                integers = type == 'i';
                if sum(integers)
                    s = {' '};
                    spaces = s(ones(1, sum(integers)));
                    tmp = [value(integers)'; spaces];
                    vals = [tmp{:}];
                    value(integers) = num2cell(int32(sscanf(vals,'%d')));
                end
                
                floats = type == 'f';
                if sum(floats)
                    s = {' '};
                    spaces = s(ones(1, sum(floats)));
                    tmp = [value(floats)'; spaces];
                    vals = [tmp{:}];
                    value(floats) = num2cell(single(sscanf(vals,'%f')));
                end
                
                hexstrings = strfind(type, 'H');
                if ~isempty(hexstrings)
                    for hex_idx = 1:numel(hexstrings)
                        value{hexstrings(hex_idx)} =  uint8(hex2dec(reshape(value{hexstrings(hex_idx)}, 2, [])')');
                    end
                end
                
                if numel(value)~=numel(tags)
                    error(message('bioinfo:samread:ErrorInFile', n_headerrows + row))
                end
                
                data(row) = cell2struct([rtext'; {cell2struct(value, tags)}], {'QueryName', 'Flag', 'ReferenceName', 'Position', ...
                    'MappingQuality', 'CigarString', 'MateReferenceName', 'MatePosition', 'InsertSize', 'Sequence', 'Quality', 'Tags'});
                
            else
                data(row) = cell2struct([rtext';{[]}], {'QueryName', 'Flag', 'ReferenceName', 'Position', 'MappingQuality', 'CigarString', ...
                    'MateReferenceName', 'MatePosition', 'InsertSize', 'Sequence', 'Quality', 'Tags'});
            end %isempty(optional_fields)
        else
            data(row) = cell2struct(rtext', {'QueryName', 'Flag', 'ReferenceName', 'Position', 'MappingQuality', 'CigarString', ...
                'MateReferenceName', 'MatePosition', 'InsertSize', 'Sequence', 'Quality'});
        end %readTagsFlag
    end %for row
catch ME
    if (strfind(ME.identifier, 'bioinfo:samread:'))
        rethrow(ME)
    elseif (strfind(ME.identifier, 'MATLAB:nomem'))
        error(message('bioinfo:samread:FileTooBig'));
    else
        error(message('bioinfo:samread:IncorrectDataFormat'));
    end
end

%--------------------------------------------------------------------------
function [blockRead,blockSize,skip,readTagsFlag,filterByReadGroup,readgroup] = parse_inputs(varargin)
% Parse input PV pairs.

% default value
blockRead = false;
blockSize = inf;
skip = 0;
readTagsFlag = true;
filterByReadGroup = false;
readgroup = [];

if rem(nargin,2) ~= 0
    error(message('bioinfo:samread:IncorrectNumberOfArguments', mfilename));
end
okargs = {'blockread', 'tags', 'readgroup'};
for j=1:2:nargin-1
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % range
            if ~isnumeric(pval) || numel(pval)> 2 || ...
                    isempty(pval) || any(pval<1) || any(rem(pval,1))
                error(message('bioinfo:samread:BadBlockRange'))
            end
            blockRead = true;
            pval = double(pval);
            
            skip = min(pval)-1;
            blockSize = max(pval)-skip;
        case 2 %tags
            readTagsFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 3 %read only a specific read group
            filterByReadGroup = true;
            readgroup = pval;
            if ~ischar(readgroup) || size(readgroup, 1)>1
                error(message('bioinfo:samread:InvalidReadGroup'));
            end
    end
end
