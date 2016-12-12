function info = saminfo(filename, varargin)
%SAMINFO Information about SAM format file.
%
%   INFO = SAMINFO(FILENAME) returns a structure whose fields contain
%   information about a SAM file. FILENAME is a string containing a file
%   name, or a path and a file name, of a SAM file. INFO is a structure
%   with the following fields:
%
%       Filename    - name of the file
%       FilePath    - path to the file
%       FileModDate - modification date of the file
%       FileSize    - size of the file in bytes
%
%   Optional information contained in header section of SAM file.  Each of
%   these will be a structure with additional fields based on the
%   information present in the SAM file:
%
%       Header - File format version, sort order and group order
%
%       SequenceDictionary - Sequence Name, sequence length, genome
%       assembly identifier, MD5 checksum of sequence, URI of sequence,
%       species
%
%       ReadGroup - Read group identifier, sample, library, description,
%       platform unit, predicted median insert size, sequencing center,
%       date, platform
%
%       Program - Program name, version, command line
%
%       NOTE: Only the optional fields present in the SAM file will appear
%       in the output structure.
%
%   SAMINFO(...,'NUMOFREADS',T/F) scans the entire SAM file to determine
%   the number of alignment records stored in the file. The count is 
%   returned in the NumReads field. Default is false.
%
%   SAMINFO(...,'SCANDICTIONARY',T/F) scans the entire SAM file to
%   determine the available reference names and a tally for the alignment
%   records stored in the file. The information is returned in the
%   ScannedDictionary and ScannedDictionaryCount fields respectively.
%   Default is false.
%
%   Example:
%
%       % Get SAM file information.
%       info = saminfo('ex1.sam');
%
%       % Along with file information, also get number of alignment records
%       info = saminfo('ex1.sam','numofreads',true)
%
%   See also SAMREAD.

%   Copyright 2009-2011 The MathWorks, Inc.


if ~ischar(filename) && ~iscellstr(filename)
    error(message('bioinfo:saminfo:InvalidInput'))
end

%=== cellstr's are valid (UNDOCUMENTED), it is changed to a string with '\n'
%    to further work with it:
if iscellstr(filename)
    filename = sprintf('%s\n',filename{:});
end

%=== Is it a PADDED string ?
if size(filename,1) > 1
    if isvector(filename)
        filename = filename(:)';
    else
        error(message('bioinfo:saminfo:InvalidInput'))
    end
end

% Check is the string has a '\n' in the first 2^12 bytes, if so, it can
% not be a file name
singleLine = isempty(regexp(filename(1:min(end,2^12)),'\n','once'));

% Read input arguments:
[numReadsFlag, dictionaryFlag] = parse_inputs(varargin{:});

info = struct;

%=== read from a file
if singleLine && numel(filename)<(2^12) && exist(filename,'file')
    dataFromFile = true;
    % get the path
    filePath = fileparts(filename);
    % get full path+name for files that are in the MATLAB path:
    if isempty(filePath)
        filename = which(filename);
        filePath = fileparts(filename);
    end
    fileInfo = dir(filename);
    info.Filename = fileInfo.name;
    info.FilePath = filePath;
    info.FileSize = uint64(fileInfo.bytes);
    info.FileModDate = fileInfo.date;
    
    % Count the lines in the file, load only the first char of every line
    fid = fopen(filename,'rt');
    c = onCleanup(@()fclose(fid));
    
    % Find number of header rows:
    header_text = textscan(fid, '@%s', 'Delimiter', '\n');
    header_text = strtrim(header_text{1});
    
else %input must be a string
    dataFromFile = false;
    header_text = textscan(filename,'@%s','delimiter','\n');
    header_text = strtrim(header_text{1});
    c1 = textscan(filename,'%c%*s','delimiter','\n');
    % empty lines at the end are removed
    nl = regexp(c1{1}','\n','once');
    if ~isempty(nl)
        c1{1} = c1{1}(1:nl-1);
    end   
    n_dataRows = numel(c1{1});
end

if ~dataFromFile && n_dataRows <= 1 && all(filename~=sprintf('\t')) && all(filename~=sprintf('\n'))
    % Identify if the string in filename was an intended filename but
    % since it was not initially found the function assumed it was a
    % a string with '\n' and counted the lines anyways, if we counted
    % one single line we assume the user intended to pass a filename
    % and we give a more appropriate error:
    error(message('bioinfo:saminfo:FileNotFound', filename));
end

n_rows = numel(header_text);

defined_type = {'HD'; 'SQ'; 'RG'; 'PG'; 'CO'};
field_name = {'Header', 'SequenceDictionary', 'ReadGroup', 'Program', 'Comments'};

defined_tags{1}={'VN', 'SO', 'GO'};
defined_tags{2}={'SN', 'LN', 'AS', 'M5', 'UR', 'SP'};
defined_tags{3}={'ID', 'SM', 'LB', 'DS', 'PU', 'PI', 'CN', 'DT', 'PL'};
defined_tags{4}={'ID', 'VN', 'CL'};

subfields{1} = {'Version', 'SortOrder', 'GroupOrder'};
subfields{2} = {'SequenceName', 'SequenceLength', 'GenomeAssemblyID', 'MD5Checksum', 'URI', 'Species'};
subfields{3} = {'ID', 'Sample', 'Library', 'Description', 'PlatformUnit', 'PredictedMedianInsertSize', ...
    'SequencingCenter', 'Date', 'Platform'};
subfields{4} = {'ID', 'Version', 'CommandLine'};

try
    for row_idx = 1:n_rows
        fields = textscan(header_text{row_idx}, '%s', 'Delimiter', '\t');
        type = fields{1}{1};
        
        index = find(strncmp(type,defined_type,numel(type)));
        if isempty(index)
            warning(message('bioinfo:saminfo:InvalidHeaderType', type, header_text{ row_idx }))
            continue
        end
        if index==5
            value = fields{1}{2};
            if isfield(info, field_name{index})
                % cat multiple lines with comments
                info.(field_name{index}) = sprintf('%s %s',info.(field_name{index}),value);
            else
                info.(field_name{index}) = value;   
            end
            continue
        end
                   
        tags = regexp(fields{1}(2:end), ':', 'split', 'once');
        
        % When the TAG is mandatory and no-other TAG is mandatory, some
        % writers omit the TAG, e.g. ID tag in @PG or @RG (G743887)
        if numel(tags{1})==1 && any(strcmp(type,{'PG','RG'}))
            tags{1} = [{'ID'} tags{1}];
        end
        
        tags = [tags{:}]';
        value = tags(2:2:end);
        tags = tags(1:2:end);
        
        if numel(value)~=numel(tags)
            error(message('bioinfo:saminfo:ErrorInFile', row_idx))
        end
       
        [~, loc] = ismember(tags, defined_tags{index});
        i = find(~loc);
        for j = 1:numel(i)
            warning(message('bioinfo:saminfo:InvalidTagField', type, tags{ i( j ) }, value{ i( j ) }))
        end
        
        tags = subfields{index}(loc(loc>0));
        value = value(loc>0);
        
        makestruct = cell2struct(value, tags);
               
        if isfield(info, field_name{index})
            n_elements = numel(info.(field_name{index}));
            current_field = info.(field_name{index});
            for tag_idx = 1:numel(tags)
                current_field = setfield(current_field,{n_elements+1},tags{tag_idx}, value{tag_idx});
            end
            info.(field_name{index})=current_field;
        else
            info.(field_name{index})=makestruct;
        end
        
    end %for row_idx
    
    if isfield(info, 'SequenceDictionary')
        for struct_idx = 1:numel(info.SequenceDictionary)
            if ~isempty(info.SequenceDictionary(struct_idx).SequenceLength)
                info.SequenceDictionary(struct_idx).SequenceLength = str2double(info.SequenceDictionary(struct_idx).SequenceLength);
            end
        end
    end
    
catch allExceptions
    if strcmpi(allExceptions.identifier,'bioinfo:saminfo:ErrorInFile')
        rethrow(allExceptions);
    else
        error(message('bioinfo:saminfo:IncorrectDataFormat'));
    end
end

% Fill in with the information resulting from scanning the file
if dataFromFile
    info.NumReads = uint64([]);
    info.ScannedDictionary = cell(0,1);
    info.ScannedDictionaryCount = uint64(zeros(0,1));
    
    if numReadsFlag || dictionaryFlag
        [~,counts,refnames] = bioinfoprivate.bamaccessmex('saminfo',filename);
        if numReadsFlag
            info.NumReads = uint64(sum(counts));
        end
        if dictionaryFlag
            if counts(1)==0
                info.ScannedDictionary = refnames(2:end);
                info.ScannedDictionaryCount = uint64(counts(2:end));
            else
                info.ScannedDictionary = [refnames(2:end);{'Unmapped'}];
                info.ScannedDictionaryCount = uint64(counts([2:end 1]));
            end
        end
    end
else
    if numReadsFlag
        info.NumReads = uint64(n_dataRows);
    end %if numReadsFlag
    if dictionaryFlag
        error(message('bioinfo:saminfo:CellstrInputUndocumented'))
    end    
end

%--------------------------------------------------------------------------
function [numReadsFlag, dictionaryFlag] = parse_inputs(varargin)
% Parse input PV pairs.

% defaults
numReadsFlag = 0;
dictionaryFlag = 0;

if rem(nargin,2) ~= 0
    error(message('bioinfo:saminfo:IncorrectNumberOfArguments', mfilename));
end

okargs = {'numofreads','scandictionary'};

for j=1:2:nargin-1
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    if k == 1
        numReadsFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    elseif k==2
        dictionaryFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end
end
