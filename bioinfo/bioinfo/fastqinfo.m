function info = fastqinfo(filename)
%FASTQINFO return a summary of the contents of a FASTQ file.
%
%   INFO = FASTQINFO(FILENAME) returns a structure whose fields contain
%   information about a FASTQ file.  FILENAME is a string containing a file
%   name, or a path and a file name, of a FASTQ file. INFO is a structure
%   with the following fields:
%
%             Filename - name of the file
%             FilePath - path to the file
%          FileModDate - modification date of the file
%             FileSize - size of the file in bytes
%      NumberOfEntries - number of sequence reads
%
%   Example:
%   % Show a summary of the contents of a FASTQ file.
%   info = fastqinfo('SRR005164_1_50.fastq')
%
%   See also FASTAINFO, FASTAREAD, FASTAWRITE, FASTQREAD, FASTQWRITE,
%   SFFINFO, SFFREAD.

%   Copyright 2009-2010 The MathWorks, Inc.


%=== initialize output structure
info = struct('Filename', '', 'FilePath', '', 'FileModDate', '', ...
              'FileSize', '', 'NumberOfEntries', []);
          
%=== check input
if ~ischar(filename) && ~iscellstr(filename)
    error(message('bioinfo:fastqinfo:InvalidInput'))
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
        error(message('bioinfo:fastqinfo:InvalidInput'))
    end
end

% Check is the string has a '\n' in the first 2^12 bytes, if so, it can
% not be a file name
singleLine = isempty(regexp(filename(1:min(end,2^12)),'\n','once'));

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
    
    numLines = 0;
    blockSize = 2^24;
    count = blockSize;
    warnFormat = false;
    determineFormat = true;
    emptyQualityHeader = false;
    sh = 0;
    while count == blockSize
        c1 = textscan(fid,'%c%*s',blockSize,'delimiter','\n');
        % empty lines at the end are removed
        nl = regexp(c1{1}','\n','once');
        if ~isempty(nl)
            c1{1} = c1{1}(1:nl-1);
        end
        count = numel(c1{1});
        numLines = numLines + count;
        if determineFormat
            if all(c1{1}(1:4:end) == '@')
                emptyQualityHeader = false;
                determineFormat = false;
            elseif all(c1{1}(1:3:end) == '@')
                emptyQualityHeader = true;
                determineFormat = false;
            else
                warnFormat = true;
                determineFormat = false;
            end
        end
           
        if ~warnFormat
          if emptyQualityHeader
            if ~(all(c1{1}(sh+1:3:end) == '@') && all(c1{1}(sh+3:3:end) == '+'))
              warnFormat = true;
            end
            sh = rem(sh+2,3);
          else
            % no interleave empty lines are expected so the current block
            % should always start with a line with a '@'
            if ~(all(c1{1}(1:4:end) == '@') && all(c1{1}(3:4:end) == '+'))
              warnFormat = true;
            end
          end
        end
    end
else % must be a string with '\n' (UNDOCUMENTED)
    dataFromFile = false;
    c1 = textscan(filename,'%c%*s','delimiter','\n');
    % empty lines at the end are removed
    nl = regexp(c1{1}','\n','once');
    if ~isempty(nl)
        c1{1} = c1{1}(1:nl-1);
    end    
    numLines = numel(c1{1});
    % Check format
    if all(c1{1}(1:4:end) == '@') && all(c1{1}(3:4:end) == '+')
        emptyQualityHeader = false;
        warnFormat = false;
    elseif all(c1{1}(1:3:end) == '@') && all(c1{1}(3:3:end) == '+')
        emptyQualityHeader = true;
        warnFormat = false;
    else
        emptyQualityHeader = false;
        warnFormat = true;
    end
end

if ~dataFromFile && numLines == 1 && all(filename~=sprintf('\t')) && all(filename~=sprintf('\n'))
    % Identify if the string in filename was an intended filename but
    % since it was not initially found the function assumed it was a
    % a string with '\n' and counted the lines anyways, if we counted
    % one single line we assume the user intended to pass a filename
    % and we give a more appropriate error:
    error(message('bioinfo:fastqinfo:InvalidFileName'));
end

if warnFormat || ((mod(numLines,4)~=0)&&~emptyQualityHeader) ||...
                 ((mod(numLines,3)~=0)&&emptyQualityHeader)
    warning(message('bioinfo:fastqinfo:InvalidFormat'))
end

if emptyQualityHeader
   info.NumberOfEntries = uint64(floor(numLines/3));
else
   info.NumberOfEntries = uint64(floor(numLines/4));
end


