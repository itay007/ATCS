function info = fastainfo(filename)
%FASTAINFO return a summary of the contents of a FASTA file.
%
%   INFO = FASTAINFO(FILENAME) returns a structure whose fields contain
%   information about a FASTA file.  FILENAME is a string containing a file
%   name, or a path and a file name, of a FASTA file.  FILENAME can also be
%   a URL or MATLAB character array that contains the text of a FASTA
%   format file. INFO is a structure with the following fields:
%
%               Filename - name of the file
%               FilePath - path to the file
%            FileModDate - modification date of the file
%               FileSize - size of the file in bytes
%        NumberOfEntries - number of sequence entries
%
%   When the FASTA file contains one sequence, FASTAINFO also returns the
%   contents of the header and the sequence length in the output structure.
%
%   Example:
%   % Show a summary of the content of a FASTA file.
%   info = fastainfo('p53nt.txt')
%
%   See also FASTAREAD, FASTAWRITE, FASTQINFO, FASTQREAD, FASTQWRITE,
%   SFFINFO, SFFREAD.

%   Copyright 2009-2010 The MathWorks, Inc.



%=== initialize output structure
info = struct('Filename', '', 'FilePath', '' , 'FileModDate', '',...
              'FileSize', '', 'NumberOfEntries', []);
          
%=== check input
if ~ischar(filename) && ~iscellstr(filename)
    error(message('bioinfo:fastainfo:InvalidInput'))
end

%=== Padded strings are deprecated, but we still accept them
if size(filename,1) > 1
    if isvector(filename)
        filename = filename(:)';
    else
        filename = cellstr(filename);
    end
end

%=== cellstr's are valid (UNDOCUMENTED), it is changed to a string with '\n'
%    to further work with it:
if iscellstr(filename)
    filename = sprintf('%s\n',filename{:});
end

%=== Check if we have an URL, is so we fetch it
if (strfind(filename(1:min(10,end)), '://'))
    if (~usejava('jvm'))
        error(message('bioinfo:fastainfo:NoJava'))
    end
    
    try
        urltext = urlread(filename);
    catch urlErr
        error(message('bioinfo:fastainfo:CannotReadURL', filename));
    end
    filename = urltext;
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
    blockSize = 2^20;
    numLinesRead = blockSize;
    while numLinesRead == blockSize
        c1 = textscan(fid,'%s',blockSize,'delimiter','\n');
        numLinesRead = numel(c1{1});
        c1 = strtrim(c1{1});
        count = sum(~cellfun(@isempty,regexp(c1(~cellfun(@isempty,c1)),'^>')));
        numLines = numLines + count;
    end
    % For FASTA files with only one sequence we also fetch the header and
    % calculate the length:
    if numLines == 1
        mmff =  bioinfoprivate.MemoryMappedFastaFile(filename);
        info.Header = mmff.Header;
        info.Length = uint64(mmff.Length);
    end
else % must be a string with '\n'
    dataFromFile = false;
    c1 = textscan(filename,'%s','delimiter','\n');
    c1 = strtrim(c1{1});
    c1 = c1(~cellfun(@isempty,c1));
    h = ~cellfun(@isempty,regexp(c1,'^>'));
    numLines = sum(h);
    % For FASTA files with only one sequence we also fetch the header and
    % calculate the length:
    if numLines == 1
        info.Header = regexprep(c1{h},'>','');
        info.Length = uint64(sum(cellfun(@length,c1(~h))));
    end
end

if ~dataFromFile && numLines == 0 && all(filename~=sprintf('\t')) && all(filename~=sprintf('\n'))
    % Identify if the string in filename was an intended filename but
    % since it was not initially found the function assumed it was a
    % a string with '\n' and counted the lines anyways, if we counted
    % one single line we assume the user intended to pass a filename
    % and we give a more appropriate error:
    error(message('bioinfo:fastainfo:InvalidFileName'));
end

if numLines == 0
        error(message('bioinfo:fastainfo:InvalidOrNotFoundInput'));
end

info.NumberOfEntries = uint64(numLines);
    
