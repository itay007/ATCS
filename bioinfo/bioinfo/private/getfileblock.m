function block = getfileblock(filename,range,delimiter)
% GETFILEBLOCK reads a block of a file
%
%   GETFILEBLOCK(FILENAME,RANGE,DELIMITER) reads a block of text from file
%   FILENAME. The function assumes that the file is made up of entries with
%   a delimiter DELIMITER marking the start or end of an entry. The block
%   will start at RANGE(1) and stop at RANGE(2).
%

%   Copyright 2002-2012 The MathWorks, Inc.

% set start and stop and make sure they are reasonable
start = floor(range(1));
if start < 1
    start = 1;
end
if numel(range) > 1
    stop = ceil(range(2));
else
    stop = start;
end

% FASTAREAD allows to open files from the path as it validates FILENAME
% using exist(), however this is not allowed by BioIndexedFile. We need to
% pre-pend the right path to a file that is expected to be found in the
% MATLAB's path:
if isempty(fileparts(filename))
    filename = which(filename);
end

% index the fasta file
tempindexfile = [tempname '.idx'];
c = onCleanup(@()delete(tempindexfile));
i = BioIndexedFile('flat',filename,tempindexfile,...
            'indexedbykeys',false,...
            'memorymappedindex',false,...
            'verbose',false,...
            'entrydelimiter',delimiter);

% test if query goes over range        
if start>i.NumEntries
    error(message('bioinfo:getfileblock:StartTooBig',i.NumEntries));    
end

try
    block = '';
    for j = start:min(stop,i.NumEntries)
        block = sprintf( '%s%s%s\n',block,delimiter,getEntryByIndex(i,j));
    end
catch theErr
    if strcmpi(theErr.identifier,'MATLAB:nomem')
        error(message('bioinfo:getfileblock:BlockTooBig'));
    else
        rethrow(theErr);
    end
end
