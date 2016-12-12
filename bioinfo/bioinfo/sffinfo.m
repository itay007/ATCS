function info = sffinfo(filename)
%SFFINFO return a summary of the contents of an SFF file.
%
%   INFO = SFFINFO(FILENAME) returns a structure whose fields contain
%   information about an SFF file.  FILENAME is a string containing a file
%   name, or a path and a file name, of an SFF file.  INFO is a structure
%   with the following fields:
%
%              Filename - name of the file
%           FileModDate - modification date of the file
%              FileSize - size of the file in bytes
%               Version - version number of the file
%          FlowgramCode - code of the format used to encode flowgram values
%         NumberOfReads - number of reads stored in the file
%  NumberOfFlowsPerRead - number of flow for each read
%             FlowChars - bases used in each flow
%           KeySequence - string of bases in the key sequence
%   
%   Example:
%   % Show information relative to the contents of a local SFF file.
%   info = sffinfo('SRR013472.sff')
%
%   See also FASTAINFO, FASTAREAD, FASTAWRITE, FASTQINFO, FASTQREAD,
%   FASTQWRITE, SFFREAD.

%   Copyright 2009-2010 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename)

%=== check input file
if (~exist(filename,'file')) && (~exist(fullfile(pwd,filename),'file'))
	error(message('bioinfo:sffinfo:FileNotFound', filename))
end

%=== check if filename contains full path, if not get it
if isempty(regexp(filename, filesep, 'once'))
    filename = which(filename);
end

%=== check endianess
try
	fid = fopen(filename);
	c = onCleanup(@()fclose(fid));
	
	magicNumber = dec2hex(fread(fid, 1, 'uint32'));
	if ~strcmp(magicNumber, '2E736666') % wrong endianess
		fclose(fid);
		[~, ~, endian] = computer;
		
		%=== reopen with correct endianess
		if endian == 'L'
			fid = fopen(filename, 'r', 'b');
		else
			fid = fopen(filename, 'r', 'l');
		end
		
		magicNumber = dec2hex(fread(fid, 1, 'uint32'));
		if ~strcmp(magicNumber, '2E736666')
			error(message('bioinfo:sffinfo:InvalidMagicNumber'));
		end
	end
	
	%=== confirm version
	versionNumber = fread(fid, 4, 'char')';
	if ~isequal(versionNumber, [0 0 0 1])
		warning(message('bioinfo:sffinfo:InvalidFileVersion'));
	end
	
	%=== skip information relative to the index
	fread(fid,1,'uint64'); 
	fread(fid,1,'uint32'); 
		
	numReads = fread(fid,1,'uint32'); 
	headerLength = fread(fid,1,'uint16'); 
	keyLength = fread(fid,1,'uint16'); 
	numFlowsPerRead = fread(fid,1,'uint16');
	flowgramFormatCode = fread(fid,1,'uint8');
	flowgramBytesPerFlow = 2 * flowgramFormatCode; % 2 is the value assigned to format code = 1;
	flowChars = fread(fid, numFlowsPerRead, '*char')';
	keySequence = fread(fid, keyLength, '*char')';

	fileInfo = dir(filename);
	info.Filename = fileInfo.name;
	info.FileModDate = fileInfo.date;
	info.FileSize = fileInfo.bytes;
	%info.MagicNumber = magicNumber;
	info.Version = versionNumber;
	info.FlowgramCode = flowgramFormatCode;
	info.NumberOfReads = numReads;
	info.NumberOfFlowsPerRead = numFlowsPerRead;
	info.FlowChars = flowChars;
	info.KeySequence = keySequence;

	
catch theErr
	
	if strfind(theErr.identifier,'sffinfo')
		rethrow(theErr)
	elseif strcmpi(theErr.identifier,'MATLAB:nomem')
            error(message('bioinfo:sffinfo:FileTooBig'));
	else
		error(message('bioinfo:sffinfo:CannotReadFile', filename));
	end
end
