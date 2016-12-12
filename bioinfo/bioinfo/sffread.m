function out = sffread(filename, varargin)
% SFFREAD Read Standard Flowgram Format (SFF) files.
%
%   OUT = SFFREAD(FILENAME) reads an SFF format file FILENAME, returning
%   the data in the file as a structure array OUT. FILENAME is a string
%   containing a file name, or a path and a file name, of a SFF file that
%   conforms to the format specifications in version 1. OUT is a MATLAB
%   structure array with the following fields: 
%           Header - universal accession number 
%         Sequence - numeric representation of nucleotide sequence
%          Quality - quality scores of bases
% 
%   Additional fields are selected using the FEATURE option. These include
%   the following fields:
%          
%         Clipping - clipping boundary positions
%    FlowgramValue - sequence of flowgram intensity values
%    FlowgramIndex - sequence of flowgram intensity indices
%
%   SFFREAD(..., 'BLOCKREAD', M) allows you to read in a single entry or
%   block of entries from a file. If M is a scalar then the Mth entry in
%   the file is read. If M is a two element vector then the block of
%   entries starting at entry M(1) and ending at entry M(2) will be read.
%   Use Inf for M(2) to read all entries in the file starting at position
%   M(1).
%
%   SFFREAD(..., 'FEATURE', FEATURESEL) allows you to specify what
%   information to include in the output structure. FEATURESEL is a string
%   from the alphabet {h,s,q,c,f,i} representing header, sequence, quality,
%   clipping, flowgram values and flowgram indices respectively. Any
%   combination of the valid characters is allowed. For example, 'sq'
%   specifies sequence and quality as fields in the output structure. By
%   default FEATURESEL is set to 'hsq'.
%
%   Examples:
%   % Read the contents of an SFF file into a MATLAB structure.
%   reads = sffread('SRR013472.sff')
%
% 	% Read the Header, Sequence and Clipping information of the first 5 reads. 
% 	reads5 = sffread('SRR013472.sff', 'block', [1 5], 'feature', 'hsc');
%
%   See also FASTAINFO, FASTAREAD, FASTAWRITE, FASTQINFO, FASTQREAD,
%   FASTQWRITE, SFFINFO.

%   Copyright 2009-2010 The MathWorks, Inc.




bioinfochecknargin(nargin,1,mfilename)
[headFlag, seqFlag, qualFlag, clipFlag, flowFlag, indexFlag, block] = parse_inputs(varargin{:});
	
if ~isempty(block)
	blockOnly = true;
else
	blockOnly = false;
end


%=== check input
if ~ischar(filename)
	error(message('bioinfo:sffread:InvalidInput'))
end

if (~exist(filename,'file') && ~exist(fullfile(pwd,filename),'file'))
	error(message('bioinfo:sffread:FileNotFound', filename))
end

%=========
% main 
%=========

try
	%=== check endianess
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
			error(message('bioinfo:sffread:InvalidMagicNumber'));
		end
	end
	
	%======================================================================
	% process common header section
	%======================================================================
	
	%=== confirm version
	versionNumber = fread(fid, 4, 'char')';
	if ~isequal(versionNumber, [0 0 0 1])
		warning(message('bioinfo:sffread:InvalidFileVersion'));
	end
	
	%=== extract info from common header section
	indexOffset = fread(fid,1,'uint64'); 
	indexLength = fread(fid,1,'uint32');  

	numReads = fread(fid,1,'uint32'); 
	headerLength = fread(fid,1,'uint16'); 
	keyLength = fread(fid,1,'uint16'); 
	numFlowsPerRead = fread(fid,1,'uint16');
	flowgramFormatCode = fread(fid,1,'uint8');
	flowgramBytesPerFlow = 2 * flowgramFormatCode; % 2 is the value assigned to format code = 1;
		
	%=== check header
	headerLengthCheck = ceil((31 + numFlowsPerRead + keyLength)/8) * 8;
	if headerLengthCheck ~= headerLength
		error(message('bioinfo:sffread:InvalidHeaderLength', filename))
	end
	
	% NOTE: using fread instead of fseek for performance's sake
	
	%=== skip flow chars
	fread(fid, numFlowsPerRead, '*char'); % flowChars
	fread(fid, keyLength, '*char');  % keySequence
	%fseek(fid, numFlowsPerRead+keyLength, 0);

	%=== skip the eight_byte_padding (see 454 manual)
	eightBytePadding = headerLength - ftell(fid);
	fread(fid, eightBytePadding);
	%fseek(fid, eightBytePadding, 0);
		
	%=== check what sequence block needs to be extracted, if any
	if blockOnly

		if numel(block) == 1
			block(2) = block(1);
		elseif block(2) == Inf
			block(2) = numReads;
		end
				
		if (block(1) <= 0) || (block(1) > numReads)
			error(message('bioinfo:sffread:BadBlockLowBoundary'));
		elseif (block(2) <= 0) || (block(2) > numReads)
			error(message('bioinfo:sffread:BadBlockUpBoundary'));
		end
		numOutputReads = diff(block) + 1;
	else
		block = [1 numReads];
		numOutputReads = numReads;
	end
	
	%=== initialize output structure
	str = struct;
	
	if headFlag
		str.Header = '';
	end
	if seqFlag
		str.Sequence = '';
	end
	if qualFlag 
		str.Quality = [];
	end
	if clipFlag 
		str.Clipping = [];
	end
	if flowFlag
		str.FlowgramValue = [];
	end
	if indexFlag
		str.FlowgramIndex = [];
	end
	
	out = repmat(str, numOutputReads, 1);
	
	%======================================================================
	% process reads
	%======================================================================
	i = 0; % index of current read
	b = 0; % index of current output read
	
	while (b < numOutputReads)
				
		%==================================================================
		% process read header section
		%==================================================================
		i = i+1;
		
		%=== skip index if present
		if indexOffset == ftell(fid)
			fread(fid, indexLength, 'uint8'); % skip index bytes
		end
					
		readHeaderLength = fread(fid,1,'uint16');
		nameLength = fread(fid,1,'uint16');
		numBases = fread(fid,1,'uint32');
		readDataLength = numFlowsPerRead * flowgramBytesPerFlow + 3 * numBases;
		readDataLengthRounded = double(ceil(readDataLength/8) * 8); % rounded to next number divisble by 8
		
		if i >= block(1) && i <= block(2) % keep the current read info
			
			b = b+1;
			
			if clipFlag 
				out(b).Clipping = fread(fid,4,'*uint16'); % qual left, qual right, adapter left, adapter right
			else
				%fseek(fid,8,0); % skip clipping info
				fread(fid,4,'*uint16');
			end
			
			if headFlag
				out(b).Header = fread(fid,nameLength, '*char')';
			else
				%fseek(fid, nameLength, 0); % skip header info
                fread(fid,nameLength, '*char');
			end
			
			eightBytePadding = readHeaderLength - 16 - nameLength;
			%fseek(fid, eightBytePadding, 0);
			fread(fid, eightBytePadding);
			
			%==============================================================
			% process read data section
			%==============================================================
			
			if flowFlag
				out(b).FlowgramValue = fread(fid, numFlowsPerRead, '*uint16')'; % depends on flowgram format code
			else
				%fseek(fid, numFlowsPerRead * 2, 0); % skip the flowgram portion
				fread(fid, numFlowsPerRead * 2);
			end
			
			if indexFlag
				out(b).FlowgramIndex = fread(fid, numBases, '*uint8')';
			else
				%fseek(fid, numBases, 0); % skip the flowgram index portion
				fread(fid, numBases);
			end
			
			if seqFlag
				out(b).Sequence = fread(fid, numBases, '*uint8')';
			else
				%fseek(fid, numBases, 0); % skip the sequence portion
				fread(fid, numBases); % skip the sequence portion
			end
			
			if qualFlag
				out(b).Quality = fread(fid, numBases, '*uint8')';
			else
				%fseek(fid, numBases, 0); % skip the quality portion
				fread(fid, numBases);
			end
			

			eightBytePadding = readDataLengthRounded - readDataLength;
			%fseek(fid, eightBytePadding, 0);
			fread(fid, eightBytePadding);
						
		else 
			%=== skip the current read (header and data)
			%fseek(fid, readHeaderLength - 8 + readDataLengthRounded, 0);
			fread(fid, readHeaderLength - 8 + readDataLengthRounded);
		end
	end

catch theErr
	if strfind(theErr.identifier, mfilename)
		rethrow(theErr)
		
	elseif strcmpi(theErr.identifier,'MATLAB:nomem')
            error(message('bioinfo:sffread:FileTooBig'));

	else
		error(message('bioinfo:sffread:CannotReadFile', filename));
	end
end
	

%--------------------------------------------------------------------------
function [headFlag, seqFlag, qualFlag, clipFlag, flowFlag, indexFlag, block] =  parse_inputs(varargin)
% Parse input PV pairs.

%=== set default values
headFlag = true;   % return the universal accession numbers (headers)
seqFlag = true;    % return sequence  
qualFlag = true;   % return quality scores
clipFlag = false;  % return clipping info
flowFlag = false;  % return flowgram intensity values
indexFlag = false; % return flowgram indices
block = []; % boundaries of the block to return (if any)

%=== parse parameters
if  nargin > 1
	if rem(nargin,2) == 1
		error(message('bioinfo:sffread:IncorrectNumberOfArguments', mfilename));
	end
	okargs = {'feature','blockread'};
	for j = 1:2:nargin-1
		pname = varargin{j};
		pval = varargin{j+1};
		k = find(strncmpi(pname,okargs,numel(pname)));
		if isempty(k)
			error(message('bioinfo:sffread:UnknownParameterName', pname));
		else
			switch(k)
				case 1  % feature
					if ischar(pval)
						m = regexpi(pval, '[hsqcfi]');
						if length(m) ~= length(pval)
							error(message('bioinfo:sffread:UnknownFeatureOption'));
						else
							pval = lower(pval);
							headFlag = any(strfind(pval, 'h'));
							seqFlag = any(strfind(pval, 's'));
							qualFlag = any(strfind(pval, 'q'));
							clipFlag = any(strfind(pval, 'c'));
							flowFlag = any(strfind(pval, 'f'));
							indexFlag = any(strfind(pval, 'i'));
						end
					else
                        error(message('bioinfo:sffread:InvalidFeatureOption'));
					end
					
				case 2  % blockread
					block = pval;
					if ~isnumeric(block) || numel(block)> 2 || isempty(block)
						error(message('bioinfo:sffread:BadBlock'))
					end
					block = sort(block);
			end
		end
	end
end
