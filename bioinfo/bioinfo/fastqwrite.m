function fastqwrite(filename, data, dataseq, dataqual)
%FASTQWRITE writes FASTQ file.
%
%   FASTQWRITE(FILENAME, DATA) writes the contents of DATA to file
%   FILENAME in the FASTQ format. DATA is a structure containing the
%   fields Header, Sequence, and Quality.
%
%   FASTQWRITE(FILENAME, HEADER, SEQ, QUAL) writes a FASTQ file with header
%   information HEADER, sequence information SEQ and quality information QUAL.
%
%
%   Examples:
%       % Read the contents of a FASTQ file and save the first 5 entries 
%       % into a different file.
%       reads = fastqread('SRR005164_1_50.fastq');
%       reads5 = reads(1:5);
%       fastqwrite('fiveReads.fastq', reads5);
%
% 		% Write a FASTQ file from separate variables.
%       h = 'MYSEQ-000_1_1_1_953_493';
%       s = 'GTTACCATGATGTTATTTCTTCATTTGGAGGTAAAA';
%       q = ']]]]]]]]]]]]]]]]]]]]]]T]]]]RJRZTQLOA';
%       fastqwrite('oneRead.fastq', h, s, q);
%
%
%   See also FASTAREAD, FASTAINFO, FASTAWRITE, FASTQINFO, FASTQREAD,
%   SFFINFO, SFFREAD.

%   Copyright 2009-2010 The MathWorks, Inc.



bioinfochecknargin(nargin,2,mfilename)
appendFile = false;

%=== check inputs
if ~ischar(filename),
	error(message('bioinfo:fastqwrite:FilenameMustBeString'));
end

%=== check whether to create a new file or appending to existing file
if exist(filename, 'file')
	appendFile = true;
end

%=== open file
fid = fopen(filename, 'a');

if fid == (-1)
	[theDir, theFile, theExtension] = fileparts(filename);
	if ~isempty(theDir)
		error(message('bioinfo:fastqwrite:CouldNotOpenFileinDir', [ theFile, theExtension ], theDir));
	else
		error(message('bioinfo:fastqwrite:CouldNotOpenFileinPwd', filename));
	end
end

%=== warn if file already exists
if appendFile
	warning(message('bioinfo:fastqwrite:AppendToFile', filename));
end

%=== main
try
	if nargin == 2 % we must have a structure
		if ~isstruct(data)
			error(message('bioinfo:fastqwrite:DataMustBeStruct'));
		end
		
		fn = fieldnames(orderfields(data));
		if ~isequal(fn', {'Header','Quality','Sequence'})
			error(message('bioinfo:fastqwrite:DataMustHaveFields'));
		end
		
	elseif nargin == 4 % we must have header, sequence and quality
		datahead = data;
		
			data = struct('Header', '', 'Sequence', '', 'Quality', '');
			if iscell(datahead) 
				if numel(datahead) == numel(dataseq) && numel(dataseq) == numel(dataqual)
					N = numel(datahead);
					data = repmat(data, 1, N);
					[data.Header] = datahead{:}; 
					[data.Sequence] = dataseq{:};
					[data.Quality] = dataqual{:};
				else
					error(message('bioinfo:fastqwrite:IncorrectInputSize'))
				end
			else
				data.Header = datahead;
				data.Sequence = dataseq;
				data.Quality = dataqual;
			end
		
	else
		error(message('bioinfo:fastqwrite:IncorrectNumberOfInputs'));
	end
	
	%=== write data
	for i = 1:numel(data)
		currSeq = data(i).Sequence;
		currHeader = data(i).Header;
		currQual = data(i).Quality;
		
		lenSeq = length(currSeq);
		lenQual = length(currQual);
		
		if lenSeq == 0 || lenQual == 0 
			error(message('bioinfo:fastqwrite:EmptyEntry'));
		end
		
		if lenSeq ~= lenQual
			error(message('bioinfo:fastqwrite:SeqQualUnequalSize'));
		end
		%=== print sequence header
		fprintf(fid,'@%s\n',currHeader);
		
		%=== print sequence
		fprintf(fid,'%s\n', currSeq);
		
		%=== print quality header
		fprintf(fid,'+%s\n',currHeader);
		
		%=== print quality
		fprintf(fid,'%s\n',currQual);
		
	end
	fclose(fid);
	
catch theErr
	
	fclose(fid);
	
	if ~appendFile
		%=== delete the file if it was not already an existing file
		delete(filename);
	end
		
	if ~isempty(strfind(theErr.identifier, mfilename)) 
		rethrow(theErr);
	else
		error(message('bioinfo:fastqwrite:CannotWriteFile', filename))
	end
end
