function [data, seq, qual] = fastqread(filename, varargin)
%FASTQREAD reads FASTQ file.
%
%   OUT = FASTQREAD(FILENAME) reads a FASTQ file FILENAME,
%   returning the data stored in the file as an array of structures OUT
%   with the following fields:
%           Header - header information
%         Sequence - sequence information
%          Quality - ASCII representation of quality score
%
%   [HEADER, SEQ] = FASTQREAD(FILENAME) reads the file into separate
%   variables HEADER and SEQ. If the file contains more than one sequence,
%   then HEADER and SEQ are cell arrays of header and sequence information.
%
%   [HEADER, SEQ, QUAL] = FASTQREAD(FILENAME) reads the file into separate
%   variables HEADER, SEQ and QUAL. If the file contains more than one
%   sequence, then HEADER, SEQ and QUAL are cell arrays of header, sequence
%   and quality information.
%
%   FASTQREAD(...,'BLOCKREAD', M) allows you to read a single entry or
%   block of entries from a file containing multiple sequences. If M is a
%   scalar then the Mth entry in the file is read. If M is a two element
%   vector then the block of entries starting at entry M(1) and ending at
%   entry M(2) are read.  Use Inf for M(2) to read all entries in the
%   file starting at position M(1).
%
%   FASTQREAD(...,'HEADERONLY', TF) allows you to specify whether to
%   return only the headers or not. Default is false.
%
%   FASTQREAD(...,'TRIMHEADERS',TF) trims the header after the first
%   whitespace when TF is true. White space characters include a horizontal
%   tab (char(9)) and a space (char(32)). Default is false.
%
%   Examples:
%       % Read the contents of a FASTQ file into a MATLAB structure.
%       reads = fastqread('SRR005164_1_50.fastq')
%
% 		% Read the contents of a FASTQ file into separate variables.
% 		[h,s,q] = fastqread('SRR005164_1_50.fastq');
%
% 		% Read a block of entries from a file.
% 		reads_5_10 = fastqread('SRR005164_1_50.fastq', 'blockread', [5 10])
%
%   See also FASTAINFO, FASTAREAD, FASTAWRITE, FASTQINFO, FASTQWRITE,
%   SFFINFO, SFFREAD.

%   Copyright 2009-2012 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename)
[range, honly, trimh] = parse_inputs(varargin{:});

if isempty(range)
	blockRead = false;
else
	blockRead = true;
end

%=== check input
if ~ischar(filename)
	error(message('bioinfo:fastqread:InvalidInput'))
end

%=== read the content
try
	%=== read from a padded string (UNDOCUMENTED)
    if size(filename,1) > 1
        dataFromFile = false;
        if blockRead
            warning(message('bioinfo:fastqread:IgnoredRange'))
        end
        ftext = textscan(filename, '%s', 'delimiter', '\n');
        ftext = ftext{:};
        
        %=== read from a file
    elseif (exist(filename,'file') || exist(fullfile(pwd,filename),'file'))
        dataFromFile = true;
        fid = fopen(filename);
        c = onCleanup(@()fclose(fid));
        
        if blockRead
            
            if range(2) == Inf % read entire file after skipping beginning
                ftext = textscan(fid, '%s', 'delimiter', '\n', ...
                    'headerlines', (range(1)-1)*4);
                
            else % read only block
                numItems = diff(range) + 1; % number of entries to consider
                ftext = textscan(fid, '%s', numItems * 4, ...
                    'delimiter', '\n', 'headerlines', (range(1)-1)*4);
            end
            
            ftext = ftext{:};
            
            if isempty(ftext)
                error(message('bioinfo:fastqread:BadBlockRange', filename));
            end
            %fclose(fid);
            
        else %=== read entire file
            ftext = textscan(fid, '%s', 'delimiter', '\n');
            ftext = ftext{:};
        end
        
        %=== error (it is not a string nor file)
    else
        dataFromFile = false;
        if blockRead
            warning(message('bioinfo:fastqread:IgnoredRange'))
        end
        ftext = strread(filename,'%s','delimiter','\n');
    end
		

	%=== check integrity of content
	N = numel(ftext);
	
	if mod(N,4) ~= 0
		error(message('bioinfo:fastqread:BadFastqFile', filename))
	end
	
	%=== extract info
	if honly
		data = repmat(struct('Header', []), 1, N/4);
	else
		data = repmat(struct('Header', [], 'Sequence', [], 'Quality',[]), 1, N/4);
		[data(1:floor(N/4)).Sequence] = ftext{2:4:N};
		[data(1:floor(N/4)).Quality] = ftext{4:4:N};
	end
	[data(1:floor(N/4)).Header] = ftext{1:4:N};
		
	%=== remove tag @ from header
	h = regexprep({data(:).Header}, '^@', '');
	[data(:).Header] = h{:};

catch theErr
	if strcmpi(theErr.identifier,'MATLAB:nomem')
		error(message('bioinfo:fastqread:FileTooBig'));
    elseif strncmpi('bioinfo:',theErr.identifier,8)
        % Identify when the actual string in filename was an intended
        % filename but since it was not initially found the function tried
        % to parse it.
        if dataFromFile || (size(filename,1)==1 && any(filename==10))
           rethrow(theErr)
        else
           error(message('bioinfo:fastqread:InvalidFileName'));
        end
    else % All other non-bioinfo errors are presumed due to an
         % incorrect format of the source.
        error(message('bioinfo:fastqread:InvalidOrNotFoundInput'));
	end
end

% trim headers
if trimh
   for i = 1:numel(data)
      data(i).Header = sscanf(data(i).Header,'%s',1);
   end
end
    
%=== in case of two or more outputs
if nargout > 1 
	
	if honly % assign output variables to empty
		if nargout == 3
			qual = [];
		end
		
		seq = '';
		
		if N == 4 % one read only
			data = data.Header;
		else
			data = {data(:).Header};
		end
		
	else
		if nargout == 3
			if N == 4 % one sequence read only
				qual = data.Quality;
			else
				qual = {data(:).Quality};
			end
		end
		
		if N == 4 % one read only
			seq = data.Sequence;
			data = data.Header;
		else
			seq = {data(:).Sequence};
			data = {data(:).Header};
		end
	end
	
elseif honly
	data = {data(:).Header};
end

%--------------------------------------------------------------------------
function [range, honly, trimh] = parse_inputs(varargin)
% Parse input PV pairs.

%=== Check for the right number of inputs
if rem(nargin,2) == 1
	error(message('bioinfo:fastqread:IncorrectNumberOfArguments', mfilename));
end

%=== Defaults
range = [];    % range of entries to return
honly = false; % headers only
trimh = false; % trim headers

%=== Allowed inputs
okargs = {'blockread', 'headeronly', 'trimheaders'};

for j = 1:2:nargin-1
	
	[k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
	
	switch(k)
		case 1  % blockread
			range = pval;
			if ~isnumeric(range) || numel(range)> 2 || isempty(range)
				error(message('bioinfo:fastqread:BadBlock'));
			end
			
			if numel(range) == 1
				range(2) = range(1);
			end
					
			range = sort(double(range));
			
		case 2 % headeronly
			honly = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
            
        case 3 % trimheaders
			trimh = bioinfoprivate.opttf(pval, okargs{k}, mfilename);    
	end
end
