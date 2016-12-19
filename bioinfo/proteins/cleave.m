function [parts, cuttingSites, lengths, missed] = cleave(seq, pattern, varargin)
%CLEAVE cuts an amino acid sequence at specified cleavage sites.
%
%   FRAGMENTS = CLEAVE(SEQ, ENZYME) cuts sequence SEQ into parts at
%   cleavage sites specific for enzyme ENZYME. ENZYME is the name or
%   abbreviation code of an enzyme or compound for which a cleavage rule
%   has been specified in the literature. Use CLEAVELOOKUP to display the
%   names of enzymes and compounds in the current cleavage rule library.
%
%   FRAGMENTS = CLEAVE(SEQ, PATTERN, POSITION) cuts sequence SEQ into parts
%   at cleavage sites specified by peptide pattern PATTERN.
%   POSITION defines the position on the PATTERN where the sequence is cut.
%   POSITION  0 corresponds to the N terminal end of the PATTERN. PATTERN
%   can be a regular expression.
%
%   CLEAVE(..., 'PARTIALDIGEST', P) simulates a partial digest where each
%   cleavage site in the sequence has a probability P of being cut.
%
%   CLEAVE(..., 'MISSEDSITES', M) returns all possible peptides that
%   can result from missing M or less cleavage sites. M must be a non-
%   negative integer. M defaults to 0, which is equivalent to an ideal
%   digestion.
%
%   CLEAVE(..., 'EXCEPTION', rule) specifies an exception rule (regular
%   expression) to the cleavage rule associated with the given enzyme. By
%   default no exception rule is applied, i.e. the exception rule is an
%   empty string.
%
%   [FRAGMENTS, CUTTING_SITES] = CLEAVE(...) returns a numeric vector
%   with the indices representing the cleave sites. A 0 (zero) will be
%   added to the list, so numel(FRAGMENTS) == numel(CUTTING_SITES).
%   CUTTING_SITES+1 can be used to point to the first amino acid of every
%   fragment respective to the original sequence.
%
%   [FRAGMENTS, CUTTING_SITES, LENGTHS] = CLEAVE(...) returns a numeric
%   vector with the lengths of every fragment.
%
%   [FRAGMENTS, CUTTING_SITES, LENGTHS, MISSED] = CLEAVE(...) returns a
%   numeric vector with the number of missed cleavage sites of every
%   fragment.
%
%
%   Some common proteases and their cleavage sites are as follows:
%
%       Trypsin:         [KR](?!P)     Position 1
%       Chymotrypsin:    [WYF](?!P)    Position 1
%       Glu C:           [ED](?!P)     Position 1
%       Lys C:           [K](?!P)      Position 1
%       Asp N:           D             Position 1
%
%   Examples:
%
%       S = getgenpept('AAA59174')
%       % Cut sequence using  Asp-N endopeptidase
%       cleave(S.Sequence, 'Asp-N')
%
%       % Cut sequence using Trypsin, which cleaves after K or R when the next residue is not P.
%       [parts, sites, lengths] = cleave(S.Sequence,'[KR](?!P)',1);
%       for i = 1:10
%           fprintf('%5d%5d   %s\n',sites(i),lengths(i),parts{i})
%       end
%
%       % Cut sequence using Trypsin and allow for 1 missed cleavage site
%       [parts, sites, lengths, missed] = cleave(S.Sequence,'trypsin','missedsites',1);
%
%
%       % Cut sequence using Trypsin's cleavage rule with exception
%       parts = cleave(S.Sequence,'trypsin','exception', 'KD');

%   See also CLEAVELOOKUP, REGEXP, RESTRICT, REBASECUTS, SEQSHOWWORDS.

%   Reference: Liebler, D. Introduction to Proteomics, Humana Press 2002.

%   Copyright 2003-2010 The MathWorks, Inc.


%=== check inputs and outputs
bioinfochecknargin(nargin,2,mfilename);

returnCuttingSites = false;
returnLengths = false;
returnMissed = false;

if nargout > 1
	returnCuttingSites = true;
end
if nargout > 2
	returnLengths = true;
end
if nargout > 3
	returnMissed = true;
end

%=== Extract sequence data from structure
if isstruct(seq)
	seq = bioinfoprivate.seqfromstruct(seq);
end

cut2 = [];

%=== retrieve pattern and position
if nargin  == 2 % must be (seq, enzyme)
	[position, pattern] = localLookupEnzyme(pattern);
elseif nargin > 2
	if isnumeric(varargin{1}) % must be (seq, patt, pos, [...])
		position = varargin{1};
		varargin = varargin(2:end);
		pattern = seq2regexp(pattern,'alphabet','aa');
	else % must be (seq, enzyme, varargins)
		[position, pattern] = localLookupEnzyme(pattern);
	end
end

nvarargin = numel(varargin);

if rem(nvarargin,2) == 1 || nargin > 7
	error(message('bioinfo:cleave:IncorrectNumberOfArguments', mfilename));    
end
 
%=== parse parameter/value pairs
[ratio, M, exceptionRule] = parse_inputs(varargin{:});

%=== validate position
if ~isnumeric(position)
	error(message('bioinfo:cleave:NonNumericPosition'))
end

%=== deal with the possibility that seq is numeric
if isnumeric(seq)
	seq = int2aa(seq);
elseif ischar(seq)
	seq = upper(seq);   % make sure it is upper case
end

%=== find cleavage sites
sites = regexp(seq, pattern);

lpos = length(position);
if lpos > 1
	if lpos > 2
		warning(message('bioinfo:cleave:OnlyTwoCuts'))
	end
	cut2 = position(2);
	position = position(1);
end

position = position - 1;  % regexp finds position of first match

%=== check that there are no overhangs at either end
seqLength = length(seq);
cutSites = sites + position;
if ~isempty(cut2)
	cut2Sites = sites + cut2 -1;
else
	cut2Sites = cutSites;
end
mask = (cutSites < 1) | (cutSites >= seqLength) | (cut2Sites >= seqLength) ;
sites = sites(~mask);


%=== remove sites to exclude because of exception rule
expSites = regexp(seq, exceptionRule);
sites = setdiff(sites, expSites);
numSites = length(sites);

%=== partial digestion
if ratio < 1    % discard some sites if ratio is not 1
	mask = rand(1,numSites)>ratio;
	sites = sites(mask);
	numSites = length(sites);
end

%=== initialize outputs
parts = cell(numSites+1,1);
if returnCuttingSites
	cuttingSites = zeros(numSites+1,1);
end
if returnLengths
	lengths = zeros(numSites+1,1);
end
if returnMissed
	missed = zeros(numSites+1,1);
end

%=== determine the peptides based on the found sites
pos = 1;
chunk = 1;
for count = 1:numSites
	chunkEnd = sites(count)+position;
	parts{chunk} =  seq(pos:chunkEnd);
	if returnCuttingSites
		cuttingSites(chunk,1) = pos - 1;
	end
	if returnLengths
		lengths(chunk,1) = chunkEnd - pos + 1;
	end
	pos = chunkEnd + 1;
	chunk = chunk + 1;
	if ~isempty(cut2)
		chunkEnd = sites(count)+cut2-1;
		parts{chunk} =  seq(pos:chunkEnd);
		if returnCuttingSites
			cuttingSites(chunk,1) = pos - 1;
		end
		if returnLengths
			lengths(chunk,1) = chunkEnd - pos + 1;
		end
		pos = chunkEnd + 1;
		chunk = chunk + 1;
	end
end

parts{chunk} = seq(pos:seqLength);
if returnCuttingSites
	cuttingSites(chunk,1) = pos - 1;
end
if returnLengths
	lengths(chunk,1) = seqLength - pos + 1;
end


%=== determine peptides resulting from missed cleavages
% The number of peptides when M cleavages are missed is equal to the sum of
% the first M+1 terms of a series with terms a_t = N-t+1, t = 1:M+1. The
% sum is equal to the 1/2 * (a_1 + a_M+1) * (M+1) = 1/2 * (M+1) * (2N-M).

if M > 0
	N = numel(parts); % number of peptides from ideal digestion
	numPep = (M+1) * (2*N-M) * 1/2; % total number of peptides with missed cleavages
	parts(N+1:numPep) = {''};
	
	i = N; % index for filling cell array parts
	k = 1; % index for prefix
	for m = 1:M
		h = m+1; % index for suffix
		next = i + 1; % index for prefix at next iteration of m
		while (h <= N)
			i = i + 1;
			parts{i} = [parts{k} parts{h}];
			
			if returnMissed
				missed(i) = m;
			end
			if returnCuttingSites
				cuttingSites(i,1) = cuttingSites(k);
			end
			if returnLengths
				lengths(i,1) = lengths(k) + lengths(h);
			end
			
			k = k + 1;
			h = h + 1;
		end
		k = next;
		
	end
	
end
%--------------------------------------------------------------------------
function [position, pattern] = localLookupEnzyme(enzyme)
		
%=== try to read enzyme
try
	s = cleavelookup('code', enzyme);
catch theErr1 %#ok<NASGU>
	try
		s = cleavelookup('name', enzyme);
	catch theErr2
		error(message('bioinfo:cleave:UnknownEnzyme', enzyme))
	end
end
	
dummy = textscan(s, '%d%s');
position = double(dummy{1});
pattern = dummy{2}{:};
		
%--------------------------------------------------------------------------
		
function [ratio, M, exceptionRule] = parse_inputs(varargin)
% Parse input PV pairs.
			
			
%=== set defaults
ratio = 1;
M = 0;
exceptionRule = '';

%=== check for the right number of inputs
if rem(nargin,2) == 1
	error(message('bioinfo:cleave:IncorrectNumberOfArguments', mfilename));
end

%=== allowed parameters
okargs = {'partialdigest', 'missedsites', 'exception'};

%=== parse inputs
for j = 1:2:nargin
	pname = varargin{j};
	pval = varargin{j+1};
	k = find(strncmpi(pname, okargs,numel(pname)));
	if isempty(k)
		error(message('bioinfo:cleave:UnknownParameterName', pname));
	elseif length(k)>1
		error(message('bioinfo:cleave:AmbiguousParameterName', pname));
	else
		switch(k)
			case 1  % partialdigest
				ratio = pval;
				if ~isnumeric(ratio) || ratio < 0 || ratio > 1
					error(message('bioinfo:cleave:BadPartialDigestRatio'));
				end
			case 2  % missedsites
				M = pval;
				if ~isscalar(M) || ~isnumeric(M) || M < 0 || rem(M, 1) ~= 0
					error(message('bioinfo:cleave:BadMissedSites'));
				end
			case 3  % exception
				if ischar(pval) && size(pval,1) <= 1
					exceptionRule = pval;
				else
					error(message('bioinfo:cleave:BadExceptionRule'));
				end
				
		end
	end
end





