function [cigars, starts] = align2cigar(aln, ref)
%ALIGN2CIGAR convert aligned sequences to the corresponding CIGAR strings.
%
%   [CIGARS, STARTS] = ALIGN2CIGAR(ALN, REF) converts the multiple
%   alignment represented in ALN into the corresponding CIGAR strings based
%   on the aligned reference sequence REF. ALN can be a cell array of
%   strings or a char array. REF must be a string of length equal to each
%   aligned sequence. CIGARS is a cell array of strings. STARTS is a vector
%   of integers containing the starting position of each aligned sequence
%   with respect to the ungapped reference REF. Lower case symbols are
%   interpreted as Soft clipped (S in CIGAR). A dot represents a skipped
%   position (N is CIGAR).
%
%   Examples:
%   aln = ['ACG-ATGC'; 'ACGT-TGC'; '  GTAT-C'];
%   ref =  'ACGTATGC';
%   [cigar, start] = align2cigar(aln, ref)
%
%   See also BIOMAP/GETALIGNMENT, BIOMAP/GETCOMPACTALIGNMENT,
%   BIOMAP/GETCOVERAGE, CIGAR2ALIGN.

%   Copyright 2010-2010 The MathWorks, Inc.



%=== Check inputs
bioinfochecknargin(nargin, 2, mfilename)

if iscell(aln)
	aln = char(aln);
end

if ~ischar(aln) || isempty(aln)
	error(message('bioinfo:align2cigar:InvalidAlignmentInput'));
end

if ~ischar(ref) || ~isvector(ref) || isempty(ref)
	error(message('bioinfo:align2cigar:InvalidReferenceInput'));
end

if length(ref) ~= size(aln,2)
	error(message('bioinfo:align2cigar:InvalidInputSize'));
end

n = size(aln, 1); % number of sequences

%=== Initialize variables
starts = ones(1, n);
cigars = cell(1, n);
soft = zeros(1,n); % keep track of first soft clipping in each aln

%=== Reconstruct CIGARS from alignment
for i = 1:n
	
	%=== determine paddings and remove from seq and ref
	[starts(i),e] = regexp(aln(i,:), '[\w-\.]+', 'start', 'end');
	currAln = aln(i,starts(i):e);
	currRef = ref(starts(i):e);
		
	%=== determine the type of operations 
	x = zeros(size(currAln));  % M (Match/Mismatch)
	x(currAln == '-') = 1; % D
	x(currRef == '-') = x(currRef == '-') + 3; % I (3) or P (4)
	x(currAln == '.') = 2; % N
	x(currAln ~= upper(currAln)) = 5; % S
	
			
	%=== run the RLE (run length encoding) algorithm 
	j = [find(x(1:end-1) ~= x(2:end)) length(x)];
	v = diff([0 j]); % number of operation
	t = x(j); % type of operation
		
	cigars{i} = '';
	%=== determine each operation symbol
	for k = 1:length(v)
		switch t(k)
			case 0 
				cigars{i} = strcat(cigars{i} , sprintf('%dM', v(k)));
			case 1
				cigars{i} = strcat(cigars{i} , sprintf('%dD', v(k)));
			case 2
				cigars{i} = strcat(cigars{i} , sprintf('%dN', v(k)));
			case 3
				cigars{i} = strcat(cigars{i} , sprintf('%dI', v(k)));
			case 4
				cigars{i} = strcat(cigars{i} , sprintf('%dP', v(k)));
			case 5
				cigars{i} = strcat(cigars{i} , sprintf('%dS', v(k)));
				if ~soft(i)
					soft(i) = v(k);
				end
		end
	end
	
	%=== get starts with respect to un-gapped reference
	refGaps = ref == '-';
	starts(i) = starts(i) - sum(refGaps(1:starts(i)));
end

starts = starts + soft; % in SAM specifications, start position does not include S positions

%--------------------------------------------------------------------------
