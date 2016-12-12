function [gs,outstarts,aln] = cigar2align(seqs, cigars, varargin)
%CIGAR2ALIGN convert unaligned sequences to aligned sequences using CIGAR strings.
%
%   ALN = CIGAR2ALIGN(S, C) converts the unaligned sequences in cell array
%   S into a matrix of aligned sequences using the CIGAR strings stored in
%   the cell array C. S and C must have the same size. If S has n cells,
%   then the output ALN is a char array with n rows.
%
%   [GAPSEQ, IDX] = CIGAR2ALIGN(S, C) with two output arguments returns a
%   cell array GAPSEQ containing the aligned sequences, without any leading
%   or trailing whitespace, and a numeric vector IDX indicating the
%   starting column for each aligned sequence in ALN. Beware that IDX is
%   not necessarily the same as the start positions of the reads in the
%   reference, since the reference may be extended to account for
%   insertions or the reads may have leading soft clipping, padding or
%   insertion symbols.
%
%   CIGAR2ALIGN(..., 'Start', POS) specifies the reference positions at
%   which the alignment for each sequence starts. By default the alignment
%   of each sequence starts at position 1 in the reference. POS is a
%   numeric array of positive integers of size equal to the number of
%   sequences to consider.
%
%   CIGAR2ALIGN(..., 'GapsInRef', TF) specifies whether positions
%   corresponding to gaps in the reference sequence should be included or
%   not. TF is a logical scalar. Default is false.
%
%   CIGAR2ALIGN(..., 'SoftClipping', TF) specifies whether positions
%   corresponding to soft clipping ends are included or not. TF is a
%   logical scalar. Default is false.
%
%   CIGAR2ALIGN(..., 'OffsetPad', TF) specifies whether padding blanks are
%   added at the beginning of each aligned sequence to represent the offset
%   of the start position with respect to the reference. TF is a logical
%   scalar. Default is false.
%
%   Examples:
%   r = {'ACGACTGC', 'ACGTTGC', 'AGGTATC'}; % unaligned sequences
%   c = {'3M1D1M1I3M', '4M1D1P3M', '5M1P1M1D1M'}; % cigar strings
%
%   % Reconstruct the alignment.
%   aln1 = cigar2align(r, c)
%
%   % Reconstruct the alignment using the corresponding offsets.
%   aln2 = cigar2align(r, c, 'start', [5 5 5], 'OffsetPad', true)
%
%   See also ALIGN2CIGAR, BIOMAP/GETALIGNMENT, BIOMAP/GETCOMPACTALIGNMENT,
%   BIOMAP/GETCOVERAGE.

%   Copyright 2009-2010 The MathWorks, Inc.


%=== Check inputs
bioinfochecknargin(nargin, 2, mfilename)

if ~iscellstr(seqs) || isempty(seqs)
	error(message('bioinfo:cigar2align:InvalidSequenceInput'));
end

if ~iscellstr(cigars)|| isempty(cigars)
	error(message('bioinfo:cigar2align:InvalidCigarInput'));
end

if numel(seqs) ~= numel(cigars)
	error(message('bioinfo:cigar2align:InvalidInputSize'));
end

%=== Parse PVP
[starts, doOffsetPad, doGapsInRef, includeSoftClipping] = parse_inputs(varargin{:});

%=== Call cigar2gappedsequence
[gs,ap,ri] = bioinfoprivate.cigar2gappedsequence(seqs,cigars,'GapsInRef',doGapsInRef,'softClipping',includeSoftClipping);

%=== Initialize variables
n = numel(cigars);

if ~isempty(starts)
    starts = double(starts); % to assure correctness of arithmetic operations downstream
    if numel(starts) ~= numel(seqs)
        error(message('bioinfo:cigar2align:InvalidStartSize'));
    end
else
    starts = ones(1, n);
end

%=== Calculate where each of the gapped sequences start in the ungapped reference
refstarts = starts(:)-ap+1;

%=== The first position required from referece
refOffset = min(refstarts);

%=== Accumulate openings needed in the reference
RI = zeros(1, max(cellfun(@length,gs))+max(refstarts)-refOffset+1); %large enough for worst case, not necessarily needed
for i = 1:n
    idx = refstarts(i)-refOffset+(1:length(ri{i}));
    RI(idx) = max(RI(idx),ri{i});
end
RI = [0 cumsum(RI)];

%=== Calculate the start for every gapped sequence in the gapped reference
if doOffsetPad
    outstarts = refstarts + RI(starts(:)-refOffset+1)';
else
    outstarts = refstarts + RI(starts(:)-refOffset+1)' - refOffset + 1;
end

%=== With two output arguments we do not need to create a padded matrix 
if nargout==2
    return
end

%=== Calculate end for the actual output alignment 
outends = outstarts + cellfun(@length,gs)-1;

%=== Now fill the output padded matrix
padSymbol = ' ';
aln = padSymbol(ones(n,max(outends),'uint8'));
for i = 1:n
    aln(i,outstarts(i):outends(i)) = gs{i};
end

%=== With one or none output arguments we return the padded matrix
if nargout<2
    gs = aln;
end

%--------------------------------------------------------------------------

function [starts, doOffsetPad, doGapsInRef, includeSoftClipping] = parse_inputs(varargin)
% Parse input PV pairs.

%=== defaults
starts = [];
doOffsetPad = false;
doGapsInRef = false;
includeSoftClipping = false;

%=== check for the right number of inputs
if rem(nargin,2) == 1
	error(message('bioinfo:cigar2align:IncorrectNumberOfArguments', mfilename));
end

%=== allowed parameters
okargs = {'start', 'offsetPad', 'gapsInRef', 'softClipping'};

%=== parse inputs
for j = 1:2:nargin
	pname = varargin{j};
	pval = varargin{j+1};
	
	if ~ischar(pname)
		error(message('bioinfo:cigar2align:InvalidParameter'));
	end
	
	k = find(strncmpi(pname, okargs, numel(pname)));
	if isempty(k)
		error(message('bioinfo:cigar2align:UnknownParameterName', pname));
	elseif length(k)>1
		error(message('bioinfo:cigar2align:AmbiguousParameterName', pname));
	else
		switch(k)
			case 1 % start
				if isnumeric(pval) && all(pval > 0)
					starts = pval;
				else
					error(message('bioinfo:cigar2align:InvalidStart'));
				end
			case 2 % offsetPad
				doOffsetPad = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
			case 3 % gapsInRef
				doGapsInRef = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
            case 4 % softClipping
                includeSoftClipping = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
                
		end
	end
end





