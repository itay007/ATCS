function [aln, idx] = getAlignment(obj, x1, x2, varargin)
%GETALIGNMENT reconstruct the alignment within a given region.
%
%   ALN = GETALIGNMENT(OBJ, X1, X2) returns the alignment of the sequences
%   represented in OBJ within the position range specified by X1 and X2. X1
%   and X2 are two valid integers such that X1 <= X2 and both are smaller
%   than the length of the reference sequence in OBJ. ALN is a char
%   array with the aligned sequence reads. 
%
%   [ALN, IDX] = GETALIGNMENT(OBJ, X1, X2) also returns the index of the
%   sequences that align within the given region.
%
%   GETALIGNMENT(OBJ, X1, X2, R) selects the reference where the alignment
%   is reconstructed. R is either an index to one of the references in the
%   SequenceDictionary property or is a string with the actual reference
%   name. 
%
%   GETALIGNMENT(..., 'OffsetPad', TF) specifies whether padding
%   blanks are added at the beginning of each aligned sequence to represent
%   the offset of the start position with respect to the reference. TF is a
%   logical scalar. Default is false.
%
%   Examples:
%   % Create a BioMap object.
%   obj = BioMap('ex1.sam');
%
%   % Reconstruct the alignment between positions 10 and 25. 
%   getAlignment(obj, 10, 25)
%
%   See also BIOMAP, BIOMAP/GETCOMPACTALIGNMENT, BIOMAP/GETBASECOVERAGE,
%   BIOMAP/GETINFO. 

%   Copyright 2009-2012 The MathWorks, Inc.

%=== Input check
bioinfochecknargin(nargin, 3, ['BioMap:' mfilename])
checkScalarInput(obj);

if nargin>3 && rem(numel(varargin),2)~=0
    % getAlignment(obj,x1,x2,R,varargin) ?
    R = varargin{1};
    varargin = varargin(2:end);
    % validate R
    if iscellstr(R) && numel(R)==1
        R = R{1};
    end
    if ischar(R) && isrow(R)
        R = find(strcmp(obj.SequenceDictionary,R));
    end
    if ~isnumeric(R) || ~isscalar(R) || R<1 || R>numel(obj.SequenceDictionary) || rem(R,1)
        error(message('bioinfo:BioMap:getAlignment:UnknownReference'))
    end
    obj = getSubset(obj,'SelectReference',R);
elseif numel(obj.SequenceDictionary)>1
    error(message('bioinfo:BioMap:getAlignment:MultipleReferences'));
end


%=== Error check
if ~isnumeric(x1) || ~isscalar(x1) || (isa(x1,'float')&&any(rem(x1,1~=0))) || ...
   ~isnumeric(x2) || ~isscalar(x2) || (isa(x2,'float')&&any(rem(x2,1~=0))) || ...
   any(x1(:)<1) || any(x2(:)<x1(:)) || numel(x1)~=numel(x2)
    error(message('bioinfo:BioMap:getAlignment:InvalidRange'));
end

%=== Parse PVP 
doOffsetPad = parse_inputs(varargin{:});

%=== Find sequences overlapping with given range
idx = getIndex(obj,x1(:),x2(:));

if isempty(idx)
    aln = [];
	return
end

%=== Reconstruct alignment
aln = cigar2align(getSequence(obj, idx), getSignature(obj, idx), ...
	'start', getStart(obj, idx), 'OffsetPad', doOffsetPad);
if ~doOffsetPad

	k1 = min(getStart(obj, idx));
	
	x1 = x1 - k1 + 1; % range w/r to the beginning of the alignment
	x2 = x2 - k1 + 1;
	
	k2 = size(aln, 2);
	
	if x1 < 0
		x1 = 1;
	end
	
	if x2 > k2
		x2 = k2;
	end
end

if size(aln,2) < x2
	x2 = size(aln,2);
end

aln =  aln(:, x1:x2);

	
%--------------------------------------------------------------------------
function doOffsetPad = parse_inputs(varargin)
% Parse input PV pairs.

%=== defaults
doOffsetPad = false;

%=== check for the right number of inputs
if rem(nargin,2) == 1
    error(message('bioinfo:BioMap:getAlignment:IncorrectNumberOfArguments'))
end

%=== allowed parameters
okargs = {'offsetPad'};

%=== parse inputs
for j = 1:2:nargin
	pname = varargin{j};
	pval = varargin{j+1};
	k = bioinfoprivate.pvpair(pname, pval, okargs, ['BioMap:' mfilename]);
	switch(k)
		case 1 % offsetPad
			doOffsetPad = bioinfoprivate.opttf(pval, okargs{k}, ['BioMap:' mfilename]);
	end
end



