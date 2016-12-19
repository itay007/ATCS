function [aln, idx, row] = getCompactAlignment(obj, x1, x2, varargin)
%GETCOMPACTALIGNMENT displays a compact alignment within a given region.
%
%   ALN = GETCOMPACTALIGNMENT(OBJ, X1, X2) returns a char array with the
%   aligned sequences in OBJ within positions X1 and X2 of the reference
%   sequence. X1 and X2 are two valid integers such that X1 <= X2 and both
%   are smaller than the length of the reference sequence in OBJ. ALN is a
%   char array with the aligned sequence reads.
%
%   ALN = GETCOMPACTALIGNMENT(OBJ, X1, X2, R) selects the reference where
%   the alignment is reconstructed. R is either an index to one of the
%   references in the SequenceDictionary property or is a string with the
%   actual reference name.
%
%   [ALN, IDX] = GETCOMPACTALIGNMENT(OBJ, X1, X2) also returns the index of
%   the sequences that align  within positions X1 and X2.
%
%   [ALN, IDX, ROW] = GETCOMPACTALIGNMENT(OBJ, X1, X2) also returns the row
%   where each sequence is best displayed in the compact alignment.
%
%   GETCOMPACTALIGNMENT(..., 'Full', TF) specifies whether only the reads
%   that fully align within the defined region should be considered. TF is
%   a logical scalar. Default is false.
%
%   GETCOMPACTALIGNMENT(..., 'TrimAlignment', TF) Specifies whether or not
%   to trim empty leading and trailing columns from the alignment. Choices
%   are true of false. Default is false, which does not trim the alignment,
%   but includes any empty leading or trailing columns, and returns an
%   alignment always of length X2 - X1 + 1.
%
%   Examples: 
%   %Create BioMap object. 
%   b = BioMap('ex1.sam');
%
%   % Display the alignment within position 1 and 100.
%   getCompactAlignment(b, 1, 100)
%
%   See also BIOMAP, BIOMAP/GETALIGNMENT, BIOMAP/GETBASECOVERAGE.

%   Copyright 2009-2012 The MathWorks, Inc. 

%=== Input check
bioinfochecknargin(nargin, 3, ['BioMap:' mfilename])
checkScalarInput(obj);

if nargin>3 && rem(numel(varargin),2)~=0
    % getCompactAlignment(obj,x1,x2,R,varargin) ?
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
        error(message('bioinfo:BioMap:getCompactAlignment:UnknownReference'))
    end
    obj = getSubset(obj,'SelectReference',R);
elseif numel(obj.SequenceDictionary)>1
    error(message('bioinfo:BioMap:getCompactAlignment:MultipleReferences'));
end

%=== Error check
if ~isnumeric(x1) || ~isscalar(x1) || (isa(x1,'float')&&any(rem(x1,1~=0))) || ...
   ~isnumeric(x2) || ~isscalar(x2) || (isa(x2,'float')&&any(rem(x2,1~=0))) || ...
   any(x1(:)<1) || any(x2(:)<x1(:)) || numel(x1)~=numel(x2)
     error(message('bioinfo:BioMap:getCompactAlignment:InvalidRange'));

end

% Parse input PV pairs.
[fullFlag, trimFlag] = parse_inputs(varargin{:});

%=== determine the most compact way to display the entries
[row, idx, start] = getRowInCompactAlignment(obj, x1, x2, fullFlag);

if isempty(idx)
	aln = [];
	return
end

%=== apply the cigar string to the sequences to find out the gapped
%    sequences that need to be shown in the alignment
gseqs = bioinfoprivate.cigar2gappedsequence(getSequence(obj,idx),getSignature(obj,idx));    

%=== cut sequences that go beyond x2
extra = start+uint32(cellfun(@length,gseqs))-1 - x2;
i = find(extra>0);
for j= 1:numel(i)
    gseqs{i(j)} = gseqs{i(j)}(1:(end-extra(i(j))));
end

%=== cut sequences that go beyond x1
extra = x1-start;
i = find(extra>0);
for j= 1:numel(i)
    gseqs{i(j)} = gseqs{i(j)}(extra(i(j))+1:end);
end

start = max(start,x1) - x1 + 1;
    
%=== initialize matrix for output alignment
space = ' ';  
aln = space(ones((x2-x1+1),max(row),'uint8'));

%=== copy sequence to alignment, columnwise to reduce jumps in memory
len = cellfun(@length,gseqs);
for i =1:numel(row)
    aln(start(i):start(i)+len(i)-1,row(i)) = gseqs{i};
end

%=== transpose alignment
aln = aln';

%=== trim empty columns
if trimFlag
    nonEmptyCol = any(aln~=' ',1);
    aln = aln(:,find(nonEmptyCol,1,'first'):find(nonEmptyCol,1,'last'));
end

%--------------------------------------------------------------------------
function [fullFlag, trimFlag] = parse_inputs(varargin)
% Parse input PV pairs.

%=== defaults
fullFlag = false; % true if consider reads fully overlapping the given range
trimFlag = false; % true if coverage base by base

%=== check for the right number of inputs
if rem(nargin,2) == 1
    error(message('bioinfo:BioMap:getCompactAlignment:IncorrectNumberOfArguments'))
end

%=== allowed parameters
okargs = {'Full', 'TrimAlignment'};

%=== parse parameter value pairs
if nargin > 1
	for j = 1:2:nargin
		pname = varargin{j};
		pval = varargin{j+1};
		k = bioinfoprivate.pvpair(pname, pval, okargs, ['BioMap:' mfilename]);
		switch(k)
			case 1  % full
				fullFlag = bioinfoprivate.opttf(pval, okargs{k}, ['BioMap:' mfilename]);
			case 2  % trim
				trimFlag = bioinfoprivate.opttf(pval, okargs{k}, ['BioMap:' mfilename]);		
		end
	end
elseif nargin ~= 0
    error(message('bioinfo:BioMap:getCompactAlignment:InvalidPVP'));
end
	
