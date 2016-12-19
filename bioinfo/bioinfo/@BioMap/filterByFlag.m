function idx = filterByFlag(obj, varargin)
%FILTERBYFLAG filter the elements of a BioMap object according to specific
%flag criteria.
%
%   IDX = FILTERBYFLAG(OBJ, FLAGNAME, FLAGVALUE) returns a logical
%   index of those entries in OBJ with the specified FLAGNAME equal to
%   FLAGVALUE. IDX is logical array of size equal to the number of elements
%   in OBJ, with a true value in correspondence of those elements that
%   satisfy the specified criteria. FLAGVALUE is a logical scalar. Valid
%   choices for FLAGNAME are the following strings:
%            'pairedInSeq': the read is paired in sequencing
%            'pairedInMap': the read is mapped in a proper pair
%          'unmappedQuery': the read is unmapped
%           'unmappedMate': the mate is unmapped
%            'strandQuery': strand of the read (0 = forward, 1 = reverse)
%             'strandMate': strand of the mate (0 = forward, 1 = reverse)
%            'readIsFirst': the read is the first in a pair
%           'readIsSecond': the read is the second in a pair
%          'alnNotPrimary': the alignment is not primary
%        'failedQualCheck': the read fails platform/vendor quality checks
%              'duplicate': the read is a PCR or optical duplicate
%
%   FILTERBYFLAG(OBJ, IDX2, FLAGNAME, FLAGVALUE) checks the flag only in a
%   subset of the entries which are indexed by IDX2. IDX2 can be a logical
%   or numeric array of non-negative integers.
%
%   FILTERBYFLAG(..., FLAGNAME1, FLAGVALUE1, FLAGNAME2, FLAGVALUE2,....)
%   applies multiple filter criteria with a single statements.
%
%   Examples:
%   % Create a BioMap object.
%   obj = BioMap('ex1.sam');
%   
%   % Determine which read is mapped in a proper pair and is the first in a
%   % pair.
%   x = filterByFlag(obj, 'pairedInMap', true, 'readIsFirst', true);
%   obj.Header(x)
%
%   See also BIOMAP, BIOMAP/BIOMAP, BIOMAP/GETQUALITY.

%   Copyright 2009-2012 The MathWorks, Inc.

%=== Input check
bioinfochecknargin(nargin, 3, ['BioMap:' mfilename])
checkScalarInput(obj);

%=== Get flags from the object
if isnumeric(varargin{1}) || islogical(varargin{1})
	flags = obj.getFlag(varargin{1});
	varargin = varargin(2:end);
else
	flags = obj.getFlag;
end

%=== Get criteria and values (see description of FLAG in SAM format)
[b, v] = parse_inputs(varargin{:});
b = find(b); % bits to consider in the flag
v = v(b); % values of bits to consider in the flag

%=== Determine elements that satisfy the criteria
idx = true(size(flags)); 
for i = 1:length(b)
	idx = idx & bitget(flags, b(i)) == v(i); 
end

%--------------------------------------------------------------------------
function [flagBit, flagValue] = parse_inputs(varargin)
% Parse input PV pairs.

%=== check for the right number of inputs
if rem(nargin,2) == 1
    error(message('bioinfo:BioMap:filterByFlag:IncorrectNumberOfArguments'))
end

%=== allowed parameters
okargs = {'pairedInSeq', 'pairedInMap', 'unmappedQuery', 'unmappedMate', ...
	'strandQuery', 'strandMate', 'readIsFirst', 'readIsSecond', ...
	'alnNotPrimary', 'failedQualCheck', 'duplicate'};

flagBit = false(1,numel(okargs));
flagValue = false(1, numel(okargs));

%=== parse parameter value pairs
for j = 1:2:nargin
	pname = varargin{j};
	pval = varargin{j+1};
	
	k = bioinfoprivate.pvpair(pname, pval, okargs, ['BioMap:' mfilename]);
	
	flagBit(k) = true;
	flagValue(k) = bioinfoprivate.opttf(pval, okargs{k}, ['BioMap:' mfilename]);
end






	
	
	
