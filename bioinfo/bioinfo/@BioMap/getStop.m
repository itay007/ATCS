function out = getStop(obj, X)
%GETSTOP compute the stop position of each alignment in a BioMap object. 
%
%   OUT = GETSTOP(OBJ) computes the stop position for the aligned
%   sequences in a BioMap object. OUT is an array of integers.
%
%   OUT = GETSTOP(OBJ, X) computes the stop position for the aligned
%   sequences indexed by X in a BioMap object. X must be a vector
%   of positive integers, a logical vector, or a cell array of strings
%   corresponding to valid 'Header' values. OUT is an array of integers.
%
%   Examples:
%
%   % Create a BioMap object.
%   s = samread('ex1.sam');
%   obj = BioMap(s);
%
%   % Compute the position where the alignment of the second sequence in
%   % the object ends (with respect to the reference sequence).
%   getStop(obj, 2)
%   getStop(obj, {s(2).QueryName})
%   getStop(obj, [false true])
%
%   % Compute the position where the alignment of the 1st and 3rd sequences
%   % in the object end (with respect to the reference sequence).
%   
%   getStop(obj, [1 3])
%   getStop(obj, {s([1 3]).QueryName})
%   getStop(obj, [true false true])
%
%   % Compute the position where each alignment in the object ends.
%   getStop(obj)
%
%   See also BIOMAP, BIOMAP/GET, BIOMAP/GETSTART.

%   Copyright 2009-2012 The MathWorks, Inc.

%=== Input check
bioinfochecknargin(nargin, 1, ['BioMap:' mfilename])
checkScalarInput(obj);
%=== Validate and clean X when it is specified
if nargin>1
	%=== transform X into numeric array
	if isempty(X)
		X = [];
	elseif islogical(X) && isvector(X)
		X = find(X);
	elseif iscellstr(X) || (ischar(X) && isvector(X))
		X = getIndexByHeader(obj, X);
	elseif ~isnumeric(X) || ~isvector(X)
        error(message('bioinfo:BioMap:getStop:InvalidInput'));
	end    
end

%=== Stop is not an object property, but in some data adaptors have this
%    information available as a field:
if obj.Index.isValidField('Stop')
    if nargin==1
        out =  obj.Index.getField('Stop');
    else
        out =  obj.Index.getField('Stop',X);
    end
    return
end

%=== If X is not specified, then all elements are to be considered
if nargin == 1 
	X = 1:obj.getNSeqs;
end

%=== return empty if X is empty
if isempty(X)
	out = [];
	return
end

%=== determine the end positions
	
start = getStart(obj, X);
cigar = getSignature(obj, X);
e = bioinfoprivate.cigar2endmex(cigar);
		
%=== compute stop
out = start + e - 1;
	
%=== reset entries with e == 0 (no informative cigar found)
out(e == 0) = start(e == 0);
