function info = getInfo(obj, X)
%GETINFO retrieve the information relative to a single element in a BioMap object.
%
%   INFO = GETINFO(OBJ,X) retrieves the information relative to the element
%   indexed by X in a BioMap object OBJ. X can be a numeric scalar, or
%   logical vector with a single true value, or a string corresponding to a
%   valid 'Header' value. INFO is a tab delimited string with the following
%   information: Header Flag Start MappingQuality Signature Sequence
%   Quality.
%
%   Examples:
%
%   % Create a BioMap object.
%   s = samread('ex1.sam');
%   obj = BioMap(s);
%
%   % Retrieve the information for the second element in the object.
%   getInfo(obj, 2)
%   getInfo(obj, {s(2).QueryName})
%   getInfo(obj, [false true])
%
%   See also BIOMAP, BIOMAP/GET, BIOMAP/BIOMAP.

%   Copyright 2009-2012 The MathWorks, Inc.

%=== Input check
bioinfochecknargin(nargin, 2, ['BioMap:' mfilename])
checkScalarInput(obj);

%=== Return empty if there are no sequences
if obj.getNSeqs == 0
	info = '';
	return
end

%=== Determine which elements are indexed by X
if iscellstr(X) || (ischar(X) && isvector(X))
	X = getIndexByHeader(obj, X);
elseif islogical(X) % we need this to use with BioIndexedFile method
	X = find(X);
elseif ~isnumeric(X)
    error(message('bioinfo:BioMap:getInfo:InvalidIndex'));
end

%=== Error checking
if isempty(X)
    error(message('bioinfo:BioMap:getInfo:EmptyIndex'));
end

if ~isscalar(X)
    error(message('bioinfo:BioMap:getInfo:IndexNotScalar'));
end

if X > obj.getNSeqs
    error(message('bioinfo:BioMap:getInfo:IndexOutOfBoundary'));
end

%=== Return info in this order: Header Flag Start MappingQuality Signature Sequence Quality
header = getHeader(obj,X);
flag = getFlag(obj,X);
start = getStart(obj,X);
mappingquality = getMappingQuality(obj,X);
signature = getSignature(obj,X); 
sequence = getSequence(obj,X);
quality = getQuality(obj,X);

info = sprintf('%s\t%d\t%d\t%d\t%s\t%s\t%s', header{1},flag,start,...
               mappingquality,signature{1},sequence{1},quality{1});
           

