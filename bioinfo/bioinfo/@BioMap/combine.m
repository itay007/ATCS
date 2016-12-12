function obj = combine(obj1, obj2, varargin)
%COMBINE combine two objects of BioMap class.
%
%   OBJ = COMBINE(OBJ1, OBJ2) combines two BioMap objects OBJ1 and OBJ2
%   into one single object OBJ. The elements in OBJ are the same as the
%   elements in OBJ1 followed by the elements in OBJ2.
%
%   OBJ = COMBINE(..., 'Name', NAME) also assigns the specified value
%   NAME to the property 'Name'  of the new object OBJ.
%
%   Example:
%
%   % Create two structures with data from a SAM file.
%   s1 = samread('ex1.sam', 'blockread', [1 10]);
%   s2 = samread('ex1.sam', 'blockread', [11 20]);
%
%   % Create two BioMap objects.
%   obj1 = BioMap(s1)
%   obj2 = BioMap(s2)
%
%   % Combine the two objects and assign a new name.
%   obj3 = combine(obj1, obj2, 'Name', 'obj1 + obj2')
%
%   See also BIOMAP, BIOMAP/GETSUBSET, BIOMAP/BIOMAP, BIOMAP/SETSUBSET.

%   Copyright 2009-2012 The MathWorks, Inc.

%=== Input check
bioinfochecknargin(nargin, 2, ['BioMap:' mfilename])

%=== Make sure both objects are of same class
if ~isequal(class(obj1), class(obj2))
    error(message('bioinfo:BioMap:combine:InvalidClassPair'));
end
checkScalarInput(obj1);
checkScalarInput(obj2);

%=== Error out if at least one of the objects has more than one reference
%    in SequenceDictionary
if numel(obj1.SequenceDictionary)>1 || numel(obj2.SequenceDictionary)>1 
    error(message('bioinfo:BioMap:combine:InvalidCombineMultipleReferences'));
end

%=== Make sure both have same reference
if ~strcmp(obj1.SequenceDictionary{1}, obj2.SequenceDictionary{1})
    error(message('bioinfo:BioMap:combine:DifferentReference'));
end

obj = combine@BioSeq(obj1,obj2,varargin{:});
