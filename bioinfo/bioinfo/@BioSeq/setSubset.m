function obj = setSubset(obj, Y, X)
%SETSUBSET set the values of a subset of elements in a BioSeq (or derived) object.
%
%   SETSUBSET(OBJ1, OBJ2, X) sets the elements of a BioSeq (or derived)
%   object OBJ1 indexed by X to the values specified by the BioSeq (or
%   derived) object OBJ2. X can be a numeric or logical vector, or a cell
%   array of strings corresponding to valid 'Header' values. OBJ2 must be a
%   valid BioSeq object of appropriate size. OBJ1 and OBJ2 must be of the
%   same class.
%
%   NOTE: SETSUBSET(OBJ1,...) without assigning to a variable does not
%   modify the object's property.  Use OBJ1 = SETSUBSET(OBJ1,...) to
%   modify OBJ1.
%
%   Examples:
%
%   % Create two BioSeq (or derived) objects.
%   obj1 = BioSeq('Header', {'H1';'H2';'H3'}, ...
%       'Sequence', {randseq(10); randseq(20); randseq(30)});
%
%   obj2 = BioSeq('Header', {'H4';'H5'}, ...
%       'Sequence', {randseq(10); randseq(20)});
%
%   % Create copies of obj1 with the first two elements equal to obj2.
%   out1 = setSubset(obj1, obj2, 1:2)
%   out2 = setSubset(obj1, obj2, [1 2])
%   out3 = setSubset(obj1, obj2, [true true false])
%   out4 = setSubset(obj1, obj2, {'H1', 'H2'})
%
%   % Set the first two elements of obj1 equal to obj2.
%   obj1 = setSubset(obj1, obj2, 1:2)
%
%   See also BIOSEQ, BIOSEQ/GETSUBSET, BIOSEQ/SET, BIOSEQ/SETHEADER,
%   BIOSEQ/SETSEQUENCE, BIOSEQ/SETSUBSEQUENCE, BIOSEQ/SETSUBSET.

%   Copyright 2009-2012 The MathWorks, Inc.

%=== Input check
bioinfochecknargin(nargin, 2, ['BioSeq:' mfilename])
checkScalarInput(obj);

if ~isa(Y, class(obj))
    error(message('bioinfo:BioSeq:setSubset:InvalidInputObject'));
end
checkScalarInput(Y);

[X,ex] = uniformizeIndexToNumericArray(obj,X);

% Indices must be unique
if numel(unique(X(:)))~=numel(X)
    error(message('bioinfo:BioSeq:setSubset:NonUniqueIndices'))
end

N = length(X);

%=== check input size
if N ~= getNSeqs(Y)
	if isempty(ex)
        error(message('bioinfo:BioSeq:setSubset:InvalidInputObjectSize'));
    else
        error(message('bioinfo:BioSeq:setSubset:DuplicateHeadersWhenSettingSubset'));
	end
end

%=== check for out-of-bounds
if any(X > obj.getNSeqs)
    error(message('bioinfo:BioSeq:setSubset:BadIndex'));
end

%=== set new values
props = properties(obj);
Y = get(Y);
for i = 1:numel(props)
    if isValidField(obj.Index,props{i})
        e1 = isempty(obj.(props{i}));
        e2 = isempty(Y.(props{i}));
        if e1~=e2
            error(message('bioinfo:BioSeq:setSubset:EmptyProperty',props{i}));
        elseif ~e1
            obj.(props{i})(X) = Y.(props{i});
        end
    end
end