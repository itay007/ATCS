function obj = combine(obj1, obj2, varargin)
%COMBINE combine two objects of BioSeq class.
%
%   OBJ = COMBINE(OBJ1, OBJ2) combines two BioSeq objects OBJ1 and OBJ2
%   into one single object OBJ. The elements in OBJ are the same as the
%   elements in OBJ1 followed by the elements in OBJ2.
%
%   OBJ = COMBINE(..., 'Name', NAME) also assigns the specified value
%   NAME to the property 'Name'  of the new object OBJ.
%
%   Example:
%
%   % Create two structures with fields 'Header' and 'Sequence'.
%   s1 = struct('Header', {'H1';'H2';'H3'}, 'Sequence', {randseq(10);
%   randseq(20); randseq(30)});
%   s2 = struct('Header', {'H4';'H5';'H6'}, 'Sequence', {randseq(10);
%   randseq(20); randseq(30)});
%
%   % Create two BioSeq objects.
%   obj1 = BioSeq(s1)
%   obj2 = BioSeq(s2)
%
%   % Combine the two objects and assign a new name.
%   obj3 = combine(obj1, obj2, 'Name', 'obj1 + obj2')
%
%   See also BIOSEQ, BIOSEQ/GETSUBSET, BIOSEQ/BIOSEQ, BIOSEQ/SETSUBSET.

%   Copyright 2009-2012 The MathWorks, Inc.

%=== Input check
bioinfochecknargin(nargin, 2, ['BioSeq:' mfilename])

%=== Parse PVP 
name = parse_inputs(varargin{:});

%=== Make sure both objects are of same class
if ~isequal(class(obj1), class(obj2))
    error(message('bioinfo:BioSeq:combine:InvalidClassPair'));
end
checkScalarInput(obj1);
checkScalarInput(obj2);

obj =  eval(class(obj1));
obj.Name = name;

props = [{'Sequence'};setdiff(properties(obj1),'Sequence')];
for i = 1:numel(props)
    if isValidField(obj1.Index,props{i})
        v1 = obj1.(props{i});
        v2 = obj2.(props{i});
        if isempty(v1)~=isempty(v2)
            error(message('bioinfo:BioSeq:combine:EmptyProperty',props{i}));
        elseif ~isempty(v1)
            obj.(props{i}) = [v1;v2];
        end
    end
end

%--------------------------------------------------------------------------
function name = parse_inputs(varargin)
% Parse input PV pairs.

%=== defaults
name = '';

%=== check for the right number of inputs
if rem(nargin,2) == 1
    error(message('bioinfo:BioSeq:combine:IncorrectNumberOfArguments'))
end

%=== allowed parameters
okargs = {'name'};

%=== parse inputs
for j = 1:2:nargin
	pname = varargin{j};
	pval = varargin{j+1};
	
	k = bioinfoprivate.pvpair(pname, pval, okargs, ['BioSeq:' mfilename]);
	switch(k)
		case 1 % name
			if ~ischar(pval) || size(pval,1) > 1
                error(message('bioinfo:BioSeq:combine:InvalidName'))
			else
				name = pval;
			end
	end
end
