function sub = getSubset(obj, X, varargin)
%GETSUBSET retrieve a subset of elements from a BioSeq (or derived) object.
%
%   OUT = GETSUBSET(OBJ, X) retrieves the elements indexed by X in a BioSeq
%   (or derived) object OBJ. X must be a vector of positive integers, a
%   logical vector, or a cell array of strings corresponding to valid
%   'Header' values. OUT is an object of the same class as OBJ.
%
%   OUT = GETSUBSET(..., 'Name', NAME) also assigns the specified value
%   NAME to the property 'Name' in the BioSeq (or derived) object OUT.
%
%   OUT = GETSUBSET(..., 'InMemory', TF) whether the new object OUT is
%   loaded to memory or reuses an existing auxiliary index file. TF
%   defaults to false. TF is ignored when the input object is already
%   loaded to memory.
%     
%   Examples:
%
%   % Create a BioSeq (or derived) object. 
%   obj = BioSeq('Header', {'H1';'H2';'H3'}, ...
%       'Sequence', {randseq(10); randseq(20); randseq(30)});
%
%   % Retrieve the information relative to the last two elements.
%   getSubset(obj, [2 3])
%   getSubset(obj, 2:3)
%   getSubset(obj, [false true true])
%   getSubset(obj, {'H2', 'H3'})
%
%   % Create a BioSeq object with the second element only.
%   out = getSubset(obj, 2, 'Name', 'OBJ_2')
%
%   See also BIOSEQ, BIOSEQ/GET, BIOSEQ/GETHEADER, BIOSEQ/GETSEQUENCE,
%   BIOSEQ/GETSUBSEQUENCE, BIOSEQ/SETSUBSET.

%   Copyright 2009-2012 The MathWorks, Inc.

%=== Input check
bioinfochecknargin(nargin, 2, ['BioSeq:' mfilename])
checkScalarInput(obj);

%=== Parse PVP 
[name,inMemory] = parse_inputs(varargin{:});

%=== Return empty if there are no sequences or the index is empty
if obj.getNSeqs == 0 || isempty(X)
    sub =  eval(class(obj));
    sub.Name = name;
    return
end

X = uniformizeIndexToNumericArray(obj,X);

% Indices must be unique
if numel(unique(X(:)))~=numel(X)
    error(message('bioinfo:BioSeq:getSubset:NonUniqueIndices'))
end

%=== Create new object reusing the Index file
sub = obj;
sub.Index = getSubset(obj.Index,X);
sub.Name = name;
if inMemory
    % Load all data into memory
    props = properties(sub)';
    props = props(~strcmp(props,'NSeqs'));
    sub = [props;get(sub,props)];
    sub = feval(class(obj),sub{:});
end

%--------------------------------------------------------------------------
function [name,inMemory] = parse_inputs(varargin)
% Parse input PV pairs.

%=== defaults
name = '';
inMemory = false;

%=== check for the right number of inputs
if rem(nargin,2) == 1
    error(message('bioinfo:BioSeq:getSubset:IncorrectNumberOfArguments'))
end

%=== allowed parameters
okargs = {'name','inmemory'};

%=== parse inputs
for j = 1:2:nargin
    pname = varargin{j};
    pval = varargin{j+1};
    k = bioinfoprivate.pvpair(pname, pval, okargs, ['BioSeq:' mfilename]);
    
    switch(k)
        case 1 % name
            if ~ischar(pval) || size(pval,1) > 1
                error(message('bioinfo:BioSeq:getSubset:InvalidName'))
            else
                name = pval;
            end
        case 2 % inmemory
            inMemory = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            
    end
end
