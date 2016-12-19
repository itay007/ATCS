function sub = getSubset(obj, X, varargin)
%GETSUBSET retrieve a subset of elements from a BioMap object.
%
%   OUT = GETSUBSET(OBJ, X) retrieves the elements indexed by X in a BioMap
%   object OBJ. X must be a vector of positive integers, a logical vector,
%   or a cell array of strings corresponding to valid 'Header' values. OUT
%   is a BioMap object.
%
%   OUT = GETSUBSET(OBJ,'SelectReference',R) creates a subset object with
%   only the short reads mapped to R. R can be a vector with indices to the
%   SequenceDictionary property or a cell string with the actual names of
%   the references.
%
%   OUT = GETSUBSET(..., 'InMemory', TF) whether the new object OUT is
%   loaded to memory or reuses an existing auxiliary index file. TF
%   defaults to false. TF is ignored when the input object is already
%   loaded to memory.
%
%   OUT = GETSUBSET(..., 'Name', NAME) also assigns the specified value
%   NAME to the property 'Name'  of the new object OUT.
%
%   Examples:
%
%   % Create a BioMap object.
%   s = samread('ex1.sam');
%   obj = BioMap(s);
%
%   % Retrieve the information relative to the 2nd and 3rd elements.
%   getSubset(obj, [2 3])
%   getSubset(obj, 2:3)
%   getSubset(obj, [false true true])
%   getSubset(obj, {s([2 3]).QueryName})
%
%   % Create a BioMap object with the first five element only.
%   out = getSubset(obj, 1:5, 'Name', 'Obj_1_5')
%
%   See also BIOMAP, BIOMAP/GET, BIOMAP/SETSUBSET.


%   Copyright 2009-2012 The MathWorks, Inc.

%=== Input check
bioinfochecknargin(nargin, 2, ['BioMap:' mfilename])
checkScalarInput(obj);

%=== If second input arg is a char then check for 
%    getSubset(obj,'SelectReference',...) signature
subsetByRef = false;
if ischar(X) && isrow(X)
    subsetByRef = true;
    varargin = [{X},varargin];
end

%=== Parse PVP 
[name,loadInMemory,references] = parse_inputs(varargin{:});

%=== Check that signatures do not collide
if subsetByRef && isempty(references)
    error(message('bioinfo:BioMap:getSubset:InvalidIndexOrReference'))
elseif ~subsetByRef && ~isempty(references)
    error(message('bioinfo:BioMap:getSubset:InvalidIndexOrReference'))
end
    
if subsetByRef %= Subset by reference names

    %=== Validate references
    if isnumeric(references)
        if any(references<1) || any(references>numel(obj.SequenceDictionary))
            error(message('bioinfo:BioMap:getSubset:UnknownReference'))
        end
        references = obj.SequenceDictionary(references(:));
    end
    references = unique(references);
    href = false(1,numel(obj.SequenceDictionary));
    for i=1:numel(references)
        k = strcmp(obj.SequenceDictionary,references{i});
        if any(k)
            href(k) = true;
        else
            error(message('bioinfo:BioMap:getSubset:UnknownReference'))
        end
    end
    
    %=== Create new object reusing the Index file
    sub = obj;
    if ~isempty(sub.DictionaryMapping)
        newIdxToSeqDic = cumsum(href(:)).* href(:);
        g = sub.DictionaryMapping>0;
        sub.DictionaryMapping(g) = newIdxToSeqDic(sub.DictionaryMapping(g));
        sub.Index = getSubsetReference(obj.Index,find(sub.DictionaryMapping>0));
        sub.SequenceDictionary = sub.SequenceDictionary(href);
    elseif obj.Index.InMemory 
        sub.Index = getSubsetReference(obj.Index,find(href));
        sub.SequenceDictionary = sub.SequenceDictionary(href);
    else 
        sub.Index = getSubsetReference(obj.Index,references);        
        sub.SequenceDictionary = sub.SequenceDictionary(href);
    end
    sub.Name = name;
    if sub.NSeqs==0
        sub = BioMap('SequenceDictionary',sub.SequenceDictionary,'Name',name);
        loadInMemory = false; % it is already in memory
    end
    
    if loadInMemory % Load all data into memory
        props = properties(sub)';
        props = props(~strcmp(props,'NSeqs'));
        sub = [props;get(sub,props)];
        sub = BioMap(sub{:});
    end

else
    % Subset by indexing: call upper class implementation
    sub = getSubset@BioSeq(obj,X,'Name',name,'InMemory',loadInMemory);    
end

%--------------------------------------------------------------------------
function [name,inMemory,references] = parse_inputs(varargin)
% Parse input PV pairs.

%=== defaults
name = '';
inMemory = false;
references = '';

%=== check for the right number of inputs
if rem(nargin,2) == 1
    error(message('bioinfo:BioMap:getSubset:IncorrectNumberOfArguments'))
end

%=== allowed parameters
okargs = {'name','inmemory','selectreference'};

%=== parse inputs
for j = 1:2:nargin
    pname = varargin{j};
    pval = varargin{j+1};
    
    k = bioinfoprivate.pvpair(pname, pval, okargs, ['BioMap:' mfilename]);
    switch(k)
        case 1 % name
            if ~ischar(pval) || size(pval,1) > 1
                error(message('bioinfo:BioMap:getSubset:InvalidName'))
            else
                name = pval;
            end
        case 2 % inmemory
            inMemory = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 3 % selectreference
            if isnumeric(pval) && isvector(pval) && all(~rem(pval,1))
                references = pval;
            elseif ischar(pval) && isrow(pval)
                references = {pval};
            elseif iscellstr(pval) && numel(pval)==numel(unique(pval))
                references = pval;
            else
                error(message('bioinfo:BioMap:getSubset:InvalidReference'))
            end
                
    end
end



