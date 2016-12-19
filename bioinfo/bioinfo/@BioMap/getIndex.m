function idx = getIndex(obj, x1, x2, varargin)
%GETINDEX finds short reads in a range.
%
%   IDX = GETINDEX(OBJ,X1,X2) returns in IDX the index of the sequences
%   that cover the alignment in the region specified by X1 and X2. X1 and
%   X2 are two non-negative integers such that X1 <= X2 and both are
%   smaller than the length of the reference sequence. X1 and X2 may also
%   be vectors of the same length representing a segmented range.
%
%   GETINDEX(OBJ,X1,X2,R) selects the reference associated to the range
%   specified by X1 and X2. R is either an index to one of the references
%   in the SequenceDictionary property or is a string with the actual
%   reference name. 
%
%   GETINDEX guaranties that the values in IDX are unique; therefore
%   short reads in overlapping segments or short reads that span over two
%   different segments in X1 and X2 are considered only once.
%
%   GETINDEX(...,'overlap',BP) specifies the minimum number of positions
%   that a short read must overlap the given range (or segmented range) in
%   order to be included in the output indices. BP is a number equal to or
%   greater than 1. BP can also contain the string 'full' to indicate that
%   the short reads must be fully contained in the range (or segmented
%   range) in order to be included in the output indices. BP can also
%   contain the string 'start' to search for mapped reads whose start
%   positions lie within the specified range. BP defaults to 1.
%
%   GETINDEX(...,'depth',D) decimates the output indices such that a
%   coverage depth at any base position is equal or less than D. D is a
%   positive integer. Defaults to Inf.
%
%   Examples:
%
%   % Create a BioMap object.
%   s = samread('ex1.sam');
%   obj = BioMap(s);
%
%   % Retrieve the starting positions for all the short reads that
%   % fully overlap the first 50 positions of the mapping.
%   idx =  getIndex(obj,1,50,'overlap','full');
%   x = getStart(obj, idx)
%   y = getStop(obj,idx)
%
%   % Retrieve the sequences for all the short reads that align at least
%   % one base to a segmented range:
%   idx =  getIndex(obj,[98;198],[100;200],'overlap',1);
%   s = getSequence(obj, idx)
%
%   See also BIOMAP, BIOMAP/GET, BIOMAP/GETALIGNMENT, BIOMAP/GETCOUNTS,
%   BIOMAP/GETBASECOVERAGE BIOMAP/GETSTART, BIOMAP/GETSTOP.

%   Copyright 2010-2012 The MathWorks, Inc.

%=== Input check
bioinfochecknargin(nargin, 3, ['BioMap:' mfilename])
checkScalarInput(obj);

Rgiven = false;
if nargin>3 && rem(numel(varargin),2)~=0
    % getIndex(obj,x1,x2,R,varargin) ?
    Rgiven = true;
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
        error(message('bioinfo:BioMap:getIndex:UnknownReference'))
    end
    
    % Calculate the relative index within the requested reference
    if isempty(obj.DictionaryMapping)
        w = find(obj.Index.getField('Reference') == R);
    else
        w = find(obj.DictionaryMapping(obj.Index.getField('Reference')) == R);
    end
        
    obj = getSubset(obj,'SelectReference',R);
    
elseif numel(obj.SequenceDictionary)>1
    error(message('bioinfo:BioMap:getIndex:MultipleReferences'));
end

%=== Error check
if ~isnumeric(x1) || ~isvector(x1) || (isa(x1,'float')&&any(rem(x1,1~=0))) || ...
   ~isnumeric(x2) || ~isvector(x2) || (isa(x2,'float')&&any(rem(x2,1~=0))) || ...
   any(x1(:)<1) || any(x2(:)<x1(:)) || numel(x1)~=numel(x2)
       error(message('bioinfo:BioMap:getIndex:InvalidRange'));
end

%=== Input parsing
[overlap, depth, fullOverlap, startOnly, checkDepth] = parse_inputs(varargin{:});

if numel(x1) == 0
    idx = zeros(0,1);
    return
end

%=== Consolidate ranges (must not overlap and be ordered)
X = [x1(:),x2(:)];
if size(X,1) > 1 
   X = sortrows(X);
   h = zeros(size(X,1),1);
   ls = -inf;
   k = 0;
   for i = 1:size(X,1)
       if X(i,1) > ls+1  % does the current range overlap any other 
                         % range in this segmented range?
           ls = X(i,2);
           k = k+1;
       else
           if X(i,2)>ls
               ls = X(i,2);
           else
               X(i,2) = ls;
           end
       end
       h(i) = k;    
   end
   h = diff(h)>0;
   X = [X([true;h],1) X([h;true],2)];
end
X = uint32(X);

%=== Load start and stop vectors from the BioMap object
start = getStart(obj);
stop = getStop(obj);

%=== Consolidate starts and stops (must be sorted by start position and be 
%    uint32)
if ~isa(start,'uint32')
    start = uint32(start);
end
if ~isa(start,'uint32')
    stop = uint32(stop);
end    

if ~issorted(start)
   [start,sort_idx] = sort(start);
   stop = stop(sort_idx);
else 
   sort_idx = [];
end

if isempty(start)
    idx = zeros(0,1);
else
    idx = bioinfoprivate.getIndexByRangemex(start,stop,X,uint32(overlap),...
                          uint32(depth),fullOverlap,startOnly,checkDepth);
end

if ~isempty(sort_idx)
   idx = sort(sort_idx(idx));
end

if Rgiven
   idx = w(idx); 
end

%--------------------------------------------------------------------------
function [overlap, depth, fullOverlap, startOnly, checkDepth] = parse_inputs(varargin)
% Parse input PV pairs.

%=== defaults
fullOverlap = false;
startOnly = false;
overlap = 1;
depth = 0;
checkDepth = false;

%=== check for the right number of inputs
if rem(nargin,2) == 1
    error(message('bioinfo:BioMap:getIndex:IncorrectNumberOfArguments'))
end

%=== parse parameter value pairs
if nargin > 1
    %=== allowed parameters
    okargs = {'overlap', 'depth'};
	for j = 1:2:nargin
		pname = varargin{j};
		pval = varargin{j+1};
		k = bioinfoprivate.pvpair(pname, pval, okargs, ['BioMap:' mfilename]);
		switch(k)
			case 1  % overlap
                if ischar(pval) && strncmpi(pval(:),'full',numel(pval))
                    fullOverlap = true;
                    startOnly = false;
                    overlap = 1;
                elseif ischar(pval) && strncmpi(pval(:),'start',numel(pval))
                    fullOverlap = false;
                    startOnly = true;
                    overlap = 1;
                elseif isnumeric(pval) && isscalar(pval) && pval>=1
                    fullOverlap = false;
                    startOnly = false;
                    overlap = pval;    
                else
                    error(message('bioinfo:BioMap:getIndex:InvalidOverlap'))
                end
			case 2  % depth
                if isnumeric(pval) && isscalar(pval) && pval>=1
                    depth = pval;
                    checkDepth = true;
                else
                    error(message('bioinfo:BioMap:getIndex:InvalidDepth'))
                end
		end
	end
end
