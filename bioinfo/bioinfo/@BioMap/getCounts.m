function count = getCounts(obj, x1, x2, varargin)
%GETCOUNTS return the raw count of mapped reads for a set of ranges
%
%   COUNT = GETCOUNTS(OBJ, X1, X2) returns a scalar with the number of
%   short reads that align to the region specified by X1 and X2. X1 and X2
%   are two non-negative integers such that X1 <= X2 and both are smaller
%   than the length of the reference sequence. X1 and X2 may also be
%   vectors of the same length representing a segmented range.
%
%   By default GETCOUNTS considers each short read once; any two or more
%   overlapping segments in X1 and X2 are merged into one, and, short reads
%   that span over two different segments are counted once. 
%
%   GETCOUNTS is equivalent to:
%         COUNT = numel(getIndex(OBJ,X1,X2));
%
%   COUNTS = GETCOUNTS(...,'independent',true) treats the segments in X1
%   and X2 independently. In this case GETCOUNTS returns a column vector
%   COUNTS with the same number of elements as X1 and X2. Default is false.
%
%   GETCOUNTS with the 'independent' option set to true is equivalent to: 
%         for i = 1:size(X,1)
%            COUNTS(i) = getCounts(OBJ,X1(i),X2(i));
%         end
%
%   GROUP_COUNTS = GETCOUNTS(OBJ, X1, X2, S) specifies a grouping variable
%   indicating which of different segmented ranges each segment in X1 and
%   X2 belongs to. Segmented ranges are treated independently. The output
%   GROUP_COUNTS is a column vector with the same number of elements as
%   unique values in S. The order of the elements in GROUP_COUNTS
%   corresponds to the ascending order of unique elements in S.
%
%   GETCOUNTS(OBJ, X1, X2, S) is equivalent to:
%         [~,h] = unique(S);
%         for i = 1:numel(h) 
%             g = S(h(i))==S;
%             GROUP_COUNTS(i) = getCounts(OBJ,X1(g),X2(g));
%         end
%
%   Also, GETCOUNTS(OBJ,X1,X2) is equivalent to
%   GETCOUNTS(OBJ,X1,X2,ones(numel(X1),1)) 
%
%   For example GETCOUNTS(OBJ,[1;30;50],[10;60;60],[2;1;2]) returns a 2x1
%   vector of counts where GROUP_COUNTS(1) contains the number of reads
%   that align to the range 30:60 and GROUP_COUNTS(2) contains the number
%   of reads that align to the segmented range [1:10 50:60]. 
%
%   GETCOUNTS(OBJ, X1, X2, S, R) specifies a reference for each one of the
%   segmented ranges defined by X1, X2 and S respectively. R can be a
%   vector with indices to the SequenceDictionary property or a cell string
%   with the actual name of the reference. R must be ordered and have the
%   same numebr of elements as the unique elements in S. R can also have
%   the same size as S, in such case, all the entries in R for each unique
%   value in S must be the same.
%
%   GETCOUNTS(...,'overlap',BP) specifies the minimum number of positions
%   that a short read must overlap the given range (or segmented range) in
%   order to be counted. BP is a number equal or greater than 1. BP can also 
%   contain the string 'full' to indicate that the short reads must be fully
%   contained in the range (or segmented range) in order to be counted. BP
%   can also contain the string 'start' to search for mapped reads whose
%   start positions lie within the specified range. BP defaults to 1.
%
%   GETCOUNTS(...,'method',M) changes the method used to measure the
%   abundance of reads for a given range (or segmented range), options are:
%
%         'raw'  - raw counts (default)
%         'rpkm' - reads per kilo-bp per million mapped reads
%         'mean' - average coverage depth computed base-by-base
%         'max'  - maximum coverage depth computed base-by-base
%         'min'  - minimum coverage depth computed base-by-base
%         'sum'  - sum of all aligned bases
%
%
%   Examples:
%
%   % Create a BioMap object.
%   s = samread('ex1.sam');
%   obj = BioMap(s);
%
%   % Get the number of reads that cover at least on base of the segmented 
%   % range [1:50 71:100]:
%
%   counts_1 = getCounts(obj,[1;71],[50;100])
% 
%   % Get the number of reads that cover at least one base of two independent 
%   % segments; [1:50] and [71:100]:
%
%   counts_2 = getCounts(obj,[1;71],[50;100],'independent',true)
%
%   % Observe that sum(counts_2)>counts_1 because there are four reads that
%   % span over the two segments and are counted twice in the second case.
%
%   % Returns a 2x1 vector of counts where counts_3(1) contains the number
%   % of reads that align to the range 30:60 and counts_3(2) contains the
%   % number of reads that align to the segmented range [1:10 50:60]:
%
%   counts_3 = getCounts(obj,[1;30;50],[10;60;60],[2 1 2])
%
%   % Return the number of mapped reads:
%
%    getCounts(obj,min(getStart(obj)),max(getStop(obj)))
%   
%   See also BIOMAP, BIOMAP/GET, BIOMAP/GETALIGNMENT, BIOMAP/GETBASECOVERAGE, 
%   BIOMAP/GETINDEX, BIOMAP/GETSTART, BIOMAP/GETSTOP. 

%   Copyright 2010-2012 The MathWorks, Inc.


%=== Input check
bioinfochecknargin(nargin, 3, ['BioMap:' mfilename])
checkScalarInput(obj);

%=== Error check
if ~isnumeric(x1) || ~isvector(x1) || (isa(x1,'float')&&any(rem(x1,1~=0))) || ...
   ~isnumeric(x2) || ~isvector(x2) || (isa(x2,'float')&&any(rem(x2,1~=0))) || ...
   any(x1(:)<1) || any(x2(:)<x1(:)) || numel(x1)~=numel(x2)
    error(message('bioinfo:BioMap:getCounts:InvalidRange'));
end

X = [x1(:),x2(:)];

if size(X,1) == 0
    count = 0;
    return
end

%=== Find out if the GETCOUNTS(OBJ,X1,X2,S,...) signature is being used
if (nargin>3) 
    if (isnumeric(varargin{1})&&isvector(varargin{1})) || iscellstr(varargin{1}) 
        S = varargin{1};
        varargin = varargin(2:end);
        Sgiven = true;
    else
        Sgiven = false;
    end
else
    Sgiven = false;
end

%=== Find out if the GETCOUNTS(OBJ,X1,X2,S,R,...) signature is being used
if Sgiven && (nargin>4)
    if (isnumeric(varargin{1})&&isvector(varargin{1})) || iscellstr(varargin{1}) 
        R = varargin{1};
        varargin = varargin(2:end);
        Rgiven = true;
    else
        Rgiven = false;
    end
else
    Rgiven = false;
end

%=== When multiple references are requested, the input object is subsetted
%    for each reference and getCounts is called with each subset:
if Rgiven
    % validate R
    if iscellstr(R) % change to numeric
        h = zeros(size(R));
        for i = 1:numel(obj.SequenceDictionary)
             k = strcmp(R,obj.SequenceDictionary{i});    
             if any(k)
                 h(k) = i;
             end
        end
        R = h(:);
    end

    if any(R<1) || any(R>numel(obj.SequenceDictionary)) || any(rem(R,1))
        error(message('bioinfo:BioMap:getCounts:UnknownReference'))
    end

    [SU,~,SUi] = unique(S);
    
    if (numel(R) ~= numel(SU)) && (numel(R)~=1) && (numel(R)~=numel(S))
        error(message('bioinfo:BioMap:getCounts:InvalidSizeReference'))
    end
    
    if (numel(R) ~= numel(SU)) && (numel(R)~=1)
        % check that for each SU there is only one R
        for i = 1:numel(SU)
            if numel(unique(R(SUi==i)))>1
                error(message('bioinfo:BioMap:getCounts:MultipleReferencesInSegment'))
            end
        end
    end
    
    RU = unique(R);
    
    if numel(RU) == 1
        count = getCounts(getSubset(obj,'SelectReference',RU),X(:,1),X(:,2),S,varargin{:});
    else
       count = zeros(numel(SU),1);
       for i=1:numel(RU)
           h = R==RU(i);
           count(h) = getCounts(getSubset(obj,'SelectReference',RU(i)),X(h,1),X(h,2),S(h),varargin{:});
       end
    end

    return
end

if numel(obj.SequenceDictionary)>1
    error(message('bioinfo:BioMap:getCounts:MultipleReferences'));
end

%=== Input parsing
[overlap, fullOverlap, startOnly, independent, method] = parse_inputs(Sgiven,varargin{:});

if Sgiven
    if numel(S)~=size(X,1)
         error(message('bioinfo:BioMap:getCounts:InvalidGroupVariable'));
    end
    % obtain unique labels (this is the order of the output)
    % S can be a numeric vector or a cellstr
    [~,~,w] = unique(S(:));
    % order groups by leftmost position, then by rightmost position
    v1 = accumarray(w(:),X(:,1),[],@min);
    v2 = accumarray(w(:),X(:,2),[],@max);
    [~,u] = sortrows([v1,v2]);
    v(u) = 1:numel(u);
    X(:,3) = v(w);
    % order segments, first by group label then by leftmost position
    X = sortrows(X,[3 1 2]);
  
else
    if independent
        % order segments by leftmost position
        [X,h] = sortrows(X);
        X(:,3) = (1:size(X,1))';
        v(h) = X(:,3); % for re-ordering counts to original order,
                       % getCountsByRangemex requires indepndent groups 
                       % to be ordered by leftmost position.
    else
        X = sortrows(X);
        X(:,3) = ones(size(X,1),1);
        v = 1;
    end
end

%=== Consolidate ranges (non independent-segments must not overlap and 
%    all segments must be ordered first by group and then by position) 
if size(X,1) > 1 
   h = zeros(size(X,1),1);
   lg = 0;
   k = 0;
   for i = 1:size(X,1)
       if X(i,3)~=lg        % is it a new segmented range?
           ls = X(i,2);
           lg = X(i,3);
           k = k+1;
       elseif X(i,1) > ls+1 % does the current range overlap any other 
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
   X = [X([true;h],1) X([h;true],[2 3])];
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
end 

if isempty(start)
    count = zeros(numel(unique(X(:,3))),1);
else
    count = bioinfoprivate.getCountsByRangemex(start,stop,X,...
              uint32(overlap),fullOverlap,startOnly,method);
end

%=== Reorder to the same order as unique group labels
count = count(v);

%=== Convert rpk to rpkm          
if strcmp(method,'rpkm')
   numMappedReads = double(sum(start~=stop));
   count = count .* (1000000./numMappedReads);
end

%--------------------------------------------------------------------------
function [overlap, fullOverlap, startOnly, independent, method] = parse_inputs(Sgiven,varargin)
% Parse input PV pairs.

%=== defaults
fullOverlap = false;
startOnly = false;
overlap = 1;
if Sgiven
    independent = true;
else
    independent = false;
end
method = 'raw';

%=== check for the right number of inputs
if rem(nargin-1,2) == 1
    error(message('bioinfo:BioMap:getCounts:IncorrectNumberOfArgumentsOrInvalidGroupVariable'))
end

%=== parse parameter value pairs
if nargin-1 > 1
    %=== allowed parameters
    okargs = {'overlap', 'independent', 'method'};
	for j = 1:2:nargin-1
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
                    error(message('bioinfo:BioMap:getCounts:InvalidOverlap'))
                end
            case 2 % independent
                if Sgiven
                    warning(message('bioinfo:BioMap:getCounts:CannotSetIndependent'))
                else
                    independent = bioinfoprivate.opttf(pval, okargs{k}, ['BioMap:' mfilename]); 
                end
            case 3 % method
                [~,method] = bioinfoprivate.optPartialMatch(pval,{'raw','rpkm','max','min','mean','sum'}, okargs{k}, ['BioMap:' mfilename]);
		end
	end
end
