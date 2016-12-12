function [coverage, bins] = getBaseCoverage(obj, x1, x2, varargin)
%GETBASECOVERAGE compute the base-by-base coverage of a BioMap object.
%
%   COV = GETBASECOVERAGE(OBJ,X1,X2) computes the base-by-base coverage of
%   the alignment in the region specified by X1 and X2. X1 and X2 are two
%   non-negative integers such that X1 <= X2 and both are smaller than the
%   length of the reference sequence.
%
%   X1 and X2 may also be vectors of the same length representing a
%   segmented range. When a segmented range is given, the length of COV is
%   numel(min(X1):max(X2)) and COV contains a NaN for the positions
%   between the segments.
%
%   COV = GETBASECOVERAGE(OBJ,X1,X2,R) select the reference where the
%   coverage is calculated. R is either an index to one of the references
%   in the SequenceDictionary property or is a string with the actual
%   reference name.
%
%   GETBASECOVERAGE(...,'binWidth',W) decimates the output using a binning
%   algorithm with bins of size W (in bp). Bins are centered within min(X1)
%   and max(X2), i.e. the first and the last bin span approximately equally
%   outside the [min(X1) max(X2)] range.
%
%   GETBASECOVERAGE(...,'numberOfBins',N) decimates the output using a
%   binning algorithm with N bins of equal size in order to span the
%   requested region. Bins are centered within min(X1) and max(X2), i.e.
%   the first and the last bin span approximately equally outside the
%   [min(X1) max(X2)] range.
%
%   GETBASECOVERAGE(...,'complementRanges',true) computes the coverage for
%   the positions between the segments. The length of COV is
%   numel(min(X1):max(X2)) and COV contains a NaN for the positions
%   inside the segments denoted by X1 and X2.
%
%   GETBASECOVERAGE(...,'binType',T) sets the type of binning algorithm to
%   'max', 'min', or 'mean'. Default is 'max'.
%
%   Note: 'binWidth' and 'numberOfBins' cannot be used at the same time.
%
%   [COV BIN] = GETBASECOVERAGE(...) returns the start position of every
%   bin in the output vector BIN. BIN has the same size as COV. If no
%   binning occurs BIN = min(X1):max(X2).
%
%   Examples:
%
%   % Create a BioMap object.
%   s = samread('ex1.sam');
%   obj = BioMap(s);
%
%   % Compute the coverage for the first 50 positions of the mapping.
%   cov = getBaseCoverage(obj, 1, 50)
%
%   % Compute the coverage for the first 2000 positions decimating the
%   % output vector by a factor of 100.
%    [cov,bin_starts] = getBaseCoverage(obj, 1, 2000, 'binWidth', 100)
%
%   See also BIOMAP, BIOMAP/GET, BIOMAP/GETALIGNMENT, BIOMAP/GETCOUNTS,
%   BIOMAP/GETINDEX, BIOMAP/GETSTART, BIOMAP/GETSTOP.

%   Copyright 2010-2012 The MathWorks, Inc.

%=== Input check
bioinfochecknargin(nargin, 3, ['BioMap:' mfilename])
checkScalarInput(obj);

if nargin>3 && rem(numel(varargin),2)~=0
    % getBaseCoverage(obj,x1,x2,R,varargin) ?
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
        error(message('bioinfo:BioMap:getBaseCoverage:UnknownReference'))
    end
    obj = getSubset(obj,'SelectReference',R);
elseif numel(obj.SequenceDictionary)>1
    error(message('bioinfo:BioMap:getBaseCoverage:MultipleReferences'));
end
   
%=== Error check
if ~isnumeric(x1) || ~isvector(x1) || (isa(x1,'float')&&any(rem(x1,1~=0))) || ...
        ~isnumeric(x2) || ~isvector(x2) || (isa(x2,'float')&&any(rem(x2,1~=0))) || ...
        any(x1(:)<1) || any(x2(:)<x1(:)) || numel(x1)~=numel(x2) || any(isnan(x1)) || ...
        any(isnan(x2)) || any(isinf(x1)) || any(isinf(x2))
    error(message('bioinfo:BioMap:getBaseCoverage:InvalidRange'));
end

X = [x1(:),x2(:)];

if size(X,1) == 0
    coverage = [];
    bins = [];
    return
end

%=== Consolidate ranges (must not overlap and be ordered)
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

Xstart = X(1);
Xend = X(end);

[binType, binWidth, numberOfBins, complementRanges] = parse_inputs(double(Xend-Xstart+1),varargin{:});

% Figure out start of and width of bins
if ~isnan(binWidth)
    binRange = ceil( double(Xend-Xstart+1) ./ binWidth ) .* binWidth;
    numberOfBins = binRange ./ binWidth;
elseif ~isnan(numberOfBins)
    binRange = ceil(double(Xend-Xstart+1) ./ numberOfBins) .* numberOfBins;
    binWidth = binRange ./ numberOfBins;
else
    binRange = double(Xend-Xstart+1);
    binWidth = 1;
    numberOfBins =  binRange;
end
leftBinExcess = floor((binRange - double(Xend-Xstart+1))./2);

if complementRanges
    X = [X(1:end-1,2)+1 X(2:end,1)-1];
    if isempty(X)
        coverage = nan(1,numberOfBins); 
        % Calculates left position of the bins
        if nargout>1
           bins = (double(Xstart)-leftBinExcess-1)+(1:binWidth:binRange);
        end
        return
    end
end

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

% Call mex implementation for computing the coverage
if isempty(start)
    coverage = zeros(1,numberOfBins);
else
    coverage = bioinfoprivate.getCoverageBaseByBasemex(start,stop,...
        X,Xstart,Xend,uint32(leftBinExcess),...
        uint32(binWidth),uint32(numberOfBins),binType);
end

% Calculates left position of the bins
if nargout>1
   bins = (double(Xstart)-leftBinExcess-1)+(1:binWidth:binRange);
end

%--------------------------------------------------------------------------
function [binType, binWidth, numberOfBins, complementRanges] = parse_inputs(Xr,varargin)
% Parse input PV pairs.

%=== defaults
binType = 'max';
binWidth = NaN;
numberOfBins = NaN;
complementRanges = false;

%=== check for the right number of inputs
if rem(numel(varargin),2) == 1
    error(message('bioinfo:BioMap:getBaseCoverage:IncorrectNumberOfArguments'))
end

%=== parse parameter value pairs
if numel(varargin) > 1
    %=== allowed parameters
    okargs = {'bintype', 'binwidth', 'numberofbins', 'complementranges'};
    for j = 1:2:numel(varargin)
        pname = varargin{j};
        pval = varargin{j+1};
        k = bioinfoprivate.pvpair(pname, pval, okargs, ['BioMap:' mfilename]);
        switch(k)
            case 1  % binType
                [~,binType] = bioinfoprivate.optPartialMatch(pval,{'max','min','mean'}, okargs{k}, ['BioMap:' mfilename]);
            case 2  % binWidth
                if isnumeric(pval) && isscalar(pval) && pval>=1 && (pval>Xr || ~rem(pval,1))
                    binWidth = double(pval);
                    if binWidth>=Xr
                        warning(message('bioinfo:BioMap:getBaseCoverage:TooLargeBinWidth'))
                        binWidth = Xr;
                    end
                else
                    error(message('bioinfo:BioMap:getBaseCoverage:InvalidBinWidth'))
                end
            case 3  % numberOfBins
                if isnumeric(pval) && isscalar(pval) && pval>=1 && (pval>Xr || ~rem(pval,1))
                    numberOfBins = double(pval);
                    if numberOfBins>=Xr
                        warning(message('bioinfo:BioMap:getBaseCoverage:TooLargeNumberOfBins'))
                        numberOfBins = Xr;
                    end
                else
                    error(message('bioinfo:BioMap:getBaseCoverage:InvalidNumberOfBins'))
                end
            case 4  % complementRanges
                complementRanges = bioinfoprivate.opttf(pval, okargs{k}, ['BioMap:' mfilename]);
        end
    end
end

if ~isnan(numberOfBins) && ~isnan(binWidth)
    error(message('bioinfo:BioMap:getBaseCoverage:SimultaneousBinWidthNumberOfBins'))
end




