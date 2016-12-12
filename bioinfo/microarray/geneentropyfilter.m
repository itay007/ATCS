function [I data labels] = geneentropyfilter(data,labels,varargin)
%GENEENTROPYFILTER filters genes with profiles with low entropy.
%
%   MASK = GENEENTROPYFILTER(DATA) identifies expression profiles in DATA
%   with entropy in the lowest 10 percent. MASK is a logical vector with
%   one element for each row in DATA. The elements of MASK corresponding to
%   rows with variance greater than the threshold have value 1, and those
%   with variance less than the threshold are 0. DATA can be a MATLAB
%   numeric matrix or a DataMatrix object.
%
%   [MASK, FDATA] = GENEENTROPYFILTER(DATA) returns the filtered data
%   matrix FDATA. FDATA can also be created using FDATA = DATA(MASK,:);  
%
%   [MASK, FDATA, FNAMES] = GENEENTROPYFILTER(DATA, NAMES) also returns the
%   filtered names array FNAMES, where NAMES is a cell array of the names
%   of the genes corresponding to each row of DATA. FNAMES can also be
%   created using FNAMES = NAMES(MASK);
%
%   GENEENTROPYFILTER(..., 'PERCENTILE', PCT) filters genes with entropy
%   levels in the lowest PCT percent of the data.
%
%   Example:
%
%       % Load the yeast workspace variables and filter out the genes
%       % with low entropy.
%       load yeastdata
%       [mask, fyeastvalues, fgenes] = geneentropyfilter(yeastvalues,genes);
%
%       % Now compare the size of the filtered set with the original set.
%       size (fgenes,1)
%       size (genes,1)
%
%   See also EXPRPROFRANGE, EXPRPROFVAR, GENELOWVALFILTER, GENEVARFILTER,
%   GENERANGEFILTER, YEASTDEMO.

%   Reference:
%     Kohane, I.S., Kho, A.T., Butte, A.J., Microarrays for an Integrative
%     Genomics, MIT Press, Cambridge, MA. 2003.

% Copyright 2003-2008 The MathWorks, Inc.


absp = 10;
numBins = ceil(size(data,2)/2);
numRows = size(data,1);

numArgs = nargin;
if numArgs < 2
    labels = [];
end

if numArgs > 1 && ischar(labels) && isvector(labels)
        varargin = {labels, varargin{:}};
        labels = [];
        numArgs = numArgs +1;
end

if numArgs > 2
    if rem(numArgs,2)== 1
        error(message('bioinfo:geneentropyfilter:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'prctile','bins','percentile'};
    for j=1:2:numArgs-2
        [k,pentropy] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
        switch(k)
            case {1,3}  % absprctile
                absp = pentropy;
            case 2  % numBins
                numBins = pentropy;
        end
    end
end 

p = zeros(numRows,numBins);

if isa(data, 'bioma.data.DataMatrix')
    for count = 1:numRows
        p(count,:) = hist(data.(count)(':'), numBins);
        p(count,:) = p(count,:)./sum(p(count,:)); %normalize hist to get probability
    end
else
    for count = 1:numRows
        p(count,:) = hist(data(count,:),numBins);
        p(count,:) = p(count,:)./sum(p(count,:)); %normalize hist to get probability
    end
end

% calculate sum(p.*log2(p))
% this will run into problems with 0*log2(0) so we trick this into doing
% the right thing.

p(p==0) = 1;
entropy = -sum(p.*log2(p),2);

% create the index
abscutoff = prctile(entropy,absp);
I = (entropy>abscutoff);

% handle labels if they were specified
if ~isempty(labels) && nargout > 2
    if ischar(labels)
        labels = cellstr(labels);
    end
    if numel(labels) ~= numRows
        warning(message('bioinfo:geneentropyfilter:LabelSizeMismatch'));
        labels = [];
    else
        labels = labels(I);
    end
end
if nargout > 1
    data = data(I,:);
end

