function [I,data,labels] = generangefilter(data,labels,varargin)
%GENERANGEFILTER filters genes with small profile ranges.
%
%   MASK = GENERANGEFILTER(DATA) identifies expression profiles in DATA
%   with range in the lowest 10 percent. MASK is a logical vector with one
%   element for each row in DATA. The elements of MASK corresponding to
%   rows with range greater than the threshold have value 1, and those with
%   range less than the threshold are 0. DATA can be a MATLAB numeric
%   matrix or a DataMatrix object.
%
%   [MASK, FDATA] = GENERANGEFILTER(DATA) returns the filtered data matrix
%   FDATA. FDATA can also be created using FDATA = DATA(MASK,:);  
%
%   [MASK, FDATA, FNAMES] = GENERANGEFILTER(DATA, NAMES) also returns the
%   filtered names array FNAMES, where NAMES is a cell array of the names
%   of the genes corresponding to each row of DATA. FNAMES can also be
%   created using FNAMES = NAMES(MASK);
%
%   GENERANGEFILTER(...,'ABSVALUE',VAL) filters genes with profile ranges
%   less than VAL. 
%
%   GENERANGEFILTER(...,'PERCENTILE',PCT) filters genes with profile ranges
%   in the lowest PCT percent of the range. 
%
%   GENERANGEFILTER(...,'LOGVALUE',VAL) filters genes with profile log
%   ranges less than VAL.
%
%   GENERANGEFILTER(...,'LOGPERCENTILE',PCT) filters genes with profile
%   ranges in the lowest PCT percent of the log range.
%
%   Example:
%
%       % Load the yeast workspace variables and filter out the genes
%       % with small profile ranges.
%       load yeastdata
%       [mask, fyeastvalues, fgenes] = generangefilter(yeastvalues,genes);
%
%       % Now compare the size of the filtered set with the original set.
%       size (fgenes,1)
%       size (genes,1)
%
%   See also EXPRPROFRANGE, EXPRPROFVAR, GENEENTROPYFILTER,
%   GENELOWVALFILTER, GENEVARFILTER.

%   Reference:
%     Kohane, I.S., Kho, A.T., Butte, A.J., Microarrays for an Integrative
%     Genomics, MIT Press, Cambridge, MA. 2003.

% Copyright 2003-2008 The MathWorks, Inc.


absPercentile = true;
absp = 10;
dorel = false;
absval = false;
logval = false;
logPercentile = false;
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
    absPercentile = false;  
    if rem(numArgs,2)== 1
        error(message('bioinfo:generangefilter:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'prctile','absvalue','logprctile','logvalue','percentile','logpercentile'};
    for j=1:2:numArgs-2
        [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);

        switch(k)
            case {1,5}  % absPercentile
                absPercentile = true;
                absp = pval;
            case 2  % absval
                absval = true;
                absv = pval;
            case {3,6}  % logPercentile
                dorel = true;
                logPercentile = true;
                logp = pval;
            case 4  % logval
                dorel = true;
                logval = true;
                logv = pval;
        end
    end
end

if dorel
    [absrange, logrange] = exprprofrange(data,'showhist',false);
    if isempty(logrange)
        dorel = false;
    end
    
else
    absrange = exprprofrange(data,'showhist',false);
end

I = true(numRows,1);

if absPercentile
    abscutoff = prctile(absrange,absp);
    I = (absrange>abscutoff);
end

if absval 
    I = I & (absrange>absv);
end

if logPercentile && dorel
    logcutoff = prctile(logrange,logp);
    I = I & (logrange>logcutoff);
end

if logval && dorel
    I = I & (logrange>logv);
end

% handle labels if defined
if ~isempty(labels) && nargout > 2
    if ischar(labels)
        labels = cellstr(labels);
    end
    if numel(labels) ~= numRows
        warning(message('bioinfo:generangefilter:LabelSizeMismatch'));
        labels = [];
    else
        labels = labels(I);
    end
end

if nargout > 1
    data = data(I,:);
end



