function [I,data,labels] = genevarfilter(data,labels,varargin)
%GENEVARFILTER filters genes with small profile variance.
%
%   MASK = GENEVARFILTER(DATA) identifies expression profiles in DATA with
%   variance in the lowest 10 percent. MASK is a logical vector with one
%   element for each row in DATA. The elements of MASK corresponding to
%   rows with variance greater than the threshold have value 1, and those
%   with variance less than the threshold are 0. DATA can be a MATLAB
%   numeric matrix or a DataMatrix object.
%
%   [MASK, FDATA] = GENEVARFILTER(DATA) returns the filtered data matrix
%   FDATA. FDATA can also be created using FDATA = DATA(MASK,:);  
%
%   [MASK, FDATA, FNAMES] = GENEVARFILTER(DATA, NAMES) also returns the
%   filtered names array FNAMES, where NAMES is a cell array of the names
%   of the genes corresponding to each row of DATA. FNAMES can also be
%   created using FNAMES = NAMES(MASK);
%
%   GENEVARFILTER(...,'PERCENTILE',PCT) filters genes with profile variances
%   in the lowest PCT percent of the range. 
%
%   GENEVARFILTER(...,'ABSVALUE',VAL) filters genes with profile variances
%   lower than VAL. 
%
%   Example:
%
%       % Load the yeast workspace variables and filter out the genes
%       % with a small profile variance.
%       load yeastdata
%       [mask, fyeastvalues, fgenes] = genevarfilter(yeastvalues,genes);
%
%       % Now compare the size of the filtered set with the original set.
%       size (fgenes,1)
%       size (genes,1)
%
%   See also CNSGENEEXPDEMO, EXPRPROFRANGE, EXPRPROFVAR, GENEENTROPYFILTER,
%   GENELOWVALFILTER, GENERANGEFILTER, YEASTDEMO.

%   Reference:
%     Kohane, I.S., Kho, A.T., Butte, A.J., Microarrays for an Integrative
%     Genomics, MIT Press, Cambridge, MA. 2003.


% Copyright 2003-2008 The MathWorks, Inc.


absPercentile = true;
absp = 10;
absval = false;
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
        error(message('bioinfo:genevarfilter:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'prctile','absvalue','percentile'};
    for j=1:2:numArgs-2
        [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
        
        switch(k)
            case {1,3}  % absPercentile
                absPercentile = true;
                absp = pval;
            case 2  % absval
                absval = true;
                absv = pval;
        end
    end
end

absvar = var(data,[],2);

% create the index
I = true(numRows,1);

if absPercentile
    abscutoff = prctile(absvar,absp);
    I = absvar>abscutoff;
end

if absval 
    I = I & (absvar>absv);
end

% handle labels if defined
if ~isempty(labels) && nargout > 2
    if ischar(labels)
        labels = cellstr(labels);
    end
    if numel(labels) ~= numRows
        warning(message('bioinfo:genevarfilter:LabelSizeMismatch'));
        labels = [];
    else
        labels = labels(I);
    end
end

if nargout > 1
    data = data(I,:);
end

