function [I, data,labels] = genelowvalfilter(data,labels,varargin)
%GENELOWVALFILTER filters genes with low absolute expression levels.
%
%   MASK = GENELOWVALFILTER(DATA) identifies expression profiles in DATA
%   where all absolute expression levels are in the lowest 10 percent of
%   the data set. MASK is a logical vector with one element for each row in
%   DATA. The elements of MASK corresponding to rows with absolute
%   expression levels greater than the threshold have value 1, and those
%   with absolute expression levels less than the threshold are 0. DATA can
%   be a MATLAB numeric matrix or a DataMatrix object.
%
%   [MASK, FDATA] = GENELOWVALFILTER(DATA) returns the filtered data matrix
%   FDATA. FDATA can also be created using FDATA = DATA(MASK,:);  
%
%   [MASK, FDATA, FNAMES] = GENELOWVALFILTER(DATA, NAMES) also returns the
%   filtered names array FNAMES, where NAMES is a cell array of the names
%   of the genes corresponding to each row of DATA. FNAMES can also be
%   created using FNAMES = NAMES(MASK);
%
%   GENELOWVALFILTER(...,'PERCENTILE',PCT) filters genes with absolute
%   expression levels in the lowest PCT percent of the range. 
%
%   GENELOWVALFILTER(...,'ABSVALUE',VAL) filters genes with absolute
%   expression levels lower than VAL. 
%
%   GENELOWVALFILTER(...,'ANYVAL',true) filters genes where any expression
%   level is less than the threshold value.
%
%   Example:
%
%       % Load the yeast workspace variables and filter out the genes
%       % with low absolute expression levels.
%       load yeastdata
%       [mask, fyeastvalues, fgenes] = genelowvalfilter(yeastvalues,genes);
%
%       % Now compare the size of the filtered set with the original set.
%       size (fgenes,1)
%       size (genes,1)
%
%   See also CNSGENEEXPDEMO, EXPRPROFRANGE, EXPRPROFVAR, GENEENTROPYFILTER, GENEVARFILTER,
%   GENERANGEFILTER, YEASTDEMO.

%   Reference:
%     Kohane, I.S., Kho, A.T., Butte, A.J., Microarrays for an Integrative
%     Genomics, MIT Press, Cambridge, MA. 2003.

% Copyright 2003-2008 The MathWorks, Inc.


absPercentile = false;  
absp = 10;
absval = false;
numRows = size(data,1);
anyvals = false;
setval = false;

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
        error(message('bioinfo:genelowvalfilter:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'prctile','absvalue','anyval','percentile'};
    for j=1:2:numArgs-2
        [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
        
        switch(k)
            case {1,4}  % absPercentile
                absPercentile = true;
                setval = true;
                absp = pval;
            case 2  % absval
                setval = true;
                absval = true;
                absv = pval;
            case 3
                anyvals = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        end
    end
end

% if no options specified, use prctile
if ~setval
    absPercentile = true;
end
%== 
if isa(data, 'bioma.data.DataMatrix')
    absdata = abs(data.(':')(':'));
else
    absdata = abs(data);
end
% all or any
if anyvals
    lowvals = min(absdata,[],2);
else
    lowvals = max(absdata,[],2);
end

% create the index
I = true(numRows,1);

if absPercentile
    abscutoff = prctile(absdata(:),absp);
    I = (lowvals>abscutoff);
end

if absval 
    I = I & (lowvals>absv);
end

% handle labels if they were specified
if ~isempty(labels) && nargout > 2
    if ischar(labels)
        labels = cellstr(labels);
    end
    if numel(labels) ~= numRows
        warning(message('bioinfo:genelowvalfilter:LabelSizeMismatch'));
        labels = [];
    else
        labels = labels(I);
    end
end

if nargout > 1
    data = data(I,:);
end
