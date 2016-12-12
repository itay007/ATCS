function [absrange,logrange] = exprprofrange(data,varargin)
%EXPRPROFRANGE calculates the range of expression profiles.
%
%   EXPRPROFRANGE(DATA) calculates the range of each expression profile in
%   dataset DATA. DATA can be a MATLAB numeric array or a DataMatrix
%   object.
%
%   [RANGE, LOGRANGE] = EXPRPROFRANGE(DATA) also calculates the log range,
%   that is log(max(prof)) - log(min(prof)), of each expression profile. If
%   DATA is a DataMatrix object, outputs RANGE and LOGRANGE are also
%   DataMatrix objects with the same row names as DATA and an empty column
%   name.
%
%   If no output arguments are specified, a histogram bar plot of the range
%   is displayed.   
%
%   EXPRPROFRANGE(...,'SHOWHIST',TF) displays a histogram of the range
%   data if TF is true.
%
%   Example:
%
%       load yeastdata
%       range = exprprofrange(yeastvalues,'showhist',true);
%
%   See also EXPRPROFVAR, GENERANGEFILTER.

% Copyright 2003-2008 The MathWorks, Inc.


showhist = false;
if nargout == 0
    showhist = true;
end
dorel = false;

if nargout >1
    dorel = true;
end

if nargin > 1
    if rem(nargin,2)== 0
        error(message('bioinfo:exprprofrange:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'showhist'};
    for j=1:2:nargin-2
        [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
        
        switch(k)
            case 1  % showhist
                showhist = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        end
    end
end

lb = min(data,[],2);
ub = max(data,[],2);
absrange = ub-lb;

% log range only makes sense if we have absolute values -- guess if this is
% the case by looking for all non-negative values.
if dorel
    if any(lb<0)
        warning(message('bioinfo:exprprofrange:CannotCalculateLogRange'))
        logrange = [];
    else
        logrange = log(ub) - log(lb);
    end
end
%== Output DataMatrix
if nargout > 0 && isa(data, 'bioma.data.DataMatrix') 
    absrange = bioma.data.DataMatrix(absrange, data.RowNames, []);
    if nargout > 1
        if isempty(logrange)
            logrange = bioma.data.DataMatrix;
        else
            logrange = bioma.data.DataMatrix(logrange, data.RowNames, []);
        end
    end
end

% show histogram -- if both abs and rel range make sense, show these on the
% same plot in two axes.
if showhist
    if isa(data, 'bioma.data.DataMatrix')
        numbuckets = max(10,numel(absrange, ':', ':')/100);
    else
        numbuckets = max(10,numel(absrange)/100);
    end
    numbuckets = ceil(min(numbuckets,100));
    if dorel
        subplot(2,1,1);
    end
    hist(absrange,numbuckets);
    title('Profile Ranges');
    if dorel
        subplot(2,1,2);
        title('Profile Log Ranges');
        hist(logrange(isfinite(logrange), :),numbuckets);
    end
end
