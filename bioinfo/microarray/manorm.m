function [X,gnorm] = manorm(X,varargin)
%MANORM normalizes microarray data
%
%   XNORM = MANORM(X), where X is a vector or a matrix, scales the values
%   in each column of X by dividing by the mean column intensity. X can be
%   a DataMatrix object. If X is a DataMatrix object, XNORM is a
%   DataMatrix object with the same properties as X.
%
%   XNORM = MANORM(MASTRUCT,FIELDNAME), where MASTRUCT is a microarray
%   structure, scales the data for field FIELDNAME for each block, or
%   print-tip, by dividing by the mean column intensity for each block. The
%   output is a matrix with each column corresponding to the normalized
%   data for each block.
%
%   [XNORM, COLVAL] = MANORM(...) returns the values used to normalize the
%   data.
%
%   MANORM(...,'METHOD',METHOD) allows you to choose the method used for
%   scaling or centering the data. METHOD can be one of 'Mean' (default),
%   'Median', 'STD' (standard deviation), 'MAD' (median absolute
%   deviation), or a function handle. If you pass a function handle then
%   the function should ignore NaNs and must return a single value per
%   column of the input data.
%
%   MANORM(...,'EXTRA_ARGS',ARGS) allows you to pass extra arguments to the
%   'METHOD' function. ARGS must be a cell array.
%
%   MANORM(...,'LOGDATA',TRUE) is used for working with log ratio data in
%   which case the mean (or METHOD value) of each column is subtracted
%   from the values in the columns, instead of dividing the column by the
%   normalizing value.
%
%   MANORM(...,'PERCENTILE',PCT) only uses the PCT percentile of the data,
%   preventing large outliers from skewing the normalization. If PCT is a
%   vector containing two values, then the range from the PCT(1) percentile
%   to the PCT(2) percentile is used. The default value is 100, that is, to
%   use all the data in the data set.
%
%   MANORM(...,'GLOBAL',TRUE) normalizes the values in the data set by the
%   global mean (or METHOD value) of the data, as opposed to normalizing
%   each column, or block, of the data independently. 
%
%   MANORM(...,'STRUCTOUTPUT',true), when the input data is a structure,
%   returns the input structure with an additional data field for the
%   normalized data.
%
%   MANORM(...,'NEWCOLUMNNAME',COLNAME), when using STRUCTOUTPUT, allows
%   you to specify the name of the column that is appended to the list of
%   ColumnNames in the structure. The default behavior is to prefix 'Block
%   Normalized' to the FIELDNAME string.
%
%   Examples:
%
%       maStruct = gprread('mouse_a1wt.gpr');
%       % Extract some data of interest.
%       Red = magetfield(maStruct,'F635 Median');
%       Green = magetfield(maStruct,'F532 Median');
%       % Create a log-log plot.
%       maloglog(Red,Green,'factorlines',true)
%       % Center the data.
%       normRed = manorm(Red);
%       normGreen = manorm(Green);
%       % Create a log-log plot of the centered data.
%       figure
%       maloglog(normRed,normGreen,'title','Normalized','factorlines',true)
%
%       % Alternatively, you can work directly with the structure
%       normRedBs = manorm(maStruct,'F635 Median - B635');
%       normGreenBs = manorm(maStruct,'F532 Median - B532');
%       % Create a log-log plot of the centered data. This includes some
%       % zero values so turn off the warning.
%       figure
%       w = warning('off','bioinfo:maloglog:ZeroValues');
%       warning('off','bioinfo:maloglog:NegativeValues');
%       maloglog(normRedBs,normGreenBs,'title',...
%               'Normalized Background-Subtracted Median Values',...
%               'factorlines',true)
%       warning(w);
%
%   See also AFFYINVARSETNORM, MABOXPLOT, MAGETFIELD, MAINVARSETNORM, MAIRPLOT,
%   MALOGLOG, MALOWESS, MOUSEDEMO, QUANTILENORM, RMASUMMARY. 

% Reference:
%       Stekel, B., Microarray Bioinformatics  
%       Cambridge University Press 2003.

% Copyright 2004-2008 The MathWorks, Inc.


% Note that the structure handling is dealt with by @struct/manorm.m.

%==Check inputs
import bioinfoprivate.*;
bioinfochecknargin(nargin,1,mfilename);

% set defaults
pct = [0,100];
globalFlag = false;
logFlag = false;
normfun = @nanmean;
optArgs = {};

% deal with the various inputs
if nargin > 2
    if rem(nargin,2) == 0
         error(message('bioinfo:manorm:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'percentile','global','logdata','method',...
        'extra_args','structoutput','newcolumnname','prctile'};
    for j=1:2:nargin-2
        [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
        switch(k)
            case {1,8}  % percentile
                if ~isnumeric(pval)
                    error(message('bioinfo:manorm:PercentileMustBeNumeric'));
                end
                if numel(pval) == 1
                    pct(2) = pval;
                elseif numel(pval) == 2  %range
                    pct = pval;
                    if pct(2)<=pct(1)
                        error(message('bioinfo:manorm:PercentileInvalidRange'));
                    end
                else
                    error(message('bioinfo:manorm:PercentileBadSize'));
                end
                if any(pct < 0) || any(pct > 100)
                    error(message('bioinfo:manorm:PercentileMustBe0To100'));
                end
            case 2 % global flag
                globalFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            case 3 % log data
                logFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            case 4 % method
                normmethod = pval;
                if ischar(normmethod)
                    okmethods = {'mean','median','std','mad'};
                    nm = strmatch(lower(normmethod), okmethods);
                    if isempty(nm)
                        error(message('bioinfo:manorm:UnknownMethodName', pval));
                    elseif length(nm)>1
                        error(message('bioinfo:manorm:AmbiguousMethodName', pval));
                    else
                        switch nm
                            case 1
                                normfun = @nanmean;
                            case 2
                                normfun = @nanmedian;
                            case 3
                                normfun = @nanstd;
                            case 4
                                normfun = @mad;
                                optArgs = {1};
                        end
                    end
                    
                elseif isa(normmethod, 'function_handle')
                    normfun = normmethod;
                else
                    error(message('bioinfo:manorm:MethodNotFunctionHandle'));
                end
            case 5 % optional arguments
                optArgs = pval;
                
                if ~iscell(optArgs)
                    optArgs = {pval};
                    
                end
            case {6,7}  % struct options
                warning(message('bioinfo:manorm:StructOnlyOptions', pname));
        end
    end
end
dmFlag = isa(X, 'bioma.data.DataMatrix');
    
% for global values, simply make the data into a single column and reshape
% later.
if globalFlag
    origSize = size(X);
    if dmFlag
        xRowNames = rownames(X);
        xColNames = colnames(X);
        xName = X.Name;
        X = X.(':')(':');
    end
    
    X = X(:);
end

% deal with percentile data
Xtemp = X;
if pct(1)>0 || pct(2) <100
    lowMask = X < repmat(prctile(X,pct(1)),size(X,1),1);
    highMask = X > repmat(prctile(X,pct(2)),size(X,1),1);
    Xtemp(highMask | lowMask) = NaN;
end
% calculate the function
gnorm = feval(normfun,Xtemp,optArgs{:});

% do we subtract or divide?
if logFlag
    X = X - repmat(gnorm,size(X,1),1);
else
    X = X./repmat(gnorm,size(X,1),1);
end

% reshape global data
if globalFlag
    if dmFlag
        X = bioma.data.DataMatrix(X, xRowNames, xColNames, 'Name', xName);
    else
        X = reshape(X,origSize);
    end
end
