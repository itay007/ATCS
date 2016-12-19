function z = malowess(x,y, varargin)
% MALOWESS smoothes microarray data using the Lowess (or Loess) method.
%
%   YSMOOTH = MALOWESS(X,Y) smoothes scatter data X,Y using the Lowess
%   smoothing method. The default window size is 5% of the length of X. X
%   and Y can be DataMatrix objects. If Y is a DataMatrix object, YSMOOTH
%   is a DataMatrix object with the same properties as Y.
%
%   MALOWESS(...,'ORDER',order) allows you to choose the order of the
%   algorithm. This can be 1 (linear fit) or 2 (quadratic fit). The default
%   order is 1. 
%   Note that the MATLAB Curve Fitting Toolbox refers to Lowess smoothing
%   of order 2 as Loess smoothing. 
%
%   MALOWESS(...,'ROBUST',TF) uses a robust fit when TF is set to true.
%   This option can take a long time to calculate.
%
%   MALOWESS(...,'SPAN',span) allows you to modify the window size for the
%   smoothing function. If span is less than 1, then the window size is
%   taken to be a fraction of the number of points in the data. If span is
%   greater than 1, then the window is of size span. The default value is
%   0.05, which corresponds to a window size equal to 5% of the number of 
%   points in X. 
%
%   Example:
%
%       maStruct = gprread('mouse_a1wt.gpr');
%       cy5data = magetfield(maStruct,'F635 Median');
%       cy3data = magetfield(maStruct,'F532 Median');
%       [x,y] = mairplot(cy5data, cy3data);
%       drawnow
%       ysmooth = malowess(x,y);
%       hold on;
%       plot(x,ysmooth,'rx');
%       ynorm = y - ysmooth;
%   
%   See also AFFYINVARSETNORM, MABOXPLOT, MAGETFIELD, MAIMAGE,
%   MAINVARSETNORM, MAIRPLOT, MALOGLOG, MANORM, QUANTILENORM, ROBUSTFIT.

% Copyright 2003-2008 The MathWorks, Inc.


% Reference: "Trimmed resistant weighted scatterplot smooth" by
% Matthew C Hutcheson.

%==Check inputs

bioinfochecknargin(nargin,2,mfilename);
span = .05;
method = 'lowess';
robustFlag = false;

% deal with the various inputs
if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:malowess:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'span','order','robust'};
    for j=1:2:nargin-2
        [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
        
        switch(k)
            case 1  % span
                if ~isnumeric(pval)
                    error(message('bioinfo:malowess:LowessSpanNumeric'));
                end
                span = pval;
                if span < 0
                    error(message('bioinfo:malowess:LowessSpanPositive'));
                end
                if span >numel(x)
                    error(message('bioinfo:malowess:LowessSpanLessThanX'));
                end
                
            case 2 %order
                if isnumeric(pval)
                    switch(pval)
                        case 1
                            method = 'lowess';
                        case 2
                            method = 'loess';
                        otherwise
                            warning(message('bioinfo:malowess:LowessBadNumericOrder'));
                    end
                    
                    
                elseif ischar(pval)
                    pval = lower(pval);
                    if ~isempty(strmatch(pval,'linear'))
                        method = 'lowess';
                    elseif ~isempty(strmatch(pval,'quadratic'))
                        method = 'loess';
                    else
                        warning(message('bioinfo:malowess:LowessBadCharOrder'));
                    end
                end
                
            case 3 % robust flag
                robustFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                
        end
    end
end

if robustFlag
    method = ['r',method];
end

if isa(x, 'bioma.data.DataMatrix')
    x = x.(':')(':');
end

dmFlag = isa(y, 'bioma.data.DataMatrix');
if dmFlag
    yRowNames = rownames(y);
    yColNames = colnames(y);
    yName = y.Name;
    y = y.(':')(':');
end

z = bioinfoprivate.masmooth(x,y,span,method);
z = reshape(z,size(x,1),size(x,2));
if dmFlag
   z = bioma.data.DataMatrix(z, yRowNames, yColNames, 'Name', yName); 
end
