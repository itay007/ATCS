function normY = mainvarsetnorm(X, Y, varargin)
%MAINVARSETNORM performs rank invariant set normalization. 
%
%   NORMY = MAINVARSETNORM(X,Y), where X and Y correspond to expression
%   values. X and Y values are converted to ranks separately. The invariant
%   ranks are selected by proportional rank difference below a given
%   threshold. Y is normalized based on the rank invariant set. 
%   
%   Note: if X or Y contain NaN values, then YNORM will also contain NaNs
%   at the corresponding positions.
% 
%   Notes: 1) The rank invariant set is a set of data points whose
%   proportional rank difference is smaller than a given threshold. The
%   threshold for each data point is determined by interpolating between a
%   threshold for the lowest average ranks and a threshold for the highest
%   average ranks. The average rank for each data point is determined by
%   averaging the two ranks in X and Y for each data point. Select these
%   two thresholds empirically to limit the spread of the invariant set,
%   but allow enough data points to determine the normalization
%   relationship. 2) The selection procedure is iteratively applied to the
%   newly selected set until the size of the new set does not decrease or
%   reaches the predetermined stop limit.
%
%   MAINVARSETNORM(...,'THRESHOLDS',THRD) sets a 1-by-2 vector of
%   thresholds for the proportional rank difference of the lowest and the
%   highest averaged ranks during invariant set selection. Default is
%   [0.03, 0.07]. The values must be between 0 to 1.
%   
%   MAINVARSETNORM(...,'EXCLUDE',N) uses an invariant set excludes the N
%   highest and lowest average ranks between X and Y. Default is 0.
% 
%   MAINVARSETNORM(...,'PERCENTILE',PCT) stops the iteration when the
%   number of points in the invariant set is PCT percentile of the total
%   number of input points. Default is 1.
% 
%   MAINVARSETNORM(...,'ITERATE',TF)iterates until the size of the
%   invariant set does not decrease or reaches the stop limit (PCT) when TF
%   is set to true. Set this option to false to compute the invariant set
%   without iteration. Default is true.
%
%   MAINVARSETNORM(...,'METHOD',METHOD) sets a smoothing method used to
%   normalize the data. METHOD can be 'Lowess' (default) or 'runmedian' for
%   running median.
% 
%   MAINVARSETNORM(...,'SPAN',SPAN) allows you to modify the window size
%   for the smoothing method. If SPAN is less than 1, then the window size
%   is taken to be a fraction of the number of points in the data. If span
%   is greater than or equal to 1, then the window is of size SPAN. Default
%   is 0.05, which corresponds to a window size equal to 5% of the number
%   of points in the invariant set.
%
%   MAINVARSETNORM(...,'SHOWPLOT',TF) displays the M-A scatter plots before
%   and after the normalization when TF is set to true. Default is false.
%
%   Examples:
%       maStruct = gprread('mouse_a1wt.gpr');
%       cy5data = magetfield(maStruct,'F635 Median');
%       cy3data = magetfield(maStruct,'F532 Median');
%       cy3norm = mainvarsetnorm(cy5data,cy3data, 'showplot', true);
%
%   See also AFFYINVARSETNORM, MAIRPLOT, MALOGLOG, MALOWESS, MANORM,
%   QUANTILENORM.

% Copyright 2003-2008 The MathWorks, Inc.


% References: 
% [1] George C. Tseng, Min-Kyu Oh, Lars Rohlin, James C. Liao, and Wing
%     Hung Wong. (2001) Issues in cDNA microarray analysis: quality
%     filtering, channel normalization, models of variations and assessment
%     of gene effects. Nucleic Acids Research. 29: 2549-2557. 
% [2] Hoffmann,R., Seidl,T., and Dugas,M. (2002) Profound effect of
%     normalization on detection of differentially expressed genes in
%     oligonucleotide microarray data analysis. Genome Biology
%     3(7):research0033.1-0033.11.


% Set defaults
%Threshold of proportion rank difference (PRD)
prdTD = [0.03, 0.07];
stopPCT = 1; % 1%
span = 0.05;
exL = 0;
iterFlag = true;
methodIsLowess = true;
dispFlag = false;
isYCol = false;

% Validate inputs
bioinfochecknargin(nargin,2,mfilename);

if ~isnumeric(X) || ~isreal(X) || ~isvector(X) ||...
        ~isnumeric(Y) || ~isreal(Y) || ~isvector(Y)
    error(message('bioinfo:mainvarsetnorm:ExpressionValuesNotNumericAndReal'))
end

% Check is Y is a column
if size(Y, 2) > 1
    isYCol = true;
end

% Columnized the vector
X = X(:);
Y = Y(:);

if numel(X)~= numel(Y)
    error(message('bioinfo:mainvarsetnorm:differentSize'))
end

N = numel(Y);

% Deal with the various input
if nargin > 2
   if rem(nargin, 2) == 1
       error(message('bioinfo:mainvarsetnorm:IncorrectNumberOfArguments', mfilename));
   end
    
   okargs = {'thresholds', 'prctile', 'exclude', 'iterate', 'method', 'span', 'showplot','percentile'};
   
   for j = 1:2:nargin-2
       pname = varargin{j};
       pval = varargin{j+1};
       
       k = find(strncmpi(pname,okargs,numel(pname)));
       
       if isempty(k)
           error(message('bioinfo:mainvarsetnorm:UnknownParameterName', pname));
       elseif length(k) > 1
            error(message('bioinfo:mainvarsetnorm:AmbiguousParameterName', pname));
       else
            switch(k)
                case 1 % prd threshold
                    if numel(pval) ~= 2 
                        error(message('bioinfo:mainvarsetnorm:InvalidSizeVectorThreshold'));
                    elseif any(pval < 0 ) || any(pval > 1)
                        error(message('bioinfo:mainvarsetnorm:badVectorThreshold'));
                    else
                        prdTD = [pval(1), pval(2)];
                    end
                case {2,8} % stop Percentile
                    if ~isnumeric(pval)
                        error(message('bioinfo:mainvarsetnorm:percentileMustBeNumeric'));
                    elseif pval < 0 || pval > 100
                        error(message('bioinfo:mainvarsetnorm:percentileMustBe0To100'));
                    end
                    
                    stopPCT = pval;
                case 3 % exclude L points of highest ranks and L points of lowest ranks
                    if ~isnumeric(pval) || pval < 0
                        error(message('bioinfo:mainvarsetnorm:toBeExcludeNumberMustBePositive'));
                    end
                    
                    exL = fix(pval);
                    if exL < 0 || exL >= floor(N/2)
                       error(message('bioinfo:mainvarsetnorm:toBeExcludeNumberLessThanX')); 
                    end
                case 4 % iterate
                    iterFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 5 % smooth method
                    if ischar(pval)
                        okmethods = {'lowess', 'runmedian'};
                        nm = strmatch(lower(pval), okmethods);
                        if isempty(nm)
                            error(message('bioinfo:mainvarsetnorm:UnknownMethodName'));
                        elseif length(nm) > 1
                            error(message('bioinfo:mainvarsetnorm:AmbiguousMethodName', pval));
                        else
                            methodIsLowess = (nm==1);
                        end
                        
                    else
                        error(message('bioinfo:mainvarsetnorm:MethodNameNotValid'));
                    end
                case 6 % span
                    if ~isnumeric(pval)
                        error(message('bioinfo:mainvarsetnorm:SmoothMethodSpanNumeric'));
                    end
                    
                    span = pval;
                    if span < 0 
                        error(message('bioinfo:mainvarsetnorm:SmoothMethodSpanPositive'));
                    end
                    if span >numel(X)
                        error(message('bioinfo:mainvarsetnorm:SmoothMethodSpanMoreThanX'));
                    end
                case 7 % display
                    dispFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            end % switch   
       end % if
       
   end % for loop
end % if

% Find nans
normY = Y;
nanids = isnan(X) | isnan(Y);
normY(nanids) = NaN;

gX = X(~nanids);
gY = Y(~nanids);

% Get rank invariant set
iset = invariantrankselect(gX, gY, prdTD, stopPCT, exL, iterFlag);
iX = gX(iset);
iY = gY(iset);

% get smooth iY
if methodIsLowess
    iYS = bioinfoprivate.masmooth(iX, iY, span);
else
    iYS = marunmed(iX, iY, span);
end

% Get normalized normY
[iYS_u, idx_u] = unique(iYS);

if numel(iYS_u) <= 1 % IYS is all the same, don't interplate
    gnormY = gY;
else
    iX_u = iX(idx_u);
    gnormY = interp1(iYS_u, iX_u, gY, 'spline', 'extrap');
    gnormY = max(gnormY,0);
end 

normY(~nanids) = gnormY;
    
% display the estimated distributions of the columns
if dispFlag
    hfig = doPlot(gX, gY, iX, iY, gnormY, iYS, gnormY(iset));
    set(hfig, 'visible', 'on');
end

% Convert back to rows if input Y is row
if isYCol
    normY = normY';
end

%***********************************************
function hfig = doPlot(X, Y, iX, iY, normY, smooth_iY, norm_iY)
units = 'normalized';
hfig = figure('Units',units,'Visible', 'off');

axes('parent', hfig, 'Position',[0 0 1 1],'Visible','off');
text(0.45, 0.95, 'M-A plots');
text(0.18, 0.85, 'Before normalization');
text(0.68, 0.85, 'After normalization');

% do MA plot
[iX_sort, idx] = sort(iX);
[Mp, Ap] = calculateMA(X, Y);
[iMp, iAp] = calculateMA(iX, iY);
[sMp, sAp] = calculateMA(iX_sort, smooth_iY(idx));
[Ma, Aa] = calculateMA(X, normY);
[iMa, iAa] = calculateMA(iX, norm_iY);

x_lim = [min(Ap) max(Ap)];
y_lim = [min(Mp), max(Mp)];
hax_1 =axes('position', [0.08 0.1100 0.4 0.8150]);
plot(Ap, Mp, '.','markersize', 2);

hold on
% plot invariant set
hi1 = plot(iAp, iMp, 'or', 'markersize', 3);
hs = plot(sAp, sMp, 'g-');
plot(x_lim, [0 0], 'k--');
hl1 = legend([hi1, hs], 'Invariant set', 'Smooth curve');
set(hl1, 'Box', 'off', 'location', 'SouthEast');
setAxesProps(hax_1, x_lim, y_lim)
hold off

hax_2 =axes('position', [0.5703  0.1100  0.4  0.8150]);
plot(Aa, Ma, '.', 'markersize', 2);
hold on
% plot invariant set
hi2 = plot(iAa, iMa, 'or', 'markersize', 3);
plot([min(Aa) max(Aa)], [0 0], 'k--');
hl2 = legend(hi2, 'Invariant set');
set(hl2, 'Box', 'off', 'location', 'NorthEast');
setAxesProps(hax_2, [min(Aa) max(Aa)], y_lim)
hold off

%***************************************
function [M, A] = calculateMA(X1, X2)
% Calculate M, A pair for MA plot. 
% M - the log2 fold change of X1 and X2
% A - the average log intensity of X1 and X1

negId = (X1 <=0 ) | ( X2<=0 );
Y1 = log2(X1(~negId));
Y2 = log2(X2(~negId));

M = Y1 - Y2;
A = (Y1 + Y2)/2;


function setAxesProps(ax, xlim, ylim)
set(ax, 'fontsize', 9 ,...
    'PlotBoxAspectRatio', [1 1 1],...
    'Xlim',xlim, 'Ylim',ylim);
xlabel(ax, 'A');
ylabel(ax, 'M');
   
