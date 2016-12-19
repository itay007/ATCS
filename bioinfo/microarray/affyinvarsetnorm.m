function  [normalizedData, medStruct] = affyinvarsetnorm(values,varargin)
%AFFYINVARSETNORM performs rank invariant normalization of probe intensities
% from multiple Affymetrix DAT or CEL files.
%
%   NORMDATA = AFFYINVARSETNORM(DATA), where the columns of DATA correspond to
%   separate chips, normalizes the probe intensities of each column to a common
%   baseline column, which has the median of overall intensities. The
%   normalization procedure ranks the intensities for each column, and then
%   finds a rank invariant set between the to-be normalized column and the
%   baseline reference column. A lowess or running median curve is fitted
%   through the invariant set and used for normalization.
% 
%   Note: if DATA contains NaN values, then NORMDATA will also contain NaNs
%   at the corresponding positions.
%
%   Notes: 
%   1) The rank invariant set is a set of data points whose proportional rank
%   difference is smaller than a given threshold. The threshold for each data
%   point is determined by interpolating between a threshold for the lowest
%   average ranks and a threshold for the highest average ranks. Select these
%   two thresholds empirically to limit the spread of the invariant set, but
%   allow enough data points to determine the normalization relationship. 2) The
%   selection procedure is iteratively applied to the newly selected set until
%   the size of the new set does not decrease or reaches the predetermined stop
%   limit.
%   
%   AFFYINVARSETNORM(...,'BASELINE',INDEX) sets column INDEX in DATA as the
%   baseline column.
%
%   AFFYINVARSETNORM(...,'THRESHOLDS',THRD) sets a 1-by-2 vector of thresholds
%   for the proportional rank difference of the lowest and the highest averaged
%   ranks. Default is [0.05, 0.005]. The values must be between 0 to 1. 
%
%   AFFYINVARSETNORM(...,'STOPPERCENTILE',SPCT) stops the iteration when the number
%   of points in the invariant set is SPCT percentile of the total number of
%   input points. Default is 1.
%
%   AFFYINVARSETNORM(...,'RAYPERCENTILE',RPCT) sets RPCT percentile of the highest
%   ranked invariant set of data points to fit a straight line through, while
%   the remaining data points are fitted to a lowess or a running median smooth
%   curve.  Default is 1.5.
% 
%   AFFYINVARSETNORM(...,'METHOD',METHOD) sets a smoothing method used to
%   normalize the data. METHOD can be 'Lowess' (default) or 'runmedian' for
%   running median.
% 
%   AFFYINVARSETNORM(...,'SHOWPLOT',SP) displays the scatter plots and M-A plots
%   before and after normalization. The plots show the baseline data versus data
%   from column with index SP. SP may also be a vector contains the column
%   indices in DATA. Use 'ALL' for SP to show plots for all the columns
%   normalized.
%
%   [NORMDATA, MEDSTRUCT] = AFFYINVARSETNORM(...) returns a structure of the
%   baseline column index, and each column's intensity medians before and after
%   the normalization.
%
%   Examples:
%       load prostatecancerrawdata;
%       normData = affyinvarsetnorm(pmMatrix,'showplot',2);
%
%   See also CELINTENSITYREAD, MAINVARSETNORM, MALOWESS, MANORM, QUANTILENORM,
%   RMABACKADJ, RMASUMMARY.

% Copyright 2003-2006 The MathWorks, Inc.


% References: 
% [1] Cheng Li and Wing Hung Wong (2001b) Model-based analysis of
%     oligonucleotide arrays: model validation, design issues and standard error
%     application, Genome Biology 2(8): research0032.1-0032.11.
% [2] http://biosun1.harvard.edu/complab/dchip/normalizing%20arrays.htm#isn

bioinfochecknargin(nargin,1,mfilename);

% Set defaults
stopPCT = 1; % 1%
rayPCT = 1; % 1.5%
plotId = 0;
baseId = -1;
methodIsLowess = true;

% Validate input data
if ~isnumeric(values) || ~isreal(values)
    error(message('bioinfo:affyinvarsetnorm:InputValuesNotNumericAndReal'))
end

% Get number of probes per column, and number of columns
if isvector(values)
    values = values(:);
end

[np, nc] = size(values);

if nc <= 1 % do this because [] is not an vector
    normalizedData = values;
    warning(message('bioinfo:affyinvarsetnorm:SingleColumnOrEmptyInputData'));
    return;
end

colId = 1:nc;
%Threshold of proportional rank difference (PRD)
prdTD = estimatethreshold(np);

% Deal with the various input
if nargin > 1
    if rem(nargin, 2) == 0
        error(message('bioinfo:affyinvarsetnorm:IncorrectNumberOfArguments', mfilename));
    end
okargs = {'thresholds', 'stopprctile', 'rayprctile', 'baseline', 'method','showplot','stoppercentile','raypercentile'};
    for j = 1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};

        k = find(strncmpi(pname,okargs,numel(pname)));

        if isempty(k)
            error(message('bioinfo:affyinvarsetnorm:UnknownParameterName', pname));
        elseif length(k) > 1
            error(message('bioinfo:affyinvarsetnorm:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % prd threshold
                    if numel(pval) ~= 2
                        error(message('bioinfo:affyinvarsetnorm:sizeVectorThreshold'));
                    elseif any(pval < 0 ) || any(pval > 1)
                        error(message('bioinfo:affyinvarsetnorm:badVectorThreshold'));
                    else
                        prdTD = [pval(1), pval(2)];
                    end
                case {2,7} % stop Percentile
                    if ~isnumeric(pval)
                        error(message('bioinfo:affyinvarsetnorm:stopPercentileMustBeNumeric'));
                    elseif pval < 0 || pval > 100
                        error(message('bioinfo:affyinvarsetnorm:stopPercentileMustBe0To100'));
                    end

                    stopPCT = pval;
                case {3,8} % ray percentile
                    if ~isnumeric(pval)
                        error(message('bioinfo:affyinvarsetnorm:rayPercentileMustBeNumeric'));
                    elseif pval < 0 || pval > 100
                        error(message('bioinfo:affyinvarsetnorm:rayPercentileMustBe0To100'));
                    end

                    rayPCT = pval;
                case 4 % baseline
                    if ~isnumeric(pval)
                        error(message('bioinfo:affyinvarsetnorm:BaselineIndexNotNumeric'));
                    elseif ~any(colId == pval)
                        error(message('bioinfo:affyinvarsetnorm:InvalidBaselineIdndex'));
                    end                 
                    baseId = pval;
                case 5 % smooth method
                    if ischar(pval)
                        okmethods = {'lowess', 'runmedian'};
                        nm = strmatch(lower(pval), okmethods);
                        if isempty(nm)
                            error(message('bioinfo:affyinvarsetnorm:UnknownMethodName'));
                        elseif length(nm) > 1
                            error(message('bioinfo:affyinvarsetnorm:AmbiguousMethodName', pval));
                        else
                            methodIsLowess = (nm==1);
                        end
                        
                    else
                        error(message('bioinfo:affyinvarsetnorm:MethodNameNotValid'));
                    end    
                case 6 % showplot
                    if isnumeric(pval)
                        if isscalar(pval)
                            pval = double(pval);
                        end
                        
                        tobeplot = ismember(pval, colId);
                        if ~all(tobeplot)
                            warning(message('bioinfo:affyinvarsetnorm:InvalidPlotIndices', num2str( plotId( ~tobeplot ) )));
                        end
                        plotId = pval;
                        
                    elseif ischar(pval) && strcmpi(pval, 'all')
                        plotId = colId;
                    else
                        warning(message('bioinfo:affyinvarsetnorm:InvalidSP'))
                    end
            end % switch
        end % if
    end % for loop
end % if

% Find array index of baseline array - the median of overall probe intensities
allMedians = nanmedian(values, 1);
if baseId == -1
    medianRank = tiedrank(allMedians);
    rankMedian = fix(median(medianRank));
    baseId = find( fix(medianRank) == rankMedian, 1);
end

medStruct.baselineID = baseId;
medStruct.baselineMedian = allMedians(baseId);
medStruct.beforeNormMedians = allMedians';
medStruct.afterNormMedians = allMedians';

normalizedData = zeros(np, nc);
% The baseline column is not changed.
normalizedData(:, baseId) = values(:, baseId);

% Loop over the columns and normalize them referenced on baseline column
colId = colId( colId ~= baseId);
baseline = values(:, baseId);

if any(plotId == baseId)
    disp(sprintf('Column #%d is baseline, no plot is created for it', baseId));
end

for iloop = colId;
    preX = values(:,iloop);
    
    % Find nans
    normalizedData(:,iloop) = preX;
    nanids = isnan(preX) | isnan(baseline);
    normalizedData(nanids,iloop) = NaN;
    
    gbaseline = baseline(~nanids);
    gpreX = preX(~nanids);
    
    % Get invariant-rank set indices
    iset = invariantrankselect(gbaseline, gpreX, prdTD, stopPCT);

    iX = gpreX(iset);
    iBaseY = gbaseline(iset);

    [normX, smooth_iX_s] = normalization(gbaseline, gpreX, iBaseY, iX, rayPCT, methodIsLowess);
    normalizedData(~nanids,iloop) = normX;
    medStruct.afterNormMedians(iloop) = median(normX);
    % The scatter plots of before and after the normalization 
    if any(plotId == iloop)
        hfig = doPlots(gpreX, gbaseline, iX, iBaseY, normX, smooth_iX_s, normX(iset), iloop, baseId);
        set(hfig, 'visible', 'on');
    end
end

%********** Helper functions *******************
function [normX, smooth_iX] = normalization(Y, X, iY, iX, rayPCT, methodIsLowess)
% Normalize X with Y using invariant set iY and iX.
% Y is the reference data, and iY is the invariant set of the reference (or
% baseline)

ws = warning;
warning('off','stats:robustfit:RankDeficient');
warning('off','stats:statrobustfit:IterationLimit'); 

% set default
span = 0.1;

% Number of points
np = numel(Y);
inp = numel(iY);

% Sort for smooth and interp
[iY_sort, iID_sort] = sort(iY);
iX_sort = iX(iID_sort);

[X_sort, ID_sort] = sort(X);

% Divide the rank invariant set into two parts: high rank sets for fit a
% straight line through due to fewer point in the high rank region. Usually
% about 1-5% of the points, and the lower rank sets for smooth curve fitting.
if inp > 1000
    % The number of points to be fitted by running median
    iN = inp - fix(rayPCT * inp / 100);

    % Sort and fitting the invariant sets   
    iID_rm = iID_sort(iN);
    iX_rm = iX(iID_rm);
    
    % Divide invariant set into two sets
    iY_1 = iY_sort(1:iN);
    iX_1 = iX_sort(1:iN);

    iY_2 = iY_sort(iN+1:inp);
    iX_2 = iX_sort(iN+1:inp);
    
    % Fit the first set with running median
    if methodIsLowess
        smooth_iX_1 = bioinfoprivate.masmooth(iY_1, iX_1, span);
    else
        smooth_iX_1 = marunmed(iY_1, iX_1, span);
    end
    % Fit the high rank set with robust fit
    try
        rf_prm = robustfit(iY_2, iX_2);
        smooth_iX_2 = rf_prm(1)+ rf_prm(2)*iY_2;
    catch
        smooth_iX_2 = marunmed(iY_2, iX_2, span);
    end
    smooth_iX = [smooth_iX_1; smooth_iX_2];
    
    % Fit all the data with smooth results from the invariant set
    % Sort and divide the whole data set into two parts
    N = find(X_sort == iX_rm, 1, 'last');
    X_1 = X(ID_sort(1:N));
    X_2 = X( ID_sort(N+1:np));

    % Fit the first set
    [iX_u, idx] = unique(smooth_iX_1);
    iBY_u = iY_1(idx);
    normX_1 = interp1(iX_u, iBY_u, X_1, 'spline');

    % Fit the high rank set
    [iX_u, idx] = unique(smooth_iX_2);
    if numel(iX_u) <= 1 % iX_u_2 are all the same, don't interplate
        normX_2 = X_2;
    else
        iBY_u = iY_2(idx);
        normX_2 = interp1(iX_u, iBY_u,  X_2, 'spline', 'extrap');
    end
    % Normalized results
    normX = [normX_1; normX_2];
    normX(ID_sort) = normX;
else % for niset <=1000, only do running median
    smooth_iX = marunmed(iY_sort, iX_sort, span);
    [iX_u, idx_u] = unique(smooth_iX);
    iBY_u = iY_sort(idx_u);
    normX = interp1(iX_u, iBY_u, X, 'spline');
end

warning(ws);

%***********************************************
function hfig = doPlots(X, Y, iX, iY, normX, smooth_iX, norm_iX, plotID, refID)
iY_sort = sort(iY);
xy_max = max(max(X), max(Y));
xy_min = min(min(X), min(Y));

dnl = linspace(xy_min, xy_max, 3);
units = 'normalized';
fw = 0.5;
fh = 0.5;
if strncmpi(computer, 'GLNX', 4)
    fw = 0.75;
    fh = 0.75;
end
hfig = figure('Units',units,...
    'Position',[.25 .10 fw fh],...
    'Visible', 'off');

axes('parent', hfig, 'Position',[0 0 1 1],'Visible','off');
labstr{1} = sprintf('X: #%d data', plotID);
labstr{2} = sprintf('Y: #%d baseline', refID);
mastr{1} = sprintf('X: A');
mastr{2} = sprintf('Y: M');
text(0.02, 0.7, labstr);
text(0.02, 0.33, 'M-A plots');
text(0.03, 0.25, mastr);
text(0.28, 0.957, 'Before normalization');
text(0.7, 0.957, 'After normalization');

hax_1 = axes('position', [0.18 0.53 0.38 0.38]); % Before normalization
plot(X, Y, '.',  dnl, dnl, 'k--', 'markersize', 2);
axis([xy_min, xy_max, xy_min, xy_max ]);
hold on
% plot invariant set
hi1 = plot(iX, iY, 'or', 'markersize', 3);

% plot smooth curve
hs1 = plot(smooth_iX, iY_sort, 'g');
hl1 = legend([hi1, hs1], 'Invariant set', 'Smooth curve');
set(hl1, 'Box', 'off', 'location', 'SouthEast');
setNatrualAxes(hax_1)
hold off

hax_2 = axes('position', [0.6 0.53 0.38 0.38]);
% plot scatter plot of input data and diagnol
plot(normX, Y, '.',  dnl, dnl, 'k--', 'markersize', 2);
hold on
% plot invariant set
hi2 = plot(norm_iX, iY, 'or', 'markersize', 3);
hl2 = legend(hi2, 'Invariant set');
set(hl2, 'Box', 'off', 'location', 'SouthEast');
axis([xy_min, xy_max, xy_min, xy_max ]);
setNatrualAxes(hax_2)
hold off

% do MA plot
[Mp, Ap] = calculateMA(Y,X);
[iMp, iAp] = calculateMA(iY, iX);
[sMp,sAp] = calculateMA(iY_sort, smooth_iX);
[Ma, Aa] = calculateMA(Y, normX);
[iMa, iAa] = calculateMA(iY, norm_iX);

x_lim = [min(Ap) max(Ap)];
y_lim = [min(Mp), max(Mp)];
hax_3 =axes('position', [0.18 0.06 0.38 0.38]);
plot(Ap, Mp, '.','markersize', 2);
hold on
% plot invariant set
plot(iAp, iMp, 'or', 'markersize', 3);
plot(sAp, sMp, '-g');
plot(x_lim, [0 0], 'k--');
set(hax_3, 'fontsize', 8, 'PlotBoxAspectRatio', [1 1 1],...
    'Xlim',x_lim, 'Ylim',y_lim)
hold off

hax_4 =axes('position', [0.6 0.06 0.38 0.38]);
plot(Aa, Ma, '.', 'markersize', 2);
hold on
% plot invariant set
plot(iAa, iMa, 'or', 'markersize', 3);
plot([min(Aa) max(Aa)], [0 0], 'k--');
set(hax_4,'fontsize', 8, 'PlotBoxAspectRatio', [1 1 1],...
    'Xlim',[min(Aa) max(Aa)], 'Ylim',y_lim)
hold off

%*************************************
function thds = estimatethreshold(npoints)
thdL = 0.05;
thdH = 0.005;
maxNH = ceil(npoints/10^5); 
thds = [thdL thdH/maxNH];

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

%*************************************
function fixticks(ax)
% This function will perform the desired ticklabel fixing
ytick=get(ax,'ytick');
% Convert ticks to a string of the desired format
ytickstr=num2str(ytick',7);
% Reset the labels to the new format
set(ax,'yticklabel',ytickstr);
xtick=get(ax,'xtick');
% Convert ticks to a string of the desired format
xtickstr=num2str(xtick',7);
% Reset the labels to the new format
set(ax,'xticklabel',xtickstr);

%************************************************
function setNatrualAxes(ax)
ticks = get(ax, 'Xtick');
if(numel(ticks) > 4)
    set(ax, 'Xtick', ticks(1:2:end));
end
fixticks(ax)
set(ax, 'fontsize', 9 ,'PlotBoxAspectRatio', [1 1 1]);











