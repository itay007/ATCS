function maxyplot(Exprs, Labels, inspectIdx)
% MAXYPLOT Pairwise comparisons of the MA and XY plots for the arrays.
% 
% MAXYPLOT(EXPRS, LABELS, ARRAYIDX) plots MA and XY scatter plots for all
% pairwise comparisons of intensities of ARRAYIDX. EXPRES is the intensity
% matrix with each column representing an array. LABELS is a cell array of
% the array labels or the column labels for the expression matrix. ARRAYIDX
% is a vector of the indices of the arrays to be compared. If EXPRS is a
% DataMatrix, LABELS can be omitted. A MA scatter plots the log2 ratio (M)
% between intensities X and Y against their log2 mean intensity (A). The MA
% plots are shown at the top right corner. A XY scatter plots the log2
% intensities of X vs. log2 intensities of Y. The XY plots are shown at the
% lower left corner of the comparisons.  
%
%   Example:
%        inspectIdx = 1:3;
%        maxyplot(Exprs, sampleLabels, inspectIdx)

%   Copyright 2009 The MathWorks, Inc.

bioinfochecknargin(nargin, 2, mfilename)

if nargin == 2
    if  ~isa(Exprs, 'bioma.data.DataMatrix')
        error(message('bioinfo:maxyplot:NeedDataMatrixInput'));
    else
        inspectIdx = Labels;
        Labels = Exprs.colnames;
        Exprs = Exprs.(':')(':');
    end
end

if isa(Exprs, 'bioma.data.DataMatrix')
   if isempty(Labels)
       Labels = Exprs.colnames;
   end
   
   Exprs = Exprs.(':')(':');
end

if ~isnumeric(Exprs) || size(Exprs,2) < 2 || size(Exprs,1)<1000
    error(message('bioinfo:maxyplot:InvalidExpressionMatrix'));
end

if  ~iscellstr(Labels)  
    error(message('bioinfo:maxyplot:LabelsNotACellArray'));
end

if ~isnumeric(inspectIdx) || ~isvector(inspectIdx)
    error(message('bioinfo:maxyplot:InspectSampleIndicesNotANumericVector'));
end

[nRows, nCols] = size(Exprs);

if ~all(ismember(inspectIdx, 1:nCols))
    error(message('bioinfo:maxyplot:InspectSampleIndicesContainInvalidIndex'));
end

dispLen = nRows;
if nRows > 10000
    dispLen = nRows;
end

X = log2(Exprs(1:dispLen, inspectIdx));
sampleToInspect = Labels(inspectIdx);
n = length(inspectIdx);
c = 0;
figure;
for i = 1:n
    for j = 1:n
        c = c+1;
        subplot(n,n,c)
        if i == j 
            text(0.5, 0.5, sampleToInspect{i},...
                'HorizontalAlignment','center', 'Interpreter', 'none');
            set(gca, 'xtick', [], 'ytick', [], 'yticklabel',[], 'box', 'on');
        elseif j-i > 0
            ym = X(:,i) - X(:,j);
            xa = (X(:,i)+ X(:,j))/2;
            x_lim = [min(xa) max(xa)];
            y_lim = [min(ym) max(ym)];
            plot(xa, ym, '.', x_lim, [0 0], 'r--');
            xlim(x_lim)
            ylim(y_lim)
            if j == n
                set(gca, 'yaxislocation', 'right');
            else
                set(gca, 'ytick', [], 'yticklabel',[]);
            end
            
            if i == 1
                set(gca, 'xaxislocation', 'top')
            else
                set(gca, 'xtick', [], 'xticklabel',[]);
            end
        else
            xy_max = max(max(X(:, i)), max(X(:, j)));
            xy_min = min(min(X(:, i)), min(X(:, j)));
            dnl = linspace(xy_min, xy_max);
            plot(X(:, i), X(:, j), '.', dnl, dnl, 'r--');
            axis([xy_min, xy_max, xy_min, xy_max ]);
            
            if i ~= n
                set(gca, 'xtick', [], 'xticklabel',[]);
            end          
            if j ~= 1
                set(gca, 'ytick', [], 'yticklabel',[]);
            end
        end
    end
end
end
