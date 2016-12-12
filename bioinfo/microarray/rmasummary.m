function expressionData = rmasummary(indices, values, varargin)
%RMASUMMARY calculates the gene expression values by using Robust
% Multiple-array Average (RMA) procedure.
% 
%   For a given probe set n with J probe pairs, let Yijn denote the
%   background adjusted, base-2 log transformed and quantile-normalized PM
%   value of chip i and probe j. Yijn follows a linear additive model 
%
%       Yijn = Uin + Ajn + Eijn, i = 1,...,I, j = 1,...,J, n=1,...,N;
%   where 
%       Uin - gene expression of the probe set n on chip i, 
%       Ajn  - probe affinity effect for the jth probe in the probe set, 
%       Eijn - residual for the jth probe on the ith chip.
%   The RMA methods assumes A1 + A2 +...+ AJ = 0 for all probe sets. A
%   robust procedure, median polish, is used to estimate Ui as the log
%   scale measure of expression. 
%   
%   EXPRESSIONDATA = RMASUMMARY(INDICES, DATA), where INDICES correspond to
%   the probe indexing, and DATA is the probe level intensity matrix in
%   natural scale and the columns correspond to separate chips. The rows of
%   DATA correspond to probes indices. The returned EXPRESSIONDATA is a
%   base-2 log scale expression measure for each probe set. The rows of
%   EXPRESSIONDATA correspond to probe sets, and the columns correspond to
%   separate chips.
% 
%   EXPRESSIONDATA = RMASUMMARY(..., 'OUTPUT', TYPE), returns output
%   expression values in the form defined by TYPE. TYPE can be one of
%   'log', 'log2' (default), 'log10', 'linear' scales, or a function
%   handle. If a function handle is passed, the output is transformed by
%   the transformation defined by the function.
% 
%   Examples:
%       load prostatecancerrawdata
%       backadjData = rmabackadj(pmMatrix);
%       normData = quantilenorm(backadjData);
%       expressionData = rmasummary(probeIndices, normData);
%
%   See also  AFFYGCRMA, AFFYINVARSETNORM, AFFYPREPROCESSDEMO, AFFYRMA,
%   CELINTENSITYREAD, GCRMA, GCRMABACKADJ, MAFDR, MANORM, MATTEST,
%   QUANTILENORM, RMABACKADJ.

% Copyright 2003-2008 The MathWorks, Inc.


% References: 
% [1] Irizarry RA, Hobbs B, Collin F, Beazer-Barclay YD, Antonellis KJ,
%     Scherf U, Speed TP. "Exploration, Normalization, and Summaries of
%     High Density Oligonucleotide Array Probe Level Data" Biostatistics
%     2003, 4, pp249-264.
% [2] Mosteller F, Tukey J. Data Analysis and Regression, Addison-Wesley
%     Publishing Company, Reading, Mass., 1977, p.165-202.


% Validate Indices, values
bioinfochecknargin(nargin,2,mfilename);

if ~isnumeric(values) || ~isreal(values)
   error(message('bioinfo:rmasummary:ProbeIntensityNotNumericAndReal')) 
end

if ~isnumeric(indices) || ~isreal(indices) || ~isvector(indices)...
        || any(indices < 0) || any(isnan(indices))
   error(message('bioinfo:rmasummary:InvalidProbeIndices')) 
end

if isvector(values)
  values = values(:);
end

% Columnized indices
indices = indices(:);

if numel(indices) ~= size(values,1)
   error(message('bioinfo:rmasummary:NotEqualNumberOfProbes'))
end

nProbes = numel(indices); % Total number of probes
psIndices = find(indices == 0); %Indices of probe sets. Probe pair always start from 0 
nProbeSets = numel(psIndices); %number of probe set 
nProbePairs = 0; %#ok Number of probe pair in a probe set %#ok
nChips = size(values,2); % Number of Chips, the first column is the probe pair number
valClass = class(values);
sumData = zeros(nProbeSets, nChips, valClass);
output = 2; % Output in log2 transformed scale
outputFun = @log2;

% deal with the various inputs

if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:rmasummary:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'output'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:rmasummary:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:rmasummary:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % output type
                    outputtype = pval;
                    if ischar(outputtype)
                        okmethods = {'log','log2','log10','natural','linear'};
                        output = strmatch(lower(outputtype), okmethods, 'exact');
                        if isempty(output)
                            error(message('bioinfo:rmasummary:UnknownOutputType', pval));
                        elseif length(output)>1
                            error(message('bioinfo:rmasummary:AmbiguousOutputType', pval));
                        else
                            switch output
                                case 1
                                    outputFun = @log;
                                case 2
                                    outputFun = @log2;
                                case 3
                                    outputFun = @log10;
                                case 4 %'natural' is deprecated (R2012b)
                                    error(message('bioinfo:rmasummary:incompatibleOutputScale'));
                                case 5
                                    outputFun = @(x)x; %No op
                            end
                        end
                    elseif isa(outputtype, 'function_handle')
                        outputFun = outputtype;
                        output = 0;
                    else
                        error(message('bioinfo:rmasummary:OutputTypeNotFunctionHandle'));
                    end
            end
        end
    end
end

for i = 1 : nProbeSets
    if i < nProbeSets
        nProbePairs = psIndices(i+1) - psIndices(i);
    else
        nProbePairs = nProbes - psIndices(i) + 1;
    end
    
    ppstart = psIndices(i); % start probe pair index in a probe set 
    ppend = ppstart+nProbePairs-1; % end prope pair index in a probe set
    
    % Get a matrix of probe values for a probe set for all the chips. Here
    % the row is the rpobe values and the column is the chips
    ppvalues = values(ppstart:ppend, 1:nChips); 
    
    sumData(i, :) = doMedianPolish(ppvalues, valClass);
end

if output == 2
    expressionData = sumData;
else
    try
        sumData = feval(@pow2, sumData); %Convert from 'log2' back to 'natural'
        expressionData = feval(outputFun, sumData);
    catch theException
        % Output in default log2 scale, and warn user about the error
        expressionData = feval(@log2, sumData);
        warning(message('bioinfo:rmasummary:OutputTypeFunctionErrored', theException.message));
    end
end
end % rmasummary


%---------------------- Median polish -------------------------------------
%  Median Polish. Median is more resistant to outliers in the data than mean, it
%  is robust to outliers. Each value in an additive model is decomposed to
%  several parts:
% 
%       data = overllEffects + rowEffects + columnEffects + residual
% 
%  The median polish algorithm consisted these steps: 
%   1. Subtract the medians of each row from the row values, record the row
%      medians in the row effect. 
%   2. Subtract the medians of each column from the column values, and record
%      the column medians in the column effect. 
%   3. Repeat these steps until it converges, i.e. the medians of both rows and
%      columns are 0.0 (+/- 0.5). The function allow for up to 10 iterations. 
%   4. Compute the median of row median and add it to the overall effect,
%      subtract it from the row median. Then do the same for the column effect.
% 
%  The function returns the overall and column effects as the expression measure
%  for a given probe set.

function avgRet = doMedianPolish(ppvalues, valclass)
[rows, cols] = size(ppvalues);

% log2 transform
z = log2(ppvalues);

% Define maximum number of iteration if the median polish procedure does not
% converge
maxIteration = 10; % This is the number RMAExpress uses
tol = 0.01;
allEffect = 0.0;
rowEffects = zeros(rows,1, valclass); % For row (chip) effects 
colEffects = zeros(1,cols, valclass); % For col (probe) effects

iloop = 1;
row_median = median(z, 2);
while iloop <= maxIteration
    % subtract the matrix by the row median
    z = bsxfun(@minus, z, row_median);
    % Add to row effect
    rowEffects = rowEffects + row_median;
    % Get column median
    col_median = median(z, 1);
    
    % Subtract the matrix by the column median
    z = bsxfun(@minus, z, col_median);
    %Add to column effect
    colEffects = colEffects + col_median;
    % Get row median
    row_median = median(z, 2);
    
    % check for convergence
    if all(abs(row_median) < tol ) && all(abs(median(z, 1)) < tol )
        break;
    end
     
    iloop = iloop+1;
end

% The commented out code below is the formal way of summing up the effects,
% since the row effect will not be added to the final results, so we can cut a
% few steps.
% Add median of the row effect to the all effect
me = median(rowEffects);

allEffect = allEffect + me;
% % Subtract it from row effect
% rowEffects = rowEffects - me;
% 
% % Add median of the col effect to the all effect
% me = median(colEffects);
% allEffect = allEffect + me;
% % Subtract it from col effect
% colEffects = colEffects - me;

% The total probe (row effect) should be zero
avgRet = colEffects + allEffect;
end % doMedianPolish







    


