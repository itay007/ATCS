function outCNV = cghcbs(cghdata, varargin)
% CGHCBS performs circular binary segmentation (CBS) on array CGH data.
%
%   S = CGHCBS(DATA) analyzes array CGH data using the circular binary
%   segmentation (CBS) algorithm. DATA can be a structure containing the
%   following fields:
%       Sample (optional)
%       Chromosome
%       GenomicPosition
%       Log2Ratio
%   or a matrix of data with the first column corresponding to chromosome
%   number, the second column corresponding to the genomic position, and
%   the third and higher columns corresponding to log2 ratio of test to
%   reference intensities. If the 'Sample' field is missing, a default
%   sample name will be assigned to output, S. S is the segment data
%   structure with the following fields:
%       Sample
%       SegmentData
%   The SegmentData field is a structure array containing segment data for
%   the sample, in the following fields:
%		Chromosome
%		Start
%		End
%		Mean
%
%   Note: The CBS algorithm recursively splits chromosomes into segments
%   based on a maximum t-statistic estimated by permutation. The
%   computation can be time consuming. If N is number of data points, the
%   computation time is ~ O(N^2).
%
%	CGHCBS(..., 'ALPHA', ALPHA) sets the significance level for the
%	statistical tests to accept change points. The default is 0.01.
%
%	CGHCBS(..., 'PERMUTATIONS', NPERM) sets the number of permutations used
%	for p-value estimation. The default NPERM is 10,000.
%
%	CGHCBS(...,'METHOD', METHOD) sets the method to estimate the p-values.
%	METHOD can be 'PERM' for full permutation or 'HYBRID' (default) for the
%	faster tail probability-based permutation. When using the 'HYBRID'
%	method, the 'PERM' method is applied automatically when segment data
%	length becomes less than 200.
%
%	CGHCBS(..., 'STOPPINGRULE', true) uses a heuristic stopping rule to
%	declare a change without performing the full number of permutations for
%	the p-value estimation. Default is false.
%
%	CGHCBS(..., 'SMOOTH', TF) smoothes outliers before segmenting. The
%	default is TRUE.
%
%	CGHCBS(..., 'PRUNE', TF) eliminates change points identified due to
%	local trends in the data that are not indicative of real copy number
%	change. The default is FALSE.
%
%	CGHCBS(...,'ERRSUM', GM) sets the allowed proportional increase in the
%	error sum of squares when eliminating change points. The value usually
%	is set to 0.05 or 0.1. Default GM = 0.05. This is ignored if the
%	'PRUNE' option is FALSE.
%
%	CGHCBS(..., 'WINDOWSIZE', WS) specifies the size of the window used to
%	divide data when using the 'PERM' method on large dataset. Default
%	WS = 200. This is ignored when using the 'HYBRID' method.
%
%	CGHCBS(..., 'SAMPLEINDEX', INDICES) analyzes only the sample(s)
%	specified by INDICES. It can be a single sample index or a vector of
%	sample indices. By default all samples in DATA will be analyzed.
%
%	CGHCBS(..., 'CHROMOSOME', CHR) analyzes only the data of chromosome(s)
%	specified by CHR. It can be a single chromosome number or a vector of
%	chromosome numbers. By default all the chromosomes in DATA will be
%	analyzed.
%
%	CGHCBS(..., 'SHOWPLOT', SP) plots segment means over the original data.
%	When SP is TRUE, all chromosomes and all samples in DATA are plotted.
%	If there are multiple samples in DATA, a separate figure will be
%	generated for each sample. If CGHCBS is called without output
%	arguments, plot(s) will be shown unless SP is FALSE. SP can specify the
%	layout of the plot. Use 'W' to plot all chromosomes in one plot or 'S'
%	to plot each chromosome in a subplot. If SP is set to TRUE, 'W' is the
%	plot layout. To plot one chromosome, SP can also be an index to one of
%	the chromosomes in DATA.
%
%   CGHCBS(...,'VERBOSE', TF) turns off the display of progress report of
%   the analysis if TF is set to FALSE. Default is TRUE.
%
%   Example:
%
%       % Analyze and plot all chromosomes from sample #3 in
%       % the Coriell cell line study
%       load coriell_baccgh
%       S = cghcbs(coriell_data,'sampleind',3,'showplot',true);
%
%       % Analyze and plot chromosome 9 from sample #32 in the pancreatic
%       % cancer study, and add chromosome 9 ideogram to the plot
%       load pancrea_oligocgh
%       PS = cghcbs(pancrea_data,'sampleind',32,'chromosome',9,'showplot',9);
%       chromosomeplot('hs_cytoBand.txt', 9, 'addtoplot', gca)
%
%       % Display CNV data for chromosomes 10 and 11 from sample #3 of the
%       % Coriell cell line study aligned to chromosomes in human ideogram
%       % The genomic position in this study is in kilo base pair unit.
%       cnvStruct = struct('Chromosome', [10 11],...
%                          'CNVType', [2 1],...
%                          'Start', [S.SegmentData(10).Start(2),...
%                                    S.SegmentData(11).Start(2)]*1000,...
%                          'End',   [S.SegmentData(10).End(2),...
%                                    S.SegmentData(11).End(2)]*1000)
%       chromosomeplot('hs_cytoBand.txt', 'cnv', cnvStruct, 'unit', 2)
%
%   See also AFFYSNPCNVDEMO, BACACGHDEMO, CHROMOSOMEPLOT, CYTOBANDREAD.

%   Copyright 2007-2010 The MathWorks, Inc.


% References:
% [1] A.B. Olshen, E.S. Venkatraman, R. Lucito, and M. Wigler. "Circular
%     binary segmentation for the analysis of array-based DNA copy number
%     data", Biostatistics 5, 4, pp. 557-572, 2004.
% [2] E.S. Venkatraman, A.B. Olshen. "A Faster Circular Binary Segmentation
%     Algorithm for the Analysis of Array CGH Data", Bioinformatics, 2007.
% [3] E.S. Venkatraman, A.B. Olshen. "DNAcopy: A Package for Analyzing DNA
%     Copy Data",
%     www.bioconductor.org/repository/devel/vignette/DNAcopy.pdf, 2006.


bioinfochecknargin(nargin,1,mfilename);

sampleNames = [];
if isstruct(cghdata)
    if ~checkInputStruct(cghdata)
        error(message('bioinfo:cghcbs:BadFieldNames'));
    end
    
    numSample = size(cghdata.Log2Ratio, 2);
    
    % More than one sample, should have Sample field
    if isfield(cghdata,'Sample')
        sampleNames = cghdata.Sample;
        
        if numSample ~= numel(sampleNames)
            error(message('bioinfo:cghcbs:SampleDataSizeNotMatch'));
        end
    else
        warning(message('bioinfo:cghcbs:MissSampleFieldName'));
    end
    
    chromIDs = unique(cghdata.Chromosome);
    CID = cghdata.Chromosome(:);
    X = cghdata.GenomicPosition(:);
    if numSample == 1
        log2Ratio = cghdata.Log2Ratio(:);
    else
        log2Ratio = cghdata.Log2Ratio;
    end
    % matrix with one sample data
elseif isnumeric(cghdata) && isreal(cghdata) && size(cghdata, 2) >= 3
    numSample = size(cghdata, 2) - 2;
    chromIDs = unique(cghdata(:,1));
    CID = cghdata(:, 1);
    X = cghdata(:, 2);
    log2Ratio = cghdata(:, 3:end);
else
    error(message('bioinfo:cghcbs:WrongInpuDataFormat'));
end

if isempty(sampleNames)
    sid = num2str((1:numSample)');
    sampleNames = strcat({'Sample'},cellstr(sid));
end

% Parse the optional inputs
[alpha, nPerm,hybridFlag, smoothFlag, winSize, pruneFlag,gamma,...
    segSampleId, segChromId, plotId, plotType, verbose,useStoppingRule] = ...
    parse_inputs(numSample, chromIDs, varargin{:});

% Hybrid p-value
if hybridFlag
    kMax = 25; % the maximum width of smaller segment for tail probability
    minHybridN = 200;
end

% Smooth outlier
if smoothFlag
    regionR = 2; % typical 2 to 5
end

% Return if no sample or chromsome is selected
if isempty(segChromId) || isempty(segSampleId)
    outCNV = [];
    return;
end

if (plotId > 0 ) && ~any(chromIDs==plotId)
    warning(message('bioinfo:cghcbs:InvalidPlotIndex'))
    plotId = -1;
end

if ~ismember(plotId, segChromId) && plotId > 0
    segChromId = [segChromId; plotId];
    
    warning(message('bioinfo:cghcbs:PlotIDNotAmongChromosomeIDs'))
end

numChroms = numel(segChromId);
numSamples = numel(segSampleId);
sampleNames = sampleNames(segSampleId(:));
log2Ratio = log2Ratio(:, segSampleId(:));


if nargout == 0
    if plotId > 0
        numChroms = 1;
        segChromId = chromIDs(chromIDs==plotId);
    elseif plotId < 0
        return;
    end
end
%-----------End of input checking---------------------------%
outCNV.Sample = [];
outCNV.SegmentData = [];
outCNV = repmat(outCNV,1,numSamples);

for smploop = 1:numSamples
    outCNV(smploop).Sample = sampleNames{smploop};
    outCNV(smploop).SegmentData.Chromosome = [];
    outCNV(smploop).SegmentData.Start = [];
    outCNV(smploop).SegmentData.End = [];
    outCNV(smploop).SegmentData.Mean = [];
    outCNV(smploop).SegmentData = repmat(outCNV(smploop).SegmentData,1,numChroms);
    
    Y = log2Ratio(:,  smploop);
    numSegs = zeros(numChroms, 1);
    for chrloop = 1:numChroms
        if verbose
            fprintf('Analyzing: %s. Current chromosome %d\n', ...
                sampleNames{smploop}, segChromId(chrloop));
        end
        
        cidx = CID == segChromId(chrloop);
        x = X(cidx);
        y = Y(cidx);
        
        % Remove NaN points
        nadidx = isnan(y);
        if any(nadidx)
            y = y(~nadidx);
            x = x(~nadidx);
        end
        
        if ~isempty(y)
            % Smooth outliers
            if smoothFlag
                y = smoothOutliers(y, regionR);
            end
            N = numel(y);
            
            % initialize segment indices
            segidx = [1, N];
            k = length(segidx);
            
            % Recursively find change point in segments
            segcount = 0;
            cp_ends = [];
            while k > 1
                % Current segment from segidx(k-1) - segidx(k)
                % Segment length
                sn = segidx(k) - segidx(k-1) + 1;
                
                if sn >= 4
                    yseg = y(segidx(k-1):segidx(k));
                    
                    if hybridFlag && sn > minHybridN
                        [cptnum, cptloc] = hybridTernaryChangePoints(yseg, sn, nPerm, alpha, kMax, useStoppingRule);
                    else
                        wlen = min(sn, winSize);
                        [cptnum, cptloc] = permuteTernaryChangePoints(yseg, sn, nPerm, alpha, wlen, useStoppingRule);
                    end
                else
                    cptnum = 0;
                end
                
                if cptnum == 0
                    segcount = segcount + 1;
                    cp_ends(segcount) = segidx(k);%#ok
                    
                    segidx = segidx(1:k-1);
                else
                    segidx = [segidx(1:k-1), segidx(k-1)+cptloc-1, segidx(k)];
                end
                
                k = length(segidx);
            end
            
            cp_segidx = cp_ends(end:-1:1);
            if pruneFlag && segcount > 1
                [cp_segidx, segcount] = prune(y, cp_segidx, gamma);
            end
            
            % compute segment means
            seglen = diff([0, cp_segidx]);
            
            sstart = [1, cp_segidx(1:end-1)];
            cp_starts = x(sstart);
            cp_ends = x(cp_segidx);
            cp_means = zeros(segcount,1);
            for s= 1:segcount
                cp_means(s) = sum(y(sstart(s):cp_segidx(s)))/seglen(s);
            end
        else % If all data points are NaNs
            cp_starts = 0;
            cp_ends = 0;
            cp_means = NaN;
            segcount = 0;
        end
        
        outCNV(smploop).SegmentData(chrloop).Chromosome = segChromId(chrloop);
        outCNV(smploop).SegmentData(chrloop).Start = cp_starts;
        outCNV(smploop).SegmentData(chrloop).End = cp_ends;
        outCNV(smploop).SegmentData(chrloop).Mean = cp_means;
        numSegs(chrloop) = segcount;
        
        if segChromId(chrloop)== plotId
            plotSingleChromosomeSegments(x, y, cp_starts, cp_ends,...
                cp_means, segChromId(chrloop),  sampleNames{smploop})
        end
    end
    
    if plotId == 0
        if plotType == 1
            plotAllChromosomes(X,Y, CID, chromIDs,segChromId, numSegs,...
                outCNV(smploop).SegmentData, sampleNames{smploop});
        elseif plotType == 2
            plotAllChromosomeSub(X,Y, CID, chromIDs,segChromId,...
                outCNV(smploop).SegmentData, sampleNames{smploop})
        end
    end
    
end % end of sample loop

end % End of main function
%------------------------------------------------------------------------%
function [cptnum, cptloc] = hybridTernaryChangePoints(x, n, nperm, alpha, kmax, useStoppingRule)
% Find the number of change points and locations by hybrid p-values
% x - data
% n - number of data points
% nperm - number of permutation
% alpha - pvalue threshold, if computed p > alpha, null hypothesis is true.
% kmax - smaller segment width

cptnum = 0;
cptloc = [0 0];
[winnum, winloc, winlens, winisodd] = getWindowParams(n, n);

% Get observed pstat and pseg
[pstat_o, pseg_o] = getMaxTStat(x,winnum, winloc, winlens, winisodd);

% Compute lower bound P(T2>b) from estimated tail probabilities
pval_2 = computeTailProb(pstat_o, n, kmax);

if pval_2 > alpha; % there is no change
    return;
end

% Compute upper bound P(T_1>b)+P(T_2>b) by permutation
nrej = 0;
maxnrej = fix((alpha-pval_2)*nperm);
bailOutStep = max(nperm/10,maxnrej*3);
bailOut = bailOutStep;
pstat = zeros(1,nperm);
for ploop = 1:nperm
    permx = x(randperm(n));
    pstat(ploop) = getHybridMaxT(permx, n, kmax);
    
    if pstat(ploop) >= pstat_o
        nrej = nrej+1;
    end
    %
    % This code is a heuristic stopping rule for the permutation test.
    %
    if useStoppingRule && ploop == bailOut
        geVparams = gevfit(pstat(1:bailOut),sqrt(eps));
        prob = gevcdf(pstat_o,geVparams(1) , geVparams(2), geVparams(3));
        if prob == 1
            break;
        end
        if prob < 1-alpha
            useStoppingRule = false;
        end
        bailOut = bailOut + bailOutStep;
    end
    
    % Stop if we have reached nrej
    if nrej > maxnrej
        return;
    end
    % Stop if we can't reach maxnrej
    if (nperm-ploop) < (maxnrej -nrej)
        break
    end
end
% If not rejected, there are change points
[cptnum, cptloc] = getSegmentEndPoints(x, pseg_o, n, alpha, nperm);
end % end of function

%------------------------------------------------------------%
function [cptnum, cptloc] = permuteTernaryChangePoints(x, n, nperm, alpha, winlen, useStoppingRule)
% Find the number of change points and locations by Permutation method
% x - data
% n - number of data points
% nperm - number of permutation
% alpha - pvalue threshold, if computed p > alpha, null hypothesis is true.
% winlen - divided data by windows of length winlen
% winoverlap - percentage of window overlap.

cptnum = 0;
cptloc = [0 0];
[winnum, winloc, winlens, winisodd] = getWindowParams(n, winlen);

% Get observed pstat and pseg
[pstat_o, pseg_o] = getMaxTStat(x, winnum, winloc, winlens, winisodd);

nrej = 0;
maxnrej = fix(alpha*nperm);
bailOutStep = max(nperm/10,maxnrej*3);
bailOut = bailOutStep;
pstat = zeros(1,nperm);
for ploop = 1:nperm
    permx = x(randperm(n));
    pstat(ploop) = getPermTStat(permx, winnum, winloc, winlens, winisodd);
    
    if pstat(ploop) >= pstat_o
        nrej = nrej + 1;
    end
    %
    % This code is a heuristic stopping rule for the permutation test.
    %
    if useStoppingRule && ploop == bailOut
        geVparams = gevfit(pstat(1:bailOut),sqrt(eps));
        prob = gevcdf(pstat_o,geVparams(1) , geVparams(2), geVparams(3));
        if prob == 1
            break;
        end
        if prob < 1-alpha
            useStoppingRule = false;
        end
        bailOut = bailOut + bailOutStep;
    end
    if nrej > maxnrej
        return;
    end
    % Stop if we can't reach maxnrej
    if (nperm-ploop) < (maxnrej -nrej)
        break
    end
end

% If not rejected, there are change points
[cptnum, cptloc] = getSegmentEndPoints(x, pseg_o, winlen, alpha, nperm);
end % end of function

%-----------------------------------------------------%
function  [cptnum, cptloc] = getSegmentEndPoints(x, seg, winlen, alpha, nperm)
% Return the change points and location index of a segment
% seg - vector of [i, j, wloc]
% winlen - window length
% nperm - number of permutations
% If not rejected, there are change points

cptnum = 0;
cptloc = [0 0];

if seg(2) == winlen
    cptnum = 1;
    cptloc = seg(1) + seg(3);
else
    if seg(1) == 1
        cptnum = 1;
        cptloc = seg(2) + seg(3);
    else
        % Do permutation test on point i
        sn1 = seg(1); % segment point 1
        sn2 = seg(2); % segment point 2
        wl = seg(3); % Window location
        tpval = computePermPvalue(x(wl:wl+sn2-1), sn1, sn2-sn1, nperm);
        
        if tpval <= alpha
            cptnum = 1;
            cptloc(1) = seg(1) + seg(3);
        end
        
        % Do permutation tests on point j
        tpval = computePermPvalue(x(wl+sn1-1:wl+winlen-1), sn2-sn1, winlen-sn2+1, nperm);
        
        if tpval <= alpha
            cptnum = cptnum + 1;
            cptloc(2) = seg(2) + seg(3);
        end
    end
end
cptnum = cptnum - sum(cptloc >= winlen);
cptloc = cptloc(cptloc ~= 0);
end

%-----------------------------------------------------%
function [maxt, maxseg] = getMaxTStat(x,  wnum, wloc, wlens, wisodd)
%Return the max t-statistic and its segment location (i, j)
% x - data
% n - data length
% wlens - window length
% wnum - number of windows
% wloc - windows locations
% wisodd - window length is odd number
% maxseg = [i;j;wloc]

% Number of windows
wt = zeros(1,wnum);
wseg = zeros(3, wnum);

for w = 1:wnum
    [wt(w), cp_i, cp_l] =  getMaxT(x(wloc(w):wloc(w+1)), wlens(w), wisodd(w));
    
    if cp_i+cp_l <= wlens(w)
        wseg(1:2, w) = [cp_i; cp_i+cp_l];
    else
        wseg(1:2, w) = [((cp_i+cp_l)-wlens(w)); cp_i];
    end
    wseg(3, w) = wloc(w);
end

[maxt, idx] = max(wt);
maxseg = wseg(:, idx);
end % end of function
%---------------------------------------------------------%
function maxt = getPermTStat(x, wnum, wloc, wlens, wisodd)
% Number of windows
wt = zeros(1,wnum);
for w = 1:wnum
    wt(w) = getMaxT(x(wloc(w):wloc(w+1)), wlens(w), wisodd(w));
end
maxt = max(wt);
end % end of function

%------------------------------------------------------%
function [maxt, maxpt, maxarclen]=  getMaxT(x, n, isodd)
% Returns the maximal t-statistic T = max_1<=i<j<=m|T_ij|
% x - data
% n - data length
% isodd - data length is odd number
%
% maxt - maximal T_ij
% maxpt - i of max(T_ij)
% arclen - j-i arc length

sumx = sum(x);
sumsqx = sum(x.^2);
meanx = sumx/n;
sx = x-meanx;

% Mean squared error
mss = sumsqx - meanx * sumx;

% Number of arc_n, i counts from 2 up to n/2
arc_num = fix(n/2);
max_x = zeros(1, arc_num-1 - 1*isodd);

max_i = zeros(size(max_x));

% Let arc_n = j-i and loop through arc_num. For each fix arc_n, find
% the max Sj-Si by loop through i. The arc_n can be up to n/2.

ssx = [sx;sx];
cx0 = [0;cumsum(ssx)];
cx1 = cx0(1:n);
weight = (2*(n-2));
overallMax = n*(max(max(cx0)-min(cx0)))^2;
for arc_n = 2:arc_num
    cx = cx0(arc_n+1:n+arc_n);
    cx = cx - cx1;
    
    if ~isodd && arc_n==arc_num
        cx(arc_n+1:end) = cx(arc_n+1:end) - ssx(n:n+arc_n-1);
    end
    
    [max_x(arc_n-1), max_i(arc_n-1)] = max(abs(cx));
    max_x(arc_n-1) = n*max_x(arc_n-1)^2/weight;
    % Break out of the loop if we are better than the best possible next
    % value given the weighting scheme.
    weight = ((arc_n+1)*(n-(arc_n+1)));
    if max_x(arc_n-1) >= overallMax/weight
        break
    end
end

% Find maximal of max_x
[max_t, maxarclen] = max(max_x);

% Compute max(T_ij) by consider the mean square error
tss = mss-max_t;
tss(tss <= 0)= 1;
maxt = sqrt(max_t/(tss/(n-2)));

% Find the i and j for the max(T_ij)
maxpt = max_i(maxarclen);
end % end of function

%------------------------------------------------------%
function pval = computePermPvalue(x, n1, n2, nperm)
% Return p values from by permutate two groups n1 and n2
% x- data segment to be tested
% n1 - number of data points in segment1
% n2 - number of data points in segment2
% n = n1 + n2

n = n1 + n2;
% Observed segment mean diff
segmd_o = sum(x(1:n1))/n1 - sum(x(n1+1:n))/n2;

% Permutation segment mean diff
segmd_p = zeros(nperm, 1);
for i = 1:nperm
    rx = x(randperm(n));
    segmd_p(i) = abs(sum(rx(1:n1))/n1 - sum(rx(n1+1:n))/n2);
end
pval = sum(segmd_p >= abs(segmd_o))/nperm;
end % end of function

%------------------------------------------------------%
function maxt = getHybridMaxT(x, n, k)
% Returns the maximal t-statistic T1 = max_1<=i<j<=m|T_ij| of A_1. A_1 =
% {i,j: j-i <=k or > m-k}
%
% x - data
% n - data length
% k - the pre-set maximum width of smaller segment (arc) for permutation
% maxt - maximal T_ij in minor arc k

sumx = sum(x);
sumsqx = sum(x.^2);
meanx = sumx/n;

% Mean suqared error
mss = sumsqx - meanx * sumx;
sx = x-meanx;

% Number of l, i counts from 2 up to n/2
max_x = zeros(1, k-1);

ssx = [sx;sx];
cx0 = [0; cumsum(ssx)];
cx1 = cx0(1:n);
weight = (2*(n-2));
overallMax = n*(max(max(cx0)-min(cx0)))^2;
for arc_n = 2:k
    cx = cx0(arc_n+1:n+arc_n);
    cx = cx - cx1;
    
    max_x(arc_n-1) = max(abs(cx));
    max_x(arc_n-1) = n*max_x(arc_n-1)^2/weight;
    % Break out of the loop if we are better than the best possible next
    % value given the weighting scheme.
    weight = ((arc_n+1)*(n-(arc_n+1)));
    if max_x(arc_n-1) >= overallMax/weight
        break
    end
end

% Find maximal of max_x
max_t = max(max_x);

% Compute max(T_ij) by considering the mean square error
tss = mss-max_t;
tss(tss <= 0)= 1;
maxt = sqrt(max_t/(tss/(n-2)));
end % end of function

%----------------------------------------------%
function pval = computeTailProb(b, m, k)
% Compute tail probability for T_2>b
% b - oberseved p value
% m - total number of data points
% k - k max

persistent fhNU

if isempty(fhNU)
    fhNU = load('cghcbshybridnu.mat');
    fhNU = fhNU.NU;
end

pr = quadgk(@getNu, 0.5, 1-(k+1)/m);

pval = (b^3)*exp(-(b^2)/2) * pr/2;
    function sqnu = getNu(t)
        xv = b./sqrt(m.*t.*(1-t));
        nu = fhNU(xv(:));
        sqnu = (nu').^2;
    end

end % end of function

%--------------------------------------------------------%
function [segidx, segcount] = prune(x, segidx, cutoff)
% Prune the some of the non-real CN changes (local trend in data)
% x - log2ratio
% segidx - segment end indices found by cbs
% cutoff - proportional increase in sum of squares allowed when eliminating
%		   segments

seglen = diff([0, segidx]);
N = numel(seglen);

% Number of change points
C = N-1;

% Square sum of x
ssx = (sum(x.^2))';

% Find the sum of each segment
cx = (cumsum(x))';
sx = diff([0 cx(segidx)]);
% Get segment cumsums
csx = cumsum(sx);

% Sum of squared deviations of data points in segments around their segment
% average
SSC = ssx - sum((sx.^2) ./seglen);

cid_o = 1:C;
for c = C-1:-1:1
    cidx0 = allnr(C, c);
    n = size(cidx0, 1);
    cidx = [cidx0, N*ones(n, 1)];
    delta = diff([zeros(n,1), cidx], 1, 2);
    
    cdelta = cumsum(delta, 2);
    ccx = diff([zeros(n,1) csx(cdelta)], 1, 2);
    cclen = diff([zeros(n,1) segidx(cdelta)], 1, 2);
    % compute the squared deviations
    ssc_c = ssx - sum(ccx.^2 ./ cclen, 2);
    
    [SSC_c, idx] = min(ssc_c);
    cptloc = cidx(idx, :);
    
    if SSC_c/SSC <= 1+cutoff % continue search next set
        cid_o = cptloc;
    else
        segcount = c+2;
        cptloc = unique([cid_o N]);
        segidx = segidx(cptloc);
        return;
    end
end

% All change points should be elminated
segidx = numel(x);
segcount = 1;
end % end of function
%----------------------------------------------------------------%
function jc = allnr(n, r)
% Find all combinations of r from n
% Reference:
% J.Gentleman, Applied Statistics, Vol.24, No.3 (1975), pp.374-376
%
% nchoosek(1:n, r) should also work, but it errored when n is over 20.

loc = 1:r;
jc = loc;
nmr = n-r;

while ~(loc(1) == nmr+1)
    i = r;
    while loc(i) == nmr+i
        i = i-1;
    end
    loc(i) = loc(i) + 1;
    loc(i+1:r) = loc(i:r-1) + 1;
    jc(end+1, :) = loc; %#ok
end
end % end of function
%------------------------------------------------------%
function [wnum, wloc, wlens, wisodd] = getWindowParams(n, wsize)
% Return parameters of the window by its length and overlapping
% n - total data length
% wsize - window length
% wnum - number of windows
% wloc - the location of windows
% wisodd - if the window length is odd number

wnum = 1;
wloc = [1, n];
if  n~= wsize
    woverlap = 0.25; % window overlap, default 25%
    wnum = ceil((n - wsize)/((1-woverlap)*wsize)) +1;
    wloc = [fix(linspace(1, max(n-wsize,1), wnum)), n];
end
wlens = diff(wloc) + 1;
wisodd = mod(wlens, 2) ~= 0;
end % end of function

%-----------------------------------------------------------------------%
function y = smoothOutliers(y, R)
% Smooth procedure (see ref[1]):
% The smoothing region for each i is given by i-R,...,i,...i+R. Let mi be
% the median of the data in the smooth region, and sigma be the standard
% deviation of the entire data. If Xi is the max or min of all the
% observation in the smoothing region, find j in th smoothing region
% closest to it. If the distance from Xi to Xj exceeds L*sigma, replace Xi
% with mi+sign(Xi-X)*M*sigma.
% R = 2; % typical 2 to 5
% L = 4;
% M = 2;

L = 4;
M = 2;
sigma = std(y);

% Construct a array of smoothing regions (rows)
yLen = numel(y);
rLen = 2*R + 1;
smoothArray = zeros(yLen, rLen);

cy = [nan(R,1); y; nan(R,1)];
for i = 1:rLen
    smoothArray(:,i) = cy(i:yLen+i-1);
end

% Get medians of each smooth region
med_smoothArray = nanmedian(smoothArray, 2);

smoothArray(:,R+1)=[];
idx = y > (max(smoothArray, [], 2) + L*sigma);
y(idx) = med_smoothArray(idx) + M*sigma;

idx = y < (min(smoothArray, [], 2) - L*sigma);
y(idx) = med_smoothArray(idx) - M*sigma;
end % end of function
%----------------------------------------------------%
function valid = checkInputStruct(dataStruct)
% Checks for input data structure field names

valid =  isfield(dataStruct,'Chromosome')&& ...
    isfield(dataStruct,'GenomicPosition') && ...
    isfield(dataStruct,'Log2Ratio');
end % end of function
%----------------------------------------------------%
function plotSingleChromosomeSegments(x, y, starts, ends, means, chrom, sample)
figure;
plot(x, y, '.', 'color', [0.6 0.6 1])

ylims = [min(min(y), -1), max(max(y), 1)];
for i = 1:numel(starts)
    line([starts(i) ends(i)], [means(i) means(i)],...
        'Color', [1 0 0],...
        'LineWidth', 1.5);
end

title(sprintf('%s - Chr %d', sample, chrom), 'Interpreter', 'none')
xlabel('Genomic Position');
ylabel('Log2(Ratio)')
ylim(gca, ylims)
end % end of function

%------------------------------------------------------%
function plotAllChromosomeSub(X,Y, CID, chromIDs, segChromId, segStruct, sample)
% X - genomic positions
% Y - log2ratio
% CID - chromosome ids
% segStruct - final structure of segments
% sample - Sample name

n = numel(chromIDs);

if n <= 2
    spn = [1 n];
else
    spn = [ceil(n/5), 5];
end

y_lims = [min(min(Y), -1.5), max(max(Y), 1.5)];

figure;
for c =1:n
    chrom = chromIDs(c);
    idx = CID == chrom;
    x= X(idx);
    y =Y(idx);
    
    idx = find(ismember(segChromId, chrom));
    x_lims = [0 max(x)];
    subplot(spn(1), spn(2), c)
    hp = plot(x, y, '.','Color', [0.3 0.3 1]);
    line(x_lims, [0 0], 'Color', [0.8 0.8 0.8], 'Linewidth', 1.5); %Zero line
    if ~isempty(idx)
        sx = [segStruct(idx).Start';segStruct(idx).End'];
        sy = repmat(segStruct(idx).Mean', 2, 1);
        line(sx, sy, 'Color', [1 0 0]);
    end
    
    ha = get(hp, 'Parent');
    set(ha, 'xtick', [],...
        'box', 'on',...
        'xlim', double(x_lims),...
        'ylim', double(y_lims));
    title(sprintf('chr %d', chrom), 'Interpreter', 'none');
end
bioinfoprivate.suptitle(sprintf('%s',sample));
end % end of function


%------------------------------------------------------%
function plotAllChromosomes(X,Y, CID, chromIDs, segChromId, numSegs, segStruct, sample)
% X - genomic positions
% Y - log2ratio
% CID - chromosome ids
% segStruct - final structure of segments
% sample - Sample name

% Plot the entile genome
n = numel(chromIDs);
chr_endIdx = zeros(1, n);
chr_data_len = zeros(1,n);
for i = 1:n
    tmp = CID == chromIDs(i);
    chr_endIdx(i) = find(tmp, 1, 'last');
    chr_data_len(i) = length(find(tmp));
end

x_lims = [0 chr_endIdx(n)];
y_lims = [min(min(Y), -2), max(max(Y), 2)];

% Draw a vertical bar at the end of a chromosome to indicate the border
x_vbar = repmat(chr_endIdx, 2, 1);
y_vbar = repmat(y_lims', 1, n);

% Label the autosome with their chromosome numbers
x_label = chr_endIdx - ceil(chr_data_len/2);
y_label = zeros(1, length(x_label))+ (y_lims(1)-0.15);
chr_labels = cellstr(num2str(chromIDs));

% A gray zero line
x_zero = [X(1); X(end)];
y_zero = [0;0];

% Plot the segments
sn = sum(numSegs);
x_seg = zeros(2, sn);
y_seg = zeros(2, sn);
count = 1;
for i = 1:numel(segChromId)
    for j = 1:numSegs(i)
        x_seg(:, count) = [find(X==segStruct(i).Start(j)&CID ==segChromId(i),1);...
            find(X==segStruct(i).End(j)&CID ==segChromId(i),1)];
        y_seg(:, count) = [segStruct(i).Mean(j);segStruct(i).Mean(j)];
        
        count = count +1;
    end
end

figure; hold on
h_ratio = plot(Y, '.');
line(x_zero, y_zero, 'Color', [0.8 0.8 0.8], 'Linewidth', 2);
line(x_vbar, y_vbar, 'color', [0.8 0.8 0.8]);
line(x_seg, y_seg, 'Color', [1 0 0], 'Linewidth', 2);
text(x_label, y_label, chr_labels,...
    'Fontsize', 8, 'HorizontalAlignment', 'Center');

h_axis = get(h_ratio, 'parent');
set(h_axis, 'xtick', [], 'ygrid', 'on', 'box', 'on',...
    'xlim', x_lims, 'ylim', y_lims)

title(sample, 'Interpreter', 'none')
xlabel({'', 'Chromosome'})
ylabel('Log2(Ratio)')
hold off
end % end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha, nPerm,hybridFlag, smoothFlag, winSize, pruneFlag,gamma,...
    segSampleId, segChromId, plotId, plotType, verbose,useStoppingRule ] =...
    parse_inputs(numSample, chromIDs, varargin)
% Initialization
alpha = 0.01; % p-value limit
nPerm = 10000; % Number of permutations

% Hybrid p-value
hybridFlag = true;

% Smooth outlier
smoothFlag = true;

% Window size for Permutation method with data length longer than 200
% points
winSize = 200;

% Prune segments
pruneFlag = false;
gamma = 0.05; % allowed error squared sum increase

% verbose option
verbose = true;

% Stopping rule
useStoppingRule = false;

% Sample indices, chromosome ids to be segmented
segChromId = chromIDs;
sampleIDs = 1:numSample;
segSampleId = sampleIDs;

if nargout == 0
    plotId = 0;
else
    plotId = -1;
end

plotType = 1; %w =1, s=2

%--------Input checking ---------------------------------%
% get input arguments
if  nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:cghcbs:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'alpha','permutations','method', 'smooth','windowsize',...
        'prune', 'errsum', 'sampleindex', 'chromosome', 'showplot','verbose','stoppingrule'};
    for j=1:2:nargin-3
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:cghcbs:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:cghcbs:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % 'alpha'
                    if isnumeric(pval) && isscalar(pval) && pval < 1 && pval > 0
                        alpha =  pval;
                    else
                        error(message('bioinfo:cghcbs:IncorrectAlpha'));
                    end
                case 2 % permutations
                    if isnumeric(pval) && isscalar(pval) && pval >= 3
                        nPerm = fix(pval);
                    else
                        error(message('bioinfo:cghcbs:IncorrectNumberOfPermutations'));
                    end
                case 3 % method
                    if ischar(pval)
                        okmethods = {'perm', 'hybrid'};
                        nm = find(strncmpi(pval,okmethods,numel(pval)));
                        if isempty(nm)
                            error(message('bioinfo:cghcbs:UnknownMethodName'));
                        elseif length(nm) > 1
                            error(message('bioinfo:cghcbs:AmbiguousMethodName', pval));
                        else
                            hybridFlag = (nm==1);
                        end
                    else
                        error(message('bioinfo:cghcbs:MethodNameNotValid'));
                    end
                case 4 % smoothFlag
                    smoothFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 5 % window size
                    if (~isscalar(pval) || pval<=1)
                        error(message('bioinfo:cghcbs:IncorrectWindowSize'))
                    end
                    winSize = pval;
                case 6 % prune
                    pruneFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 7 % error sum of square
                    if ~isnumeric(pval) || ~isscalar(pval) || pval <= 0
                        error(message('bioinfo:cghcbs:ErrSumMustBeNumeric'));
                    elseif pval >= 1
                        error(message('bioinfo:cghcbs:ErrSumMustBeLessThan1'));
                    end
                    gamma = pval;
                case 8 % sample indices
                    if ~isnumeric(pval) || ~isvector(pval)
                        error(message('bioinfo:cghcbs:SampleIndicesNotNumericVector'));
                    else
                        segSampleId = unique(pval(:));
                        
                        if ~all(ismember(segSampleId, sampleIDs))
                            inidx = (ismember(segSampleId, sampleIDs))';
                            warning(message('bioinfo:cghcbs:InvalidSampleIndices', segSampleId( ~inidx )));
                            segSampleId = segSampleId(inidx);
                        end
                    end
                    
                case 9 % chromosomes
                    if ~isnumeric(pval) || ~isvector(pval)
                        error(message('bioinfo:cghcbs:ChromosomeIdNotNumericVector'));
                    else
                        segChromId = unique(pval(:));
                        
                        if ~all(ismember(segChromId, chromIDs))
                            inidx = (ismember(segChromId, chromIDs))';
                            warning(message('bioinfo:cghcbs:InvalidChromosomeId', segChromId( ~inidx )));
                            segChromId = segChromId(inidx);
                        end
                    end
                    
                case 10 % showplot
                    if bioinfoprivate.opttf(pval)
                        if isnumeric(pval)
                            if isscalar(pval)
                                plotId = double(pval);
                            else
                                plotId = 0;
                                warning(message('bioinfo:cghcbs:SPNoScalar'))
                            end
                        else
                            plotId = 0;
                        end
                    else
                        if ischar(pval)
                            types = {'w', 's'};
                            plotType = find(strncmpi(pval, types, 1));
                            
                            if isempty(plotType)
                                plotId = -1;
                                warning(message('bioinfo:cghcbs:SPNoValid'))
                            else
                                plotId = 0;
                            end
                        else
                            plotId = -1;
                        end
                    end
                case 11 % verbose flag
                    verbose = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 12 % stopping rule flag
                    useStoppingRule = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            end
        end
    end
end
end
