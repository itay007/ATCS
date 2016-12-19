function [affpm, affmm, A, stats] = affyprobeaffinities(seqMatrix, mms, varargin)
%AFFYPROBEAFFINITIES computes Affymetrix probe affinities from their
% sequences and mismatch probe intensities from a non-specific binding
% experiment.
% 
%   [ALPHAPM, ALPHAMM] = AFFYPROBEAFFINITIES(SEQMATRIX, MMS) returns
%   affinities ALPHAPM for PM probes and the affinities ALPHAMM for MM
%   probes from a non-specific binding (NSB) experiment. SEQMATRIX is an
%   Nx25 matrix of the probe (PM) sequence bases in number, with rows for
%   the probes in a chip and every column is the position of 25. MMS are
%   the MM probe intensities of a chip. The affinity of a probe is the sum
%   of position-dependent base affinities. NaN is returned for probes with
%   no sequence information.
% 
%   [ALPHAPM, ALPHAMM, A] = AFFYPROBEAFFINITIES(...) returns the base
%   affinity coefficients A of a polynomial of degree 3. For a given type
%   of base, the positional effect is modeled as a polynomial of degree 3.
% 
%   [ALPHAPM, ALPHAMM, A, STATS] = AFFYPROBEAFFINITIES(...) returns a
%   regression statistics vector STATS containing, in the following order,
%   the R-square statistic, the F statistic and p-value, and an estimate of
%   the error variance of the estimate of base affinity coefficients A. 
% 
%   AFFYPROBEAFFINITIES(...,'PROBEINDICES', INDICES) uses a probe indices
%   vector to normalize the probe intensities with the median of their
%   probe set intensities. 
% 
%   Note: This option is recommended only to compute probe affinities from
%   data that are not from a non-specific binding (NSB) experiment. See DOC
%   AFFYPROBEAFFINITIES for more details.
% 
%   AFFYPROBEAFFINITIES(...,'SHOWPLOT',TF) displays the probe affinity base
%   profile if TF is TRUE. If AFFYPROBEAFFINITIES is called without output
%   arguments, a plot will be shown unless TF is FALSE. 
% 
%   Example:
% 
%       load prostatecancerrawdata
%       
%       % Since this dataset is not from an NSB experiment, the
%       % 'PROBEINDICES' option is used to compute the probe affinities.  
%       [apm, amm] = affyprobeaffinities(seqMatrix, mmMatrix(:,1),...
%                          'ProbeIndices', probeIndices, 'showplot', true);
% 
%   See also AFFYGCRMA, AFFYPREPROCESSDEMO, AFFYPROBESEQREAD,
%   CELINTENSITYREAD, GCRMA, GCRMABACKADJ, REGRESS.
% 
%   Affymetrix is a registered trademark of Affymetrix, Inc.

% Copyright 2006-2008 The MathWorks, Inc.


% References: 
% [1] Naef F, Magnasco MO. "Solving the Riddle of the Bright Mismatches:
%     Labeling and Effective Binding in Oligonucleotide Arrays" Physical
%     Review E 68, 011906 (2003).
% [2] Zhijin Wu, Rafael A. Irizarry, Robert Gentleman, Francisco
%     Martinez Murillo, and Forrest Spencer, "A Model Based Background
%     Adjustment for Oligonucleotide Expression Arrays", J Amer Stat Assoc
%     Vol.99, No.468, 2004, pp.909-917.

bioinfochecknargin(nargin,2,mfilename);

indices = [];

if ~isnumeric(mms) || ~isreal(mms) || ~isvector(mms)
   error(message('bioinfo:affyprobeaffinities:ProbeIntensityNotNumericAndRealVector')) 
end

if ~isa(seqMatrix, 'uint8')
    error(message('bioinfo:affyprobeaffinities:SequenceMatrixNotUINT8')) 
end

N = numel(mms);

if size(seqMatrix,1) ~= N || size(seqMatrix,2) ~= 25
    error(message('bioinfo:affyprobeaffinities:SequenceMatrixWrongSize', N)); 
end

% Columnize the input if is a vector
mms = mms(:);

% Initialize
if nargout > 0
    showplot = false;
else
    showplot = true;
end

% Deal with the various input
if nargin > 2
    if rem(nargin, 2) == 1
        error(message('bioinfo:affyprobeaffinities:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'probeindices', 'showplot'};
    for j = 1:2:nargin-3
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs, numel(pname)));
        
        if isempty(k)
            error(message('bioinfo:affyprobeaffinities:UnknownParameterName', pname));
        elseif length(k) > 1
            error(message('bioinfo:affyprobeaffinities:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % probe indices
                    if ~isa(pval, 'uint8')
                        error(message('bioinfo:affyprobeaffinities:InvalidProbeIndices'))
                    end
                    
                    if numel(pval) ~= N
                        error(message('bioinfo:affyprobeaffinities:SizeProbeIndices'));
                    end
                    indices = pval(:);
                case 2 % showplot
                    showplot = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            end % switch
        end % if
    end % for loop
end % if

if nargout == 0 && ~showplot
    return;
end

% Normalize probe intensities with the probeset median intensity if indices are
% provided
if ~isempty(indices)
    psIndices = find(indices == 0); %Indices of probeset.Probe pair always start at 0
    nProbeSets = numel(psIndices); %number of probe set
    nProbePairs = 0; %#ok Number of probe pair in a probe set %#ok

    Y = zeros(N, 1);
    psIndices = [psIndices; N+1];
    
    for i = 1 : nProbeSets
        nProbePairs = psIndices(i+1) - psIndices(i);

        ppstart = psIndices(i); % start probe pair index in a probe set
        ppend = ppstart+nProbePairs-1; % end prope pair index in a probe set

        ppvalues = mms(ppstart:ppend);
        ps_median = median(ppvalues);
        Y(ppstart:ppend) = ppvalues /ps_median;
    end
    mms = Y;
end

% Filter out the probes that does not have sequences
noseq_idx = seqMatrix(:,1)==0;
affpm = nan(N,1);
affmm = nan(N,1);

seqMatrix_seq = seqMatrix(~noseq_idx, :);
intensities_seq = mms(~noseq_idx);

% Using MM intensities
mmseqMatrix = seqMatrix_seq;
mmseqMatrix(:, 13) = seqcomplement(seqMatrix_seq(:, 13));

if nargout == 0 && showplot
    affyprobebaseprofile(mmseqMatrix, intensities_seq, 'showplot', showplot);
    return;
end

% Get base profile coefficients
[A, stats, mmX] = affyprobebaseprofile( mmseqMatrix, intensities_seq, 'showplot', showplot);

%---------Compute affinities ---------------------------
% Find the affinities at position 13 for all bases, T13 is set to 0
k = 13;
poly13 = [ones(size(k)) k k.^2 k.^3];
alpha_A_13 = poly13*A(2:5);
alpha_C_13 = poly13*A(6:9);
alpha_G_13 = poly13*A(10:13);
alpha_T_13 = 0;

affmm(~noseq_idx) = mmX * A(2:13);
affpm(~noseq_idx) = affmm(~noseq_idx);

idx = seqMatrix(:, 13) == 1; % base A
affpm(idx) = affmm(idx) + alpha_T_13 - alpha_A_13;

idx = seqMatrix(:, 13) == 2; % base C
affpm(idx) = affmm(idx) + alpha_G_13 - alpha_C_13;

idx = seqMatrix(:, 13) == 3; % base G
affpm(idx) = affmm(idx) + alpha_C_13 - alpha_G_13;

idx = seqMatrix(:, 13) == 4; % base T
affpm(idx) = affmm(idx) + alpha_A_13 - alpha_T_13;
