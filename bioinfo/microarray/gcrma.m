function expression = gcrma(pms, mms, indices, varargin)
%GCRMA performs GC Robust Multi-array Average (GCRMA) background
% adjustment, quantile normalization and median-polish summarization.
% 
%   EXPRESSION=GCRMA(PMS,MMS,INDICES,ALPHAPM,ALPHAMM) performs GCRMA
%   background adjustment with probe affinity information ALPHAPM and
%   ALPHMM, quantile normalization and RMA summarization. EXPRESSION are
%   the log2 based expression measure for each probe set. PMS are the
%   perfect match(PM) probe intensities. MMS are the mismatch(MM) probe
%   intensities. INDICES correspond to the probe indexing. ALPHAPM and
%   ALPHMM are affinities for PM and MM probes. 
% 
%   GCRMA(PMS,MMS,INDICES,SEQS) performs GCRMA background adjustment with
%   probe sequence information when PM and MM probe affinities are not
%   provided. The probe affinities will be computed from sequence
%   information SEQS.  SEQS is an Nx25 matrix of the probe (PM) sequence
%   bases in unsigned integers.
% 
%   Note: When there are not APM and AMM affinity information and sequence
%   information, the parameters for background correction are estimated
%   from PMS and MMS intensities.
% 
%   GCRMA(...,'CHIPINDEX', CID) computes affinities with sequence
%   information and mismatch probe intensities from the chip with column
%   index CID if there is not affinity APM and AMM information and need to
%   be computed from the sequence information SEQS. Default CID = 1. If the
%   affinity information is provided, this option is ignored.
% 
%   GCRMA(...,'OPTICALCORR', TF) performs optical background correction on
%   the input probe intensities if TF is set to true. Default is set to
%   TRUE.
% 
%   GCRMA(...,'CORRCONST', RHO) sets the correlation coefficient constant
%   RHO of the bivariate-normal distribution of the log non-specific
%   binding background intensity for a PM and MM probe pair. The value
%   should be between 0 and 1. Default RHO = 0.7.
% 
%   GCRMA(...,'METHOD', METHOD) sets an estimate method to estimate the
%   signal. METHOD can be the faster but ad hoc Maximum Likelihood Estimate
%   'MLE' (default) or 'EB' for the slower and more formal Empirical-Bayes
%   type estimator. 
%
%   GCRMA(...,'TUNINGPARAM', TP) sets a tuning parameter for the estimate
%   method. This tuning parameter sets the lower bound of signal values
%   with positive probability. The value must be greater than or equal to
%   0. Default is 5 for MLE or 0.5 for EB. See DOC GCRMA for more details.
% 
%   GCRMA(...,'GSBCORR', TF) performs gene specific binding (GSB)
%   correction with probe affinity data if TF is set to true. Default is
%   set to TRUE. If there is not probe affinity information this option is
%   ignored.
%
%   GCRMA(...,'NORMALIZE', TF) performs quantile normalization on
%   background adjusted data if TF is et to true. Default is set to TRUE.
% 
%   GCRMA(...,'VERBOSE', TF) turns off verbose output if TF is set to
%   FALSE. Default is TRUE.
%  
%   Example:
%    
%       load prostatecancerrawdata
% 
%       % Compute the probe affinities using mismatch probe intensity
%       [apm, amm] = affyprobeaffinities(seqMatrix, mmMatrix(:,1),...
%                          'ProbeIndices', probeIndices);
% 
%       % Compute probeset gene expression data
%       expdata = gcrma(pmMatrix, mmMatrix, probeIndices, apm, amm);
% 
%       % With probe sequence matrix
%       expdata = gcrma(pmMatrix, mmMatrix, probeIndices, seqMatrix);
% 
%   See also AFFYGCRMA, AFFYPREPROCESSDEMO, AFFYPROBEAFFINITIES,
%   AFFYPROBESEQREAD, AFFYRMA, CELINTENSITYREAD, GCRMABACKADJ, MAFDR,
%   MATTEST, QUANTILENORM, RMABACKADJ, RMASUMMARY.

%   GeneChip and Affymetrix are registered trademarks of Affymetrix, Inc.

%   Copyright 2006-2012 The MathWorks, Inc.


% References:
% [1] Zhijin Wu, Rafael A. Irizarry, Robert Gentleman, Francisco Martinez
%     Murillo, and Forrest Spencer, "A Model Based Background Adjustment
%     for Oligonucleotide Expression Arrays", J Amer Stat Assoc Vol.99,
%     No.468, 2004, pp.909-917.
% [2] Zhijin Wu, Rafael A. Irizarry. "Stochastic Models Inspired by
%     Hybridization Theory for Short Oligonucleotide Arrays", Proceedings
%     of RECOMB 2004. J Comput Biol. 2005 Jul-Aug;12(6):882-93.
% [3] Zhijin Wu and Rafael A. Irizarry, "A Statistical Framework for the
%     Analysis of Microarray Probe-Level Data", Johns Hopkins University,
%     Biostatistics Working Papers 73 (2005).
% [4] T.Speed, "Background models and GCRMA", Lecture 10, Statistics
%     246, UC Berkeley, Spring 2006
%     http://www.stat.berkeley.edu/users/terry/Classes/s246.2006/Week
%     10/Week10L1.pdf.

%   GCRMA(...,'ALPHA', A) sets the prior signal exponential distribution
%   mean alpha A. It is used for compute posterior mean for the EB method
%   in GCRMABACKADJ. Default is 1.  See DOC RMABACKADJ for more details.
% 
%   GCRMA(...,'STEPS', S) sets number of steps when dividing the signal
%   distribution in natural log space into steps. This value is used for
%   compute posterior mean for the EB method in GCRMABACKADJ. Default is
%   128. 
% 
%   GCRMA(...,'STRETCH', S) sets stretch correction constant. It can be
%   used for correcting results after MLE or EB background adjustments. S
%   can be a logical or a numerical number. No stretch correction If
%   S=FALSE (default). If S=TRUE, the correction factor is set to 1.15 for
%   MLE method and 1 for EB method. A numerical number will be the stretch
%   correction factor for the method used in background adjustment. 

% Validate input data
bioinfochecknargin(nargin,3,mfilename);

if ~isnumeric(pms) || ~isreal(pms)
   error(message('bioinfo:gcrma:PMProbeIntensityNotNumericAndReal')) 
end

if ~isnumeric(mms) || ~isreal(mms)
    error(message('bioinfo:gcrma:MMProbeIntensityNotNumericAndReal'))
end

if ~isequal(size(pms), size(mms))
    error(message('bioinfo:gcrma:IntensitiesOfPMMMSizeNotMatch'));
end

% Check probe indices are ok
if ~isa(indices, 'uint8') || ~isvector(indices)
    error(message('bioinfo:gcrma:InvalidProbeIndices'))
end

if isvector(pms)
    pms = pms(:);
    mms = mms(:);
end

nchips = size(pms,2);
nprobes = size(pms,1);

% initialization
apm = [];
amm = [];
seqs = [];
input_num = 3;

% determining required input variable
if nargin > 3
   if isvector(varargin{1}) && isnumeric(varargin{1})
       apm = varargin{1};
       if nargin >= 5 && isvector(varargin{2}) && isnumeric(varargin{2})
           amm = varargin{2};
           input_num = 5;
       else
           error(message('bioinfo:gcrma:NeedBothPMAndMMAffinities'))
       end
       
   elseif isa(varargin{1}, 'uint8')
       seqs = varargin{1};
       
       if size(seqs,1) ~= nprobes || size(seqs,2) ~= 25
           error(message('bioinfo:gcrma:SequenceMatrixWrongSize', nprobes));
       end
       input_num = 4;
   end       
end

if isempty(apm) || isempty(amm)
    useseqsflag = true;
    alphaId = 1;
else
    if ~isnumeric(apm) || ~isreal(apm) || ~isvector(apm)
        error(message('bioinfo:gcrma:PMAffinityNotNumericAndRealVector'))
    end

    if ~isnumeric(amm) || ~isreal(amm) || ~isvector(amm)
        error(message('bioinfo:gcrma:MMAffinityNotNumericAndRealVector'))
    end
    
    if nprobes ~= numel(apm) || nprobes ~= numel(amm)
        error(message('bioinfo:gcrma:IntensityandAffinityLengthNotMatch'));
    end
    
    % Columnize the data
    apm = apm(:);
    amm = amm(:);
    useseqsflag = false; 
    alphaId = 0;
end


% Initialize
opticalflag = true;
gsbflag = true;
normflag = true;
rho = 0.7;
method = 'mle';
tuningparam = [];
verbose = true;
alpha_eb = 1;
steps_eb = 128;
stretch = false;

% Deal with the various input
if nargin > input_num
    
    if input_num == 4 && rem(nargin, 2) == 1
        error(message('bioinfo:gcrma:IncorrectNumberOfArguments', mfilename));
    end
    
    if input_num ~=4 && rem(nargin, 2) == 0
        error(message('bioinfo:gcrma:IncorrectNumberOfArguments', mfilename));
    end
    
    okargs = {'chipindex','opticalcorr','corrconst','method',...
                'tuningparam','gsbcorr', 'normalize', 'verbose',...
                'alpha', 'steps','stretch'};
    for j = 1:2:nargin-input_num
        pname = varargin{j+(input_num-3)};
        pval = varargin{j+(input_num-3)+1};
        k = find(strncmpi(pname, okargs, numel(pname)));

        if isempty(k)
            error(message('bioinfo:gcrma:UnknownParameterName', pname));
        elseif length(k) > 1
            error(message('bioinfo:gcrma:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % chipindex
                    if isnumeric(pval)
                        if isscalar(pval)
                            alphaId = double(pval);
                        else
                            alphaId = 1;
                            warning(message('bioinfo:gcrma:CIDNotScalar'))
                        end
                    else
                        alphaId = 1;
                        warning(message('bioinfo:gcrma:CIDNotANumber'))
                    end
                    
                    if alphaId > nchips || alphaId <= 0
                        error(message('bioinfo:gcrma:NotExistChipIndex'))
                    end  
                case 2 % Optical correction
                   opticalflag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 3 % Correlation coefficient constant
                   if ~isnumeric(pval)
                        error(message('bioinfo:gcrma:CorrContMustBeNumeric'));
                    elseif pval < 0 || pval > 1
                        error(message('bioinfo:gcrma:CorrContMustBe0To1'));
                    end

                    rho = pval;               
                case 4 % estimate method
                    if ischar(pval)
                        okmethods = {'mle', 'eb'};                    
                        nm = find(strncmpi(pval, okmethods, numel(pval))); %#ok
                                          
                        if isempty(nm)
                            error(message('bioinfo:gcrma:UnknownMethodName'));
                        else
                            method = pval;
                        end
                        
                    else
                        error(message('bioinfo:gcrma:MethodNameNotValid'));
                    end    
                case 5 % turning parameter
                    if ~isnumeric(pval) || ~isscalar(pval)
                        error(message('bioinfo:gcrma:TuningParamNotSingleNumericValue'));
                    elseif (pval < 0 )
                        error(message('bioinfo:gcrma:badTuningParameterValue'));
                    else
                        tuningparam = pval;
                    end
                case 6 % GSB correction
                     gsbflag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 7 % normalize
                    normflag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 8 % verbose flag
                    verbose = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                 case 9 % Alpha - EB
                   if ~isnumeric(pval)
                        error(message('bioinfo:gcrma:AlphaMustBeNumeric'));
                    elseif pval < 0
                        error(message('bioinfo:gcrma:AlphaMustBeGreaterThan0'));
                    end

                    alpha_eb = pval;
                case 10 % Steps - EB
                   if ~isnumeric(pval)
                        error(message('bioinfo:gcrma:StepsMustBeNumeric'));
                    elseif pval < 1
                        error(message('bioinfo:gcrma:StepsMustBeGreaterThan1'));
                    end

                    steps_eb = pval;  
                 case 11 % stretch
                     stretch = pval;
            end % switch
        end % if
    end % for loop
end % if

if isempty(seqs) && useseqsflag
    gsbflag = false;
end

% Compute affinities for PM and MM probes
if useseqsflag    
    try
        if verbose
            disp('Calculating probe affinities.');
        end
        [apm, amm] = affyprobeaffinities(seqs, mms(:,alphaId), 'probeindice', indices);
    catch theException
        gsbflag = false;
        warning(message('bioinfo:gcrma:UnableComputeProbeAffinity', theException.message));
    end
end

% GCRMA background adjustment
try
    pms_adj = gcrmabackadj(pms, mms, apm, amm,...
        'opticalcorr', opticalflag,...
        'corrconst', rho,...
        'method', method,...
        'tuningparam', tuningparam,...
        'gsbcorr', gsbflag,...
        'verbose', verbose,...
        'alpha', alpha_eb,...
        'steps', steps_eb,...
        'stretch', stretch);
catch theException
    % catch error throw by gcrmabackadj and rethrow 
    msgId = 'bioinfo:gcrma:gcrmabackadjError';
    newException = MException(msgId,getString(message(msgId)));
    throw(addCause(newException,theException))
 end

% Normalize with quantilenorm and summarize with rmasummary
if normflag
    try 
        if verbose
            disp('Normalizing.');
        end
        pms_adj = quantilenorm(pms_adj);
    catch theException
        % catch error throw by quantilenorm and rethrow 
        msgId = 'bioinfo:gcrma:quantilenormError';
        newException = MException(msgId,getString(message(msgId)));
        throw(addCause(newException,theException))
    end
end

try 
    if verbose
        disp('Calculating expression.');
    end
    expression = rmasummary(indices, pms_adj);
catch theException
    % catch error throw by rmasummary and rethrow 
    msgId = 'bioinfo:gcrma:rmasummaryError';
    newException = MException(msgId,getString(message(msgId)));
    throw(addCause(newException,theException))    
end
