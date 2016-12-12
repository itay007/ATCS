function [pms_adj, nsbstruct] = gcrmabackadj(pms, mms, apm, amm, varargin)
%GCRMABACKADJ performs GC Robust Multi-array Average (GCRMA) background
% adjustment for Affymetrix GeneChip probe-level data.
% 
%   The GCRMA background model is described briefly as, for any particular
%   probe-pair,
%       PM = Opm + Npm + S; MM = Omm + Nmm;  
%   O represents optical noise and is an array-dependent constant, N
%   represents non-specific binding (NSB) noise, and S is a quantity
%   proportional to RNA expression (signal). log(Npm) and log(Nmm) follow a
%   bivariate-normal distribution with correlation coefficient constant RHO
%   across probes. Maximum Likelihood estimate (MLE) and Empirical Bayes
%   (EB) approaches are the two estimators used to estimate the signal S. 
% 
%   PMS_ADJ=GCRMABACKADJ(PMS,MMS,APM,AMM) performs GCRMA background
%   adjustment with probe sequence information. PMS_ADJ are the background
%   adjusted perfect match(PM) probe intensities PMS. MMS are the
%   mismatch(MM) probe intensities. PMS and MMS can be matrices with
%   several chips of the same type, each column corresponds with each chip,
%   and all have the same probe indices. APM and AMM are the probe
%   affinities for PM and MM probes with known sequences. 
% 
%   Note: When there is not APM and AMM affinity information, the
%   parameters to estimate NSB background are computed from PM and MM
%   intensities.
% 
%   [PMS_ADJ, NSBSTRUCT]=GCRMABACKADJ(PMS,MMS,APM,AMM) returns non-specific
%   binding (NSB) parameters used for NSB background correction for all the
%   chips. The fields of NSBSTRUCT are: sigma, and median of mu of the
%   bivariate-normal distribution of NSB.   
% 
%   GCRMABACKADJ(...,'OPTICALCORR', TF) performs optical background
%   correction on the input probe intensities if TF is set to true. Default
%   is set to TRUE.
% 
%   GCRMABACKADJ(...,'CORRCONST', RHO) sets the correlation coefficient
%   constant RHO of the bivariate-normal distribution of the log
%   non-specific binding background intensity for a PM and MM probe pair.
%   The value should be between 0 and 1. Default RHO = 0.7.
% 
%   GCRMABACKADJ(...,'METHOD', METHOD) sets an estimate method to estimate
%   the signal. METHOD can be the faster but ad hoc Maximum Likelihood
%   Estimate 'MLE' (default) or 'EB' for the slower and more formal
%   Empirical-Bayes-type estimator. 
%
%   GCRMABACKADJ(...,'TUNINGPARAM', TP) sets a tuning parameter for the
%   estimate method. This tuning parameter sets the lower bound of signal
%   values with positive probability. The value must be greater than or
%   equal to 0. Default is 5 for MLE or 0.5 for EB. See DOC GCRMABACKADJ
%   for more details.
% 
%   GCRMABACKADJ(...,'ADDVARIANCE', TF) adds the signal variance to the
%   weight function for smoothing low signal edge when set true. Default is
%   set to FALSE.
%   
%   GCRMABACKADJ(...,'GSBCORR', TF) performs gene specific binding (GSB)
%   correction with probe affinity data if TF is set to true. Default is
%   set to FALSE. If there is not probe affinity information this option
%   is ignored.
% 
%   GCRMABACKADJ(...,'SHOWPLOT',SP) displays optical adjusted log2 based MM
%   intensities vs. its affinities and LOWESS fit for computing NSB data of
%   the chip with column index SP. When SP is TRUE, the plot of the first
%   chip in PMS is displayed. If GCRMABACKADJ is called without output
%   arguments, a plot will be shown unless SP is FALSE. In the case with
%   no affinity information, this option is ignored.
% 
%   GCRMABACKADJ(...,'VERBOSE',false) turns off verbose output. Default is
%   true.
%  
%   Example:
% 
%       load prostatecancerrawdata
% 
%       % Compute the probe affinities using mismatch probe intensity
%       [apm, amm] = affyprobeaffinities(seqMatrix, mmMatrix(:,1),...
%                          'ProbeIndices', probeIndices);
% 
%       % Do background adjustment with MLE estimator
%       pms_adj = gcrmabackadj(pmMatrix, mmMatrix, apm, amm,'showplot',3);
%
%       % Adjust background with EB estimator (slower)
%       pms_adj = gcrmabackadj(pmMatrix, mmMatrix, apm, amm,'method','EB');
% 
%   See also AFFYGCRMA, AFFYPREPROCESSDEMO, AFFYPROBEAFFINITIES,
%   AFFYPROBESEQREAD, AFFYRMA, CELINTENSITYREAD, GCRMA, QUANTILENORM,
%   RMABACKADJ, RMASUMMARY.

%   Copyright 2006-2009 The MathWorks, Inc.


% References: 
% [1] Zhijin Wu, Rafael A. Irizarry, Robert Gentleman, Francisco
%     Martinez Murillo, and Forrest Spencer, "A Model Based Background
%     Adjustment for Oligonucleotide Expression Arrays", J Amer Stat Assoc
%     Vol.99, No.468, 2004, pp.909-917.
% [2] Zhijin Wu, Rafael A. Irizarry. "Stochastic Models Inspired by
%     Hybridization Theory for Short Oligonucleotide Arrays", Proceedings
%     of RECOMB 2004. J Comput Biol. 2005 Jul-Aug;12(6):882-93.
% [3] Zhijin Wu and Rafael A. Irizarry, "A Statistical Framework for the
%     Analysis of Microarray Probe-Level Data", Johns Hopkins University,
%     Biostatistics Working Papers 73 (2005).
% [4] Zhijin Wu, Rafael A. Irizarry. "A Model Based Background Adjustment
%     for Oligonucleotide Expression Arrays", RSS workshop on Gene
%     expression, Wye, England (2003).
% [5] T.Speed, "Background models and GCRMA", Lecture 10, Statistics
%     246, UC Berkeley, Spring 2006
%     http://www.stat.berkeley.edu/users/terry/Classes/s246.2006/Week
%     10/Week10L1.pdf.

%   GCRMABACKADJ(...,'ALPHA', A) sets the prior signal exponential
%   distribution mean alpha A. It is used for compute posterior mean for
%   the EB method. Default is 1.  See DOC RMABACKADJ for more details.
% 
%   GCRMABACKADJ(...,'STEPS', S) sets number of steps when dividing the
%   signal distribution in natural log space into steps. This value is used
%   for compute posterior mean for the EB method. Default is 128. 
% 
%   GCRMABACKADJ(...,'STRETCH', S) sets stretch correction constant. It can
%   be used for correcting results after MLE or EB background adjustments.
%   S can be a logical or a numerical number. No stretch correction If
%   S=FALSE (default). If S=TRUE, the correction factor is set to 1.15 for
%   MLE method and 1 for EB method. A numerical number will be the stretch
%   correction factor for the method used in background adjustment. 
%
%   Note: The stretch correction is done to the output from fast algorithm
%   (MLE method) in the gcrma package in Bioconductor by default. No
%   reference was given.

bioinfochecknargin(nargin,4,mfilename);

if ~isnumeric(pms) || ~isreal(pms)
   error(message('bioinfo:gcrmabackadj:PMProbeIntensityNotNumericAndReal')) 
end

if ~isnumeric(mms) || ~isreal(mms)
    error(message('bioinfo:gcrmabackadj:MMProbeIntensityNotNumericAndReal'))
end

if ~isequal(size(pms), size(mms))
    error(message('bioinfo:gcrmabackadj:IntensitiesOfPMMMSizeNotMatch'));
end

if isvector(pms)
    pms = pms(:);
    mms = mms(:);
end

nchips = size(pms,2);
nprobes = size(pms,1);

if isempty(apm) || isempty(amm)
    mmadjustflag = true;
else
    if ~isnumeric(apm) || ~isreal(apm) || ~isvector(apm)
        error(message('bioinfo:gcrmabackadj:PMAffinityNotNumericAndRealVector'))
    end


    if ~isnumeric(amm) || ~isreal(amm) || ~isvector(amm)
        error(message('bioinfo:gcrmabackadj:MMAffinityNotNumericAndRealVector'))
    end
    
    if nprobes ~= numel(apm) || nprobes ~= numel(amm)
        error(message('bioinfo:gcrmabackadj:IntensityandAffinityLengthNotMatch'));
    end
    
    % Columnize the data
    apm = apm(:);
    amm = amm(:);
    mmadjustflag = false; % using only mm intensities for estimate parameters
    
    % Now with pms, mms and affinities for pm and mm probes.
    % Find the indices that apm and amm are not NANs.
    alpha_idx = ~isnan(amm);
end

rho = 0.7; % See Ref[2]
% Tuning parameter for MLE and EB
tuningparam = [];
tparam_mle = 5; % See Ref.[3-4]
tparam_eb = 0.5; % See Ref.[3-4]
opticalflag = true;
addvarflag = false;
gsbflag = true;
verbose = true;
alpha_eb = 1;
steps_eb = 128;
stretch = false;

if mmadjustflag
    mleflag = false;
else
    mleflag = true;
end

if nargout > 0
    plotId = 0; 
else
    plotId = 1;
end

% Deal with the various input
if nargin > 4
    if rem(nargin, 2) == 1
        error(message('bioinfo:gcrmabackadj:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'corrconst', 'method','tuningparam', 'opticalcorr',...
              'addvariance', 'gsbcorr', 'showplot', 'verbose',...
              'alpha', 'steps', 'stretch'};
            
    for j = 1:2:nargin-4
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs, numel(pname)));

        if isempty(k)
            error(message('bioinfo:gcrmabackadj:UnknownParameterName', pname));
        elseif length(k) > 1
            error(message('bioinfo:gcrmabackadj:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % Correlation coefficient constant
                   if ~isnumeric(pval)
                        error(message('bioinfo:gcrmabackadj:CorrContMustBeNumeric'));
                    elseif pval < 0 || pval > 1
                        error(message('bioinfo:gcrmabackadj:CorrContMustBe0To1'));
                    end

                    rho = pval;               
                case 2 % estimate method
                    if ischar(pval)
                        okmethods = {'mle', 'eb'};
                        nm = strmatch(lower(pval), okmethods);
                        if isempty(nm)
                            error(message('bioinfo:gcrmabackadj:UnknownMethodName'));
                        elseif length(nm) > 1
                            error(message('bioinfo:gcrmabackadj:AmbiguousMethodName', pval));
                        else
                            mleflag = (nm==1);
                        end
                        
                    else
                        error(message('bioinfo:gcrmabackadj:MethodNameNotValid'));
                    end    
                case 3 % turning parameter
                    if ~isempty(pval)
                        if ~isnumeric(pval) || ~isscalar(pval)
                            error(message('bioinfo:gcrmabackadj:TuningParamNotSingleNumericValue'));
                        elseif (pval < 0 )
                            error(message('bioinfo:gcrmabackadj:badTuningParameterValue'));
                        else
                            tuningparam = pval;
                        end
                    else
                        tuningparam = pval;
                    end
                case 4 % Optical correction
                    opticalflag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 5 % Add variance to weight function
                    addvarflag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 6 % GSB correction
                     gsbflag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 7 % showplot
                    if bioinfoprivate.opttf(pval)
                        if isnumeric(pval)
                            if isscalar(pval)
                                plotId = double(pval);
                            else
                                plotId = 1;
                                warning(message('bioinfo:gcrmabackadj:SPNoScalar'))
                            end
                        else
                            plotId = 1;
                        end
                    else
                        plotId = 0;
                    end
                case 8 % verbose flag
                    verbose = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 9 % Alpha - EB
                   if ~isnumeric(pval)
                        error(message('bioinfo:gcrmabackadj:AlphaMustBeNumeric'));
                    elseif pval < 0
                        error(message('bioinfo:gcrmabackadj:AlphaMustBeGreaterThan0'));
                    end

                    alpha_eb = pval;
                case 10 % Steps - EB
                   if ~isnumeric(pval)
                        error(message('bioinfo:gcrmabackadj:StepsMustBeNumeric'));
                    elseif pval < 1
                        error(message('bioinfo:gcrmabackadj:StepsMustBeGreaterThan1'));
                    end

                    steps_eb = pval;
                case 11 % stretch
                    if bioinfoprivate.opttf(pval)
                        if isnumeric(pval)
                            if isscalar(pval)
                                stretch = double(pval);
                            else
                                error(message('bioinfo:gcrmabackadj:StretchNotScalar'))
                            end
                        else
                            stretch = true;
                        end
                    else
                        stretch = false;
                    end
                  
            end % switch
        end % if
    end % for loop
end % if

% Plot without any output
if nargout == 0 && plotId == 0
    return;
end

% Tuning parameter 
if isempty(tuningparam) 
    if mleflag
        tuningparam = tparam_mle;
    else
        tuningparam = tparam_eb;
    end
end

% Optical correct the PM, and MM intensities
if opticalflag
    [pms, mms] = opticalbackadj(pms, mms);
end

if nargout == 0 && plotId > 0 && plotId <= nchips && ~mmadjustflag
    gcrmaestimatensb( mms(alpha_idx, plotId), amm(alpha_idx), apm(alpha_idx), 'showplot', true);
    return;
end

if verbose
    if mleflag
        method = 'MLE';
    else
        method = 'EB';
    end
end

if islogical(stretch) && stretch
    if mleflag
        stretch = 1.15;
    else
        stretch = 1;
    end
end

if nargout == 2
    sigmas = zeros(1,nchips);
    medianmu = zeros(1, nchips);
end

%== GSB correction getting the linear regression coefficient.
if gsbflag
    fitval = estimategsb(pms(alpha_idx, :), apm(alpha_idx));
end

%== Background correction for each chip
varS = [];
pms_adj = pms;
nsbstruct = [];

for n = 1:nchips
    if verbose
        fprintf('Adjusting background for chip # %d of %d using %s method.\n',...
            n, nchips, method);
    end
    
    %== Estimate NSB parameters
    if mmadjustflag % with no affinity information,
        gsbflag = false;
        mu = log(mms(:,n)); % mu is computed directly from mms
        % Find the ones that is NSB
        idx = pms(:,n) < mms(:,n);
        sigma = sqrt(mean((log(pms(idx, n)) - mu(idx)).^2));

        if addvarflag && mleflag
            varS = exp(2*mu + sigma^2) .* (exp(sigma^2) - 1);
        end
    else
        % Get NSB parameters from MM intensities, mm and pm affinities
        mu_pm = zeros(numel(pms(:,n)), 1);
        mu_mm = mu_pm;
        showplot = false;
        if plotId == n
            showplot = true;
        end

       [sigma_0, mu_mm(alpha_idx), mu_pm(alpha_idx)] =gcrmaestimatensb(...
           mms(alpha_idx, n), amm(alpha_idx), apm(alpha_idx), 'showplot', showplot);

        % Fill in the mu for which don't have affinities
        mu_pm(~alpha_idx) = median(mu_pm(alpha_idx));
        mu_mm(~alpha_idx) = median(mu_mm(alpha_idx));
        
        mu = mu_pm + rho*(log(mms(:, n))-mu_mm);
        sigma = sqrt((1-rho^2))*(sigma_0);
     
        if addvarflag && mleflag
            varS = exp(2*mu_pm + sigma_0.^2) * (exp(sigma_0^2) - exp(sigma_0^2 * rho^2));
        end
    end
    
    if mleflag
        pms_adj(:,n) = mlestimator(pms(:,n), mu, sigma, tuningparam, varS);
    else
        pms_adj(:,n) = ebestimator(pms(:,n), mu, sigma, tuningparam, alpha_eb, steps_eb);
    end
    
    %== GSB adjustment, one chip at a time    
    if gsbflag
        pms_log2 = log2(pms_adj(alpha_idx,n));
        pms_adj(alpha_idx,n) = 2.^(pms_log2 - fitval(2)*apm(alpha_idx) + mean(fitval(2)*apm(alpha_idx)));
    end
    
    %== Stretch for fast algorithm - this is commented out because no ref.
    if ~islogical(stretch) && stretch ~= 1
        mu_s = mean(log(pms_adj(:,n)));
        pms_adj(:,n) = exp(mu_s + stretch*(log(pms_adj(:, n)) - mu_s));
    end
    
    %== output
    if nargout == 2
        sigmas(n) = sigma;
        medianmu(n) = median(mu);
    end
end

if nargout == 2
    nsbstruct.sigma = sigmas;
    nsbstruct.median_mu = medianmu;
end

%--------------------Helper functions ------------------------
function [pm, mm] = opticalbackadj(pm, mm)
% Optical background correction. O is a constant and is estimated with the
% minimum intensity observed on each array. Subtract 1 from O to avoid
% negatives when adjusting by subtracting O from PM and MM intensities. 
min_const = 1;
min_tmp = min(min(pm), min(mm)) - min_const;
min_int = repmat(min_tmp, size(pm,1), 1);
pm = pm - min_int;
mm = mm - min_int;

%--------------------------------------------------------------------------
function S1 = mlestimator(S1, mu, sigma, tuningparam, varS)
% MLE of Signal is PM-N* if PM-N* > tuningparam, tuningparam otherwise.
% Compute S1 expected values given S2, and N2. N1 and N2 follow the
% bivariate lognormal distribution. Var(log(N1))=Var(log(N2))=sigma^2. The
% correlation rho is a constant.

% References:
% [1] N.A.Abd Rabbo, H. M. Barakat "Estimation Problems in Bivariate
% Lognormal Distribution". Indian J. pure appl. Math., 10(7): 815-825, July
% 1979.

% The expected N is
N1 = exp(mu + (sigma^2)/2); 

% MLE of S for S >0
S = S1 - N1;

% Work S > 0 points only
idx = S > 0;
S(~idx) = 0; % to avoid Inf in log

% For microarray experiment it is important to estimate the fold-change of
% expression. It is better to minimize the mean square error (MSE):
%  E[{log(S*/S)}^2 l{S>0}|PM, MM]
% For log ratio we should exclude cases for which S = 0. The tuning
% parameter tuningparam is added to S. Adjust tuningparam to find better
% estimation. By default it is set to 5 for most of the Affymetrix arrays.
% A percentile weight function is used to smooth out the edge for small S
% after adding the tuning parameter. Let x = (S+tuningparam)/tuningparam,
% w(x) = xlog(x)/(xlog(x)+1). w(x) monotonously increase with x. If varflag
% is set, then incorporate var(S) in the weight function.

% Create a weight coefficient.
w = zeros(numel(S),1);

% S plus tuningparam
x = (S(idx) + tuningparam)/tuningparam;
xlogx = x.*log(x);
if ~isempty(varS)
    w(idx) = xlogx ./((xlogx + 1) + varS(idx) ./((tuningparam^2).*(xlogx +1)));
else
    w(idx) = xlogx ./(xlogx + 1);
end

S1 = exp(w.*log(S+tuningparam) + (1-w).*log(tuningparam));

%********************************************************
function S = ebestimator(S, mu, sigma, tuningparam, alpha, steps)
% Empirical Bayes Estimate - Consider S as a random variable,and minimizing 
%     E[{log(S'/S)}^2|S>0, PM, MM]
% The solution is the posterio mean estimate: s' = E[s|S>0, PM, MM], with
% s=log(S). A prior distribution of f(s)=1/(PM-exp(b)). mu = rho*(log(MM) -
% MUmm) + MUpm, and sigma* = sqrt(1-rho^2)*sigma. The value tuningparam is
% the smallest value of S with positive probability.

idx = S > tuningparam;
S_post = postmean(S(idx),  mu(idx), sigma, tuningparam, alpha, steps);

S(idx) = S_post;
S(~idx) = tuningparam; 

%-------------------------Helper Functions -------------------------------
function y = postmean(S, mu, sigma, lowerbnd, alpha, steps)
% Compute the posterior mean for log(S) with prior distribution of an
% exponential with rate a lowerbnd - smallest values of S with positive
% probability. The prior distribution of S is exponential with mean alpha. 

% Initialize
% Imposing a uniform distribution U[lowerbnd, log(2^16)] for a heavy-tailed
% prior distribution of S. Numbers based on optical scanner properties.
% Divided the distrubition into small bins. The bin width is set to an
% constant.
step_width = log(2^16)/steps;

% We are integrating from lowerbnd to S in exponentially-spaced steps.
log_lowerbnd = log(lowerbnd);
% The number of steps for each S
N = ceil((log(S)-log_lowerbnd)/step_width);

% Find the largest steps
[maxS, maxid] = max(S);
log_maxS = log(maxS);

% Largest number of steps or bins
bvec = exp(linspace(log_lowerbnd, log_maxS, N(maxid)));
bvec(1) = lowerbnd;
bvec(N(maxid)) = maxS;

% Numerically integrate over the steps for each S. See Ref.[4].
sumT = zeros(size(S));
sumW = zeros(size(S));
b = zeros(size(S));
normdn = zeros(size(S));
normdn0 = zeros(size(S));

% A constant for the erfc function
C = 1/(sqrt(2)*sigma);

% The first bin
b0 = bvec(1)*ones(size(S));
idx = S - b0 >=0;
z = (mu(idx) - log(S(idx)-b0(idx)))*C;
normdn0(idx) = 0.5*erfc(z);

for i = 2:numel(bvec)
    % smaller than N
    idx = (i < N);
    b(idx) = bvec(i);
    
    % for N
    idx = (N == i);
    b(idx) = S(idx); 
       
    % Cumulate S >= b only
    idx = ((S-b) >= 0);
    z = (mu(idx) - log(S(idx)-b(idx)))*C;
    normdn(idx) = 0.5*erfc(z);

    delW = (normdn - normdn0) .* (b0 .^(-alpha) + b .^(-alpha))./2;
    
    sumW = sumW + delW;
    sumT = sumT + delW .*(log((b0 + b)./2));
    
    normdn0 = normdn;
    b0 = b;
end

% Return the intensity in natural scale
y = exp(sumT ./ sumW);

%----------------------------------------------------------
function fitval=estimategsb(pms, apm)
% Fit a stright line through PM vs PM affinities plot
nrows = size(pms,1);
y = log2(pms(:));
% Create a RandStream
rstream = RandStream('shr3cong', 'Seed', 1);
M = numel(y);
N = min(M, 25000);
rp = randperm(rstream, M);
sample = rp(1:N);
asample = mod((sample-1), nrows)+1;
apm_sm = apm(asample);
y_sm = y(sample);
ws = warning;            
warning('off', 'stats:statrobustfit:IterationLimit')
fitval = robustfit(apm_sm, y_sm);
warning(ws);



