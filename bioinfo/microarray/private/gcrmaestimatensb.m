function [sigma,mu,mu2,mu3] = gcrmaestimatensb(Y,affi,affi2, varargin)
%GCRMAESTIMATENSB estimates non-specific binding background parameters.
% 
%   [SIGMA, MU, MU2] = GCRMAESTIMATENSB(INTENSITIES, ALPHA1, ALPHA2)returns
%   non-specific binding (NSB) parameters SIGMA and MU, estimated from MM
%   INTENSITIES and probe affinities ALPHA1 with known sequences.
%   INTENSITIES are the MM intensity values after optical background
%   correction, but before NSB correction. Affiniteis ALPHA2 are for the PM
%   probes whose expected NSB intensity is to be predicted from the
%   returned NSB parameters MU2.
% 
%   [SIGMA, MU, MU2, MU3] = GCRMAESTIMATENSB(..., 'EXTRAGROUP', ALPHA3)
%   returns NSB parameters MU3 for extra group of probes with affinites
%   ALPHA3 whose expected NSB intensity is to be predicted from the
%   returned parameters.
% 
%   [...] = GCRMAESTIMATENSB(...,'SHOWPLOT',TF) displays log2(INTENSITIES)
%   vs. AFFINITIES for NSB data if TF is true.
% 
%   Example:
%       load affy95a_spikein
%       [sigma, mu_mm, mu_pm] = gcrmaestimatensb(mm_oc(1,:), affi_mm, affi_pm,...
%                     'showplot', true);
% 
%   See also AFFYPROBESEQREAD, AFFYREAD, CELLINTENSITYREAD, PROBELIBRARYINFO.

%   References:
%       [1] W Zhijin Wu, Rafael A. Irizarry, Robert Gentleman, Francisco
%           Martinez Murillo, and Forrest Spencer, "A Model Based
%           Background Adjustment for Oligonucleotide Expression Arrays", J
%           Amer Stat Assoc Vol.99, No.468, 2004, pp.909-917.
%       [2] Zhijin Wu and Rafael A. Irizarry, "A Statistical Framework for
%           the Analysis of Microarray Probe-Level Data", Johns Hopkins
%           University, Biostatistics Working Papers 73 (2005).
%       [3] T.Speed, "Background models and GCRMA", Lecture 10, Statistics
%           246, UC Berkeley, Spring 2006
%           http://www.stat.berkeley.edu/users/terry/Classes/s246.2006/Week
%           10/Week10L1.pdf.

%   Copyright 2006 The MathWorks, Inc.


if nargin < 3
    error(message('bioinfo:gcrmaestimatensb:NotEnoughInputs', mfilename));
end

if ~isnumeric(Y) || ~isreal(Y) || ~isvector(Y)
   error(message('bioinfo:gcrmaestimatensb:ProbeIntensityNotNumericAndRealVector')) 
end

if ~isnumeric(affi) || ~isreal(affi) || ~isvector(affi)
   error(message('bioinfo:gcrmaestimatensb:AffinityNotNumericAndRealVector')) 
end

if numel(Y) ~= numel(affi)
    error(message('bioinfo:gcrmaestimatensb:IntensityandAffinityLengthNotMatch')); 
end

if ~isempty(affi2) && (~isnumeric(affi2) || ~isreal(affi2) ||~isvector(affi2))
    error(message('bioinfo:gcrmaestimatensb:AffinityNotNumericAndRealVector'))
else
    affi2 = affi2(:);
end

Y = Y(:);
affi = affi(:);

affi3 = [];
showplot = false;
mu2 = [];
mu3 = [];

% Deal with the various input
if nargin > 3
    if rem(nargin, 2) == 0
        error(message('bioinfo:gcrmaestimatensb:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'extragroup', 'showplot'};
    for j = 1:2:nargin-3
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs, numel(pname)));

        if isempty(k)
            error(message('bioinfo:gcrmaestimatensb:UnknownParameterName', pname));
        elseif length(k) > 1
            error(message('bioinfo:gcrmaestimatensb:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % extra group
                    affi3 = pval;
                    if ~isnumeric(affi3) || ~isreal(affi3) || ~isvector(affi3)
                        error(message('bioinfo:gcrmaestimatensb:ExtraAffinityNotNumericAndRealVector'))
                    end
                    affi3 = affi3(:);
                case 2 % showplot
                   showplot = bioinfoprivate.opttf(pval);
                    if isempty(showplot)
                        error(message('bioinfo:gcrmaestimatensb:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
            end % switch
        end % if
    end % for loop
end % if

nanIndices = isnan(affi);
y = log(Y(~nanIndices));
affi_good = affi(~nanIndices);

% The number of data point in an array over 100,000 is too big for
% malowess, try to reduce the number by random permutation of the probes,
% and pick the first 5000 points to do the smoothing.

M = numel(y);
N =  min(M, 5000);

sample = randsample(M, N);
aff_sm = affi_good(sample);
y_sm = y(sample);

% Smoothing with lowess
y_sm_smooth = malowess(aff_sm, y_sm, 'span', 0.175); 

% Find mu for affinity. Need to find unique aff_sm, otherwise interp1 will
% error.
[aff_sm_u, idx_u] = unique(aff_sm);
y_sm_smooth_u = y_sm_smooth(idx_u);

% Interpolate to get mu
mu = interp1(aff_sm_u, y_sm_smooth_u, affi, 'linear', 'extrap');

if ~isempty(affi2)
    mu2 = interp1(aff_sm_u, y_sm_smooth_u, affi2, 'linear', 'extrap');
end

if ~isempty(affi3)
    mu3 = interp1(aff_sm_u, y_sm_smooth_u, affi3, 'linear', 'extrap');
end

% Compute sigma estimate using MAD of negative residuals.
res = y_sm - y_sm_smooth;
sigma = madsigma(res);

if showplot
    plotsmoothfit(affi, log2(Y), mu)
end

%-----------------------------------------------------------------------%
function plotsmoothfit(x,y, mu)
% Plot the optical noise adjusted log2(MM) vs affinity for NSB data. The
% solid line is a lowess fit.

% sort x for display a line
[x_sort, id_sort] = sort(x); 
y_interp = log2(exp(mu(id_sort)));
% y_interp = mu(id_sort);
figure
plot(x, y, '.');

hold on
l1 = plot(x_sort, y_interp, 'r', 'linewidth', 2);
legend(l1, 'LOWESS smooth fit');
xlabel('Affinity');
ylabel('Log2(Optically adjusted MM probe intensities)');
hold off

function sigma = madsigma(r)
% Compute sigma estimate using MAD of negative residuals.
nr = r(r < 0);
r = [nr; -nr];
sigma = mad(r, 1)/ norminv(0.75);
