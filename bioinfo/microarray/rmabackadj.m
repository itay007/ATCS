function y = rmabackadj(x, varargin)
%RMABACKADJ performs the background adjustment of Affymetrix microarray
% probe-level data using the Robust Multi-array Average (RMA) procedure. 
%
%   The perfect match (PM) probe intensities are modeled as the sum of a normal
%   noise component N (normal distribution with mean mu and variance sigma^2)
%   and an exponential signal component S (exponential with mean alpha). To
%   avoid any possibility of negatives, the normal to the-left-of-the-mode is
%   truncated to zero. Given X as the observed intensity, the adjustment to the
%   signal is
% 
%       E(y|X=x) = a + b *(phi(a/b) - phi((x-a)/b))/(Phi(a/b)-Phi((x-a)/b)-1)
% 
%   where a = x - mu - sigma^2 * alpha and b = sigma. phi and Phi are the
%   standard normal distribution pdf and cdf respectively.
% 
%   BACKADJDATA = RMABACKADJ(DATA) returns the background adjusted PM
%   intensities of a chip. DATA can be a matrix with several chips, each column
%   corresponds with each chip, and all have the same probe indices. 
%
%   BACKADJDATA = RMABACKADJ(...,'METHOD',METHOD) sets the estimation method for
%   the background adjustment model parameters. METHOD can be 'RMA' (default) or
%   'MLE'. The details of the 'RMA' estimation method are described by Bolstad
%   [2]. 'MLE' method estimates the parameters using maximum likelihood.
% 
%   BACKADJDATA = RMABACKADJ(...,'TRUNCATE',TF) uses nontruncated Gaussian as
%   background noise model when TF is set to false. By default truncated
%   Gaussian is used.
% 
%   RMABACKADJ(...,'SHOWPLOT',SP) plots the histogram of the log2(PM) values,
%   and the convoluted probability distribution function with estimated mu,
%   sigma and alpha for the chip with column index SP. SP may also be a vector
%   contains the column indices in DATA. Use 'all' for SP to show plots for all
%   the chips. 
%       
%   Example:
%       load prostatecancerrawdata
% 
%       % Background adjust data from one chip
%       backadjData = rmabackadj(pmMatrix(:,1));
% 
%       % Background adjust data from multiple chips, and plot the histogram of
%       % log2(PM) and the estimated background density of chip 1
%       backadjData = rmabackadj(pmMatrix, 'showplot', 1); 
% 
%       % Background adjust data from multiple chips using the MLE estimation
%       % method
%       backadjData = rmabackadj(pmMatrix, 'method', 'mle', 'showplot', 1);
%
%   See also AFFYGCRMA, AFFYINVARSETNORM, AFFYPREPROCESSDEMO, AFFYREAD,
%   AFFYRMA, CELINTENSITYREAD, GCRMA, GCRMABACKADJ, PROBESETLOOKUP,
%   PROBESETVALUES, QUANTILENORM, RMASUMMARY.

% References: 
% [1] Irizarry RA, Hobbs B, Collin F, Beazer-Barclay YD, Antonellis KJ, Scherf
%     U, Speed TP. "Exploration, Normalization, and Summaries of High Density
%     Oligonucleotide Array Probe Level Data" Biostatistics 4, pp249-264,
%     2003.
% [2] B. Bolstad. "affy: Built-in Processing Methods",
%     http://www.bioconductor.org/packages/2.1/bioc/vignettes/affy/inst/doc/builtinMethods.pdf
% [3] Best, C.J.M., Gillespie, J.W., Yi, Y., Chandramouli, G.V.R.,
%     Perlmutter, M.A., Gathright, Y., Erickson, H.S., Georgevich, L., Tangrea,
%     M.A., Duray, P.H., Gonzalez, S., Velasco, A., Linehan, W.M., Matusik,
%     R.J., Price, D.K., Figg, W.D., Emmert-Buck, M.R., and Chuaqui, R.F.
%     Molecular alterations in primary prostate cancer after androgen
%     ablation therapy. Clinical Cancer Research 11, 6823?834, 2005.

% Copyright 2006-2010 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);
% Set defaults
rmaFlag = true; % use strictly RMA bg-adjustment implementation
truncFlag = true; % use to use truncated or nontruncated Gaussian.
plotId = 0;

% Validate input data

if ~isnumeric(x) || ~isreal(x)
   error(message('bioinfo:rmabackadj:ProbeIntensityNotNumericAndReal')) 
end

% Columnzie the input if is a vector
if isvector(x)
  x = x(:);
end

% Get number of probes per column, and number of columns
[N, nChips] = size(x);
colId = 1:nChips;

% deal with the various inputs
if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:rmabackadj:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'method', 'truncate', 'showplot'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:rmabackadj:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:rmabackadj:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % method
                    if ischar(pval)
                        okmethods = {'rma','mle'};
                        nm = strmatch(lower(pval), okmethods);
                        if isempty(nm)
                            error(message('bioinfo:rmabackadj:MethodNameNotValid'));
                        elseif length(nm)>1
                            error(message('bioinfo:rmabackadj:AmbiguousMethodName', pval));
                        else
                            rmaFlag = nm == 1;
                        end
                    else
                        error(message('bioinfo:rmabackadj:MethodNameNotValid'));
                    end

                case 2 % truncate flag
                    truncFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 3 % show plot
                    if isnumeric(pval)
                        if isscalar(pval)
                            pval = double(pval);
                        end
                        
                        tobeplot = ismember(pval, colId);
                        if ~all(tobeplot)
                            warning(message('bioinfo:rmabackadj:InvalidPlotIndices', num2str( plotId( ~tobeplot ) )));
                        end
                        plotId = pval;
                        
                    elseif ischar(pval) && strcmpi(pval, 'all')
                        plotId = colId;
                    else
                        warning(message('bioinfo:rmabackadj:InvalidSP'))
                    end
            end
        end
    end
end

y = x;

% Iterate background adjustment for every sample
for ns = 1:nChips
    if nargout > 0 || (ns == plotId)
        o = x(:,ns);

        % Too few data points, the model will not work properly
        if N < 1000
            o = randsample(o,20000,true);
        end
        
        % estimate parameters
        if rmaFlag
            % estimate initial guess for mu
            mu = findMaxDensity(o);           
            % estimate mu from left-of-the-mode data
            mu = findMaxDensity( o(o < mu));
            
            sigma = estimateSigma(o, mu, rmaFlag);
            alpha = estimateAlpha(o, mu, rmaFlag);
        else %MLE
            % estimate initial guess for mu
            mu = findMaxDensity(o);
            
             % estimate alpha from the right tail
            os = sort(o,'descend');
            midpoint = floor(numel(os)/2);
            alpha = 4/(mean(os(1:midpoint))-os(midpoint));
            
            % estimate sigma
            sigma = 1.2*estimateSigma(o, mu, rmaFlag);
            
            % safe shift to avoid the effects of the truncated normal density
            shift = mean(o);
            
            % Define the join pdf function
            pdf = @(x,mu) alpha*exp(-(x-mu)*alpha+(alpha^2)*(sigma^2)/2) .*...
                (normcdf((x-mu-alpha*sigma^2)/sigma)+normcdf((mu+alpha*sigma^2)/sigma)-1) + eps;

            % Given a formally defined PDF, we can use Maximum likehood estimation
            % function to get better estimates of the mean.
            ws = warning;
            warning('off','MATLAB:optimset:InvalidParamName');
            try
                mu = mle(randsample(o,400)+shift,'pdf',pdf,'start',mu+shift)-shift;
            catch allExceptions %#ok<NASGU>
                warning(message('bioinfo:rmabackadj:MLEErrored'));
            end
            warning(ws);
        end

        a = o - mu - alpha .* (sigma.^2);

        if ~truncFlag || rmaFlag
            % RMA implementation uses untruncated Gaussian
            y(:,ns) = a + sigma .* normpdf(a./sigma) ./ normcdf(a./sigma);
        elseif truncFlag
            y(:,ns) = a + sigma .* ( normpdf(a./sigma) - normpdf((o-a)./sigma) )./...
                ( normcdf(a./sigma) + normcdf((a-o)./sigma) );
        end

        % Show the estimated distributions of the signal and noise
        if any(ns == plotId)
            % Define the join pdf function
            pdf = @(x,mu,sigma,alpha)...
                alpha*exp(-(x-mu)*alpha+(alpha^2)*(sigma^2)/2) .*...
                (normcdf((x-mu-alpha*sigma^2)/sigma)+normcdf((mu+alpha*sigma^2)/sigma)-1);
            figure
            % Plot histogram of log2 PM density
            np = length(o);
            l_o = log2(o);
            bw = .2;
            edges = 0:bw:20;
            hgt = histc(l_o, edges) ./ (np*bw);
            bar(edges,hgt,'histc');

            hh = findobj(gca,'Type','patch');
            set(hh,'FaceColor', [0.5 0.5 0.8], 'EdgeColor', [0.1 0.1 0.5])
            hold on;
            
            % plot estimated background density
            ord = floor(log10(np)); % order ofthe number of points
            P = ceil(np/10^ord)*10^ord;

            Y = pdf(1:P,mu,sigma,alpha);
            yscale = max(hgt)/max(Y);
            hp = plot(log2(1:P),yscale*Y,'r','linewidth',2);
            
            title(sprintf('Histogram of log2(PM) and estimated BG density of sample: #%d', ns));
            xlabel('log2(PM Intensities)')
            ylabel('Frequency')
            hl = legend(hp, 'Estimated background density');
            set(hl, 'box', 'off', 'location','Northeast','Interpreter','none');

            xlim([min(l_o)-2, max(l_o)]);
            hold off;
        end % if ns
    end % if nargout
end % for ns
    
% Using kernel density function to find the intensity position of the mode
% (maximum density. Irizarry et al. assumed the normal distribution of MM
% density. The mode of the histogram is a natural estimate of the mean
% background level. This is mu. But this maybe replace just by mode
% function.
function varargout = findMaxDensity(z)
npoints = min(numel(z)/100, 16384);
[f, x] = ksdensity(z,  min(z):(max(z)-min(z))/npoints:max(z), 'kernel', 'epanechnikov');
[~,h] = max(f); 

varargout{1} = x(h);

% Find maximum density of the sample in z and set the number of
% equally-spaced points in z to 2^14 points.
function maxX = findMaxDensity2(z)
[f, x] = ksdensity(z, 'npoints', 16384, 'kernel', 'epanechnikov');
[~,h] = max(f); 
maxX = x(h);


% Estimate the normal distribution sigma parameter from the
% left-of-the-mode data
function sigma = estimateSigma(z, mu, rmaFlag)
leftz = z(z < mu);

if rmaFlag
    n = numel(leftz);
    sigma = sqrt(sum((leftz - mu).^2) / (n - 1)) * sqrt(2.0); 
else
    sigma = std(leftz);
end

% Estimate the exponantial distribution alpha parameter from the
% right-of-the-mode data
function alpha = estimateAlpha(z, mu, rmaFlag)
% subtract mu from z, use only z-mu > 0
if rmaFlag
    p = z - mu;
    p(p<0)=[];
    alpha = findMaxDensity2(p);
    alpha = 1/alpha;
else 
    p = z - mean(z);
    p(p<0)=[];
    % ML estimator of the decay exponential, observe that exponetial is
    % memoryless, so we can subtract any quantiy and still have a good
    % estimate of the decay, we go to the right enough to avoid measuring
    % the normal noise.
    alpha = 1/mean(p);
end





















