function [X,Y] = msppresample(LP,N,varargin)
%MSPPRESAMPLE peak preserving resample
%
%   [X,Y] = MSPPRESAMPLE(P,N) resamples a list of peaks such that the
%   resulting signal preserves the peak information. P is a matrix with two
%   columns with peak information, the first column has the location of the
%   peaks in the separation axis and the second column has the intensity of
%   the peaks. N is the number of equally spaced points in the resampled 
%   signal. X and Y are column vectors where paired values represent points
%   in the reconstructed signal.
%
%   Units in the X scale (separation units or s.u.) may quantify
%   wavelength, frequency, distance, time or m/z depending on the type of 
%   instrument that generates the signal. 
%
%   P can also be a cell array with multiple peak lists for different
%   scans, in such case, X is a common scale vector and Y is a matrix with
%   the reconstructed signals arranged columnwise. 
%
%   MSPPRESAMPLE reverses the centroid operation performed by bioanalytics
%   instruments in order to have semi-continuous raw signals. When applied
%   to a group of signals MSPPRESAMPLE uniformizes the X vector to a common
%   grid. MSPPRESAMPLE uses a Gaussian kernel to reconstruct the signal,
%   the intensity at any given value in X is given by taking the maximum
%   contribution of all the close peaks.  
%
%   MSPPRESAMPLE is useful for: 1) image representation of a set of
%   signals, 2) perform comparative analysis between different sets of
%   signals, or, 3) extract region profiles in multidemensional datastes,
%   for example the chromatographic time profile in LCMS or GCMS datasets. 
%
%   MSPPRESAMPLE(...,'RANGE',R) uses a 1-by-2 vector with the range for the
%   desired output signal. R(1) defaults to the earliest peak found in P
%   and R(2) defaults to latest peak found in P. 
%
%   MSPPRESAMPLE(...,'FWHH',F) sets the full width at half height (FWHH) in
%   s.u. The FWHH is used to convert each peak to a Gaussian shaped curve.
%   F is a scalar value and defaults to median(diff(X))/2. To ensure that
%   the resolution of the peaks is preserved, set FWHH value to half the
%   distance between the two peaks of interest that are closest to each
%   other. 
%
%   MSPPRESAMPLE(...,'SHOWPLOT',SP) plots the original peaks and the
%   resampled signal. When SP is TRUE, the first list of peaks in P is
%   used. If MSPPRESAMPLE is called without output arguments, a plot will
%   be shown unless SP is false. SP may also contain an index to one of the
%   list of peaks in P.
%
%   Examples:
%
%     load lcmsdata
%
%     % Resample to create a heat map of the whole LCMS sample:
%     [MZ,Y] = msppresample(ms_peaks,5000);
%     msheatmap(MZ,ret_time,log(Y))
%
%     % Plot the reconstructed profile spectra between two times:
%     figure
%     t1 = 3370;
%     t2 = 3390;
%     h = find(ret_time>t1 & ret_time<t2);
%     [MZ,Y] = msppresample(ms_peaks(h),10000);
%     plot3(repmat(MZ,1,numel(h)),repmat(ret_time(h)',10000,1),Y)
%     xlabel('Mass/Charge (M/Z)')
%     ylabel('Retention Time')
%     zlabel('Relative Intensity')
%
%     % Resample to plot the Total Ion Chromatogram (TIC):
%     figure
%     [MZ,Y] = msppresample(ms_peaks,5000);
%     plot(ret_time,sum(Y))
%     title('Total Ion Chromatogram (TIC)')
%     xlabel('Retention Time')
%     ylabel('Relative Intensity')
%     
%     % Resample to plot the eXtracted Ion Chromatogram (XIC) in the
%     % 450-500 M/Z range: 
%     figure
%     [MZ,Y] = msppresample(ms_peaks,5000,'Range',[450 500]);
%     plot(ret_time,sum(Y))
%     title('Extracted Ion Chromatogram (XIC) from 450 to 500 M/Z')
%     xlabel('Retention Time')
%     ylabel('Relative Intensity')
%
%   See also DIFFPROTDEMO, LCMSDEMO, MSPEAKS, MSPALIGN, MSRESAMPLE,
%   MZXML2PEAKS, MZXMLREAD.

%   Copyright 2006-2008 The MathWorks, Inc.


% References:
%  [1] Lars Linsen, Julia Locherbach, Matthias Berth, Dorte Becher, Jorg 
%      Bernhardt, "Visual Analysis of Gel-Free Proteome Data," IEEE
%      Transactions on Visualization and Computer Graphics ,vol. 12, no. 4,
%      pp. 497-508, July/August, 2006.   

% check inputs
bioinfochecknargin(nargin,2,mfilename);

if ~iscell(LP) 
    LP = {LP}; 
end

if nargout > 0 
    plotId = 0; 
else
    plotId = 1;
end

% setting some defaults
KT = 0.05; % default kernel threshold, kernel values below this threshold 
           % (relative to the peak height) are not considered  
FWHH = []; % default peak resolution 
Xmin = max(0,min(cellfun(@(x) min(x(:,1)),LP))); % default min(X)
Xmax = max(cellfun(@(x) max(x(:,1)),LP)); % default max(X)
M = numel(LP);

if  nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:msppresample:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'range','fwhh','showplot'};
    for j=1:2:nargin-3
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:msppresample:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:msppresample:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % range
                    if numel(pval)~=2 || diff(pval)<=0 
                        error(message('bioinfo:msppresample:badRange'))
                    end
                    Xmin = max(pval(1),Xmin);
                    Xmax = min(pval(2),Xmax);
                case 2 % fwhh
                    if ~isnumeric(pval) || ~isscalar(pval) || pval<=0 
                        error(message('bioinfo:msppresample:InvalidFWHH'))
                    end
                    FWHH = pval;
                case 3 % show
                    if bioinfoprivate.opttf(pval) 
                        if isnumeric(pval)
                            if isscalar(pval)
                                plotId = double(pval); 
                            else
                                plotId = 1;
                                warning(message('bioinfo:msppresample:SPNoScalar'))
                            end
                        else
                            plotId = 1;
                        end
                    else
                        plotId = 0;
                    end        
            end
        end
    end
end

OPD = min(cellfun(@(x) median(diff(x(x(:,1)>Xmin & x(:,1)<Xmax,1))),LP)/2); % observed peak distance
if isempty(FWHH)
    FWHH = OPD; % default peak resolution 
end

% calculating some constants
Xstep = (Xmax-Xmin)/(N-1); 
X = (Xmin:Xstep:Xmax)';
SIGMA = -log(.5)/(FWHH/2)^2;
SPHG = ceil(sqrt(log(KT)/-SIGMA)/Xstep);
SPG = 2*SPHG+1;

if FWHH > OPD
   warning(message('bioinfo:msppresample:LowResolution')) 
end
if Xstep * 2 > OPD
    warning(message('bioinfo:msppresample:TooFewPoints'))
end

% do the resampling
Y = zeros(N,M,class(LP{1}));
for i = 1:M
if nargout>0 || (i == plotId)    
    xx = LP{i}(:,1);
    t = repmat(round((xx-Xmin)/Xstep),1,SPG)+repmat(-SPHG:SPHG,numel(xx),1);
    y = exp(-SIGMA*(Xmin+t*Xstep - repmat(xx,1,SPG)).^2).*repmat(LP{i}(:,2),1,SPG);
    Y(:,i) = accumarray(min(numel(X),max(1,t(:)+1)),y(:),[numel(X),1],@max);
    if (i == plotId)
       figure
       plot(X,Y(:,i));
       hold on
       plot(LP{i}(:,1),LP{i}(:,2),'xr','linewidth',2);
       title(sprintf('Signal ID: %d',i));
       xlabel('Separation Units')
       ylabel('Relative Intensity')
       legend('Resampled signal','Peaks')
       axis([min(X) max(X) min(Y(:,i)) max(Y(:,i))])
       setAllowAxesRotate(rotate3d(gcf),gca,false)
       grid on
       hold off       
    end
end    
end

if nargout == 0 
    clear X
end




