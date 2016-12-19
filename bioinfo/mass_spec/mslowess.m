function Y = mslowess(X,Y,varargin)
%MSLOWESS provides nonparametric smoothing of signals with peaks
%
%   YOUT = MSLOWESS(X,Y) smoothes the signal Y using the 'Lowess' method
%   with a default span of 10 samples. 
%
%   X and Y are column vectors where paired values represent points in the
%   signal. Y can be a matrix with several signals, all sharing the same X
%   scale. Units in the X scale (separation units or s.u.) may quantify
%   wavelength, frequency, distance, time or m/z depending on the type of
%   instrument that generates the signal. 
%
%   Notes:
%   1) MSLOWESS assumes a X vector that is not uniformly spaced; therefore,
%   the sliding window for smoothing is centered using the closest samples
%   in terms of the X value and not in terms of the X index. 2) When the X
%   vector does not have repeated values or NaNs, the algorithm will be
%   approximately two times faster.
%
%   MSLOWESS(...,'ORDER',O) sets the order O of the 'lowess' smoother.
%   Order can be 1 (linear fit or Lowess) or 2 (quadratic fit or Loess).
%   The default order is 1. Order can also be 0, which is equivalent to a
%   weighted local mean estimator, presumably faster, because only a mean
%   computation is performed instead of a least squares regression. 
%
%   Note: the MATLAB Curve Fitting Toolbox refers to Lowess smoothing of
%   order 2 as Loess smoothing. 
%
%   MSLOWESS(...,'SPAN',S) modifies the window size for the smoothing
%   kernel. If SPAN is greater than 1, then the window is of size span (in
%   samples, i.e., independently of the X vector). The default value is 10
%   samples. Higher values smooth the signal more at the expense of
%   computation time.  If S is less than 1, then the window size is taken
%   to be a fraction of the number of points in the data; e.g., for S =
%   0.005, the window size is equal to 0.50% of the number of points in X.
%
%   MSLOWESS(...,'KERNEL',K) selects the kernel function. The kernel
%   function is used to weight the observed ion intensities such that
%   those samples close to the X location being smoothed have most weight
%   in determining the estimate. Options are
%
%       'tricubic' (default)    (1 - (dist/dmax).^3).^3  
%       'gaussian'              exp(-(2*dist/dmax).^2)
%       'linear'                1-dist/dmax
%
%   MSLOWESS(...,'ROBUSTITERATIONS',I) sets the number of iterations of a
%   robust fit. If I is 0 (default), no robust fit is performed. For robust
%   smoothing, small residual values at every span are weighted to
%   improve the new estimate. 1 or 2 robust iterations are usually
%   adequate; larger values can be computationally expensive.
%
%   Note: For a uniformly spaced X vector a nonrobust smoothing with
%   order 0 is equivalent to filtering the signal with the KERNEL vector.
%
%   MSLOWESS(...,'SHOWPLOT',SP) plots the smoothed signal over the
%   original. When SP is TRUE, the first signal in Y is used. If MSLOWESS
%   is called without output arguments, a plot will be shown unless SP is
%   FALSE. SP can also contain an index to one of the signals in Y.
%
%   Example: 
%
%       load sample_lo_res
%
%       % Smooth a group of spectrograms.
%       YS = mslowess(MZ_lo_res,Y_lo_res);
%
%       % Plot the third spectrogram in Y_lo_res and its smoothed signal.
%       mslowess(MZ_lo_res,Y_lo_res,'SHOWPLOT',3);
%
%   See also MSALIGN, MSBACKADJ, MSBATCHPROCESSING, MSHEATMAP, MSNORM,
%   MSPREPRODEMO, MSRESAMPLE, MSSGOLAY, MSVIEWER. 

% Copyright 2003-2008 The MathWorks, Inc.


% check inputs
bioinfochecknargin(nargin,2,mfilename);
% set default parameters
kernel = 'tricubic';
span = 10;
order = 1; % 'lowess'
robustIter = 0;
if nargout == 1
    plotId = 0; 
else
    plotId = 1;
end


% deal with the various inputs
if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:mslowess:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'kernel','robustiterations','span','order','showplot'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:mslowess:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:mslowess:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % method
                    kernels = {'tricubic','gaussian','linear'};
                    kernel = strmatch(lower(pval),kernels); 
                    if isempty(kernel)
                        error(message('bioinfo:mslowess:NotValidKernel'))
                    end
                    kernel=kernels{kernel};
                case 2 % robust iterations
                    if ~isnumeric(pval) || ~isscalar(pval)
                        error(message('bioinfo:mslowess:LowessSpanNumeric'));
                    end
                    robustIter = pval;
                    if robustIter < 0
                        error(message('bioinfo:mslowess:IterationsPositive'));
                    end
                case 3  % span
                    if ~isnumeric(pval) || ~isscalar(pval)
                        error(message('bioinfo:mslowess:LowessSpanNumeric'));
                    end
                    span = pval;
                    if span < 0
                        error(message('bioinfo:mslowess:LowessSpanPositive'));
                    end
                    if span >numel(X)
                        error(message('bioinfo:mslowess:LowessSpanLessThanX'));
                    end
                case 4 % order
                    if ~isnumeric(pval) || ~isscalar(pval)
                        error(message('bioinfo:mslowess:OrderNumeric'));
                    end
                    order = pval;
                    if ~ismember(order,[0 1 2])
                        error(message('bioinfo:mslowess:LowessBadNumericOrder'));
                    end
                 case 5 % show
                    if bioinfoprivate.opttf(pval) 
                        if isnumeric(pval)
                            if isscalar(pval)
                                plotId = double(pval); 
                            else
                                plotId = 1;
                                warning(message('bioinfo:mslowess:SPNoScalar'))
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

% set the method to call masmooth
methods = {'mean','lowess','loess'};
method = methods{order + 1};
if robustIter > 0
    method = ['r',method];
end

% validate X and Y
if size(X,2) ~= 1
   error(message('bioinfo:mslowess:OnlyOneXScale'))
end
if ~isnumeric(Y) || ~isreal(Y)
   error(message('bioinfo:mslowess:IntensityNotNumericAndReal')) 
end
if ~isnumeric(X) || ~isreal(X)
   error(message('bioinfo:mslowess:XNotNumericAndReal')) 
end
if size(X,1) ~= size(Y,1)
   error(message('bioinfo:mslowess:NotEqualNumberOfSamples'))
end

numSignals = size(Y,2);

if (plotId~=0) && ~any((1:numSignals)==plotId)
    warning(message('bioinfo:mslowess:InvalidPlotIndex'))
end

% iterate for every signal
for ns = 1:numSignals
if nargout>0 || (ns == plotId)
    
    if (ns == plotId)
        yo = Y(:,ns);
    end
    
    Y(:,ns) = bioinfoprivate.masmooth(X,Y(:,ns),span,method,robustIter,kernel);
    
    if (ns == plotId)
        figure
        plot(X,[yo Y(:,ns)])
        title(sprintf('Signal ID: %d',ns));
        xlabel('Separation Units')
        ylabel('Relative Intensity')
        legend('Original signal','Smoothed signal')
        axis([min(X) max(X) min(Y(:,ns)) max(Y(:,ns))])
        grid on
        hold off
        setAllowAxesRotate(rotate3d(gcf),gca,false)
    end
    
end % if nargout>0 || (ns == plotId)    
end % for ns = 1:numSignals 

if nargout == 0 
    clear Y
end
