function Y = msbackadj(X,Y,varargin)
%MSBACKADJ provides background correction for a signal with peaks
%
%   YOUT = MSBACKADJ(X,Y) adjusts the variable background (baseline) of a
%   signal with peaks by following three steps: 1) estimates the background
%   within multiple shifted windows of width 200 separation units along the
%   x-axis, 2) regresses the varying baseline to the window points using a
%   spline approximation, and 3) adjusts the background of the input signal Y.  
%
%   X and Y are column vectors where paired values represent points in the
%   signal. Y can be a matrix with several signals, all sharing the same X
%   scale. Units in the X scale (separation units or s.u.) may quantify
%   wavelength, frequency, distance, time or m/z depending on the type of
%   instrument that generates the signal. 
%
%   MSBACKADJ(...,'WINDOWSIZE',WS) sets the width for the shifting window.
%   The default is 200 s.u., which means a background point is estimated
%   for windows of 200 s.u. in width. WS may also be a function handle.
%   The referenced function is evaluated at each X value to compute 
%   a variable width for the windows. This option is useful for cases where
%   the resolution of the signal is dissimilar at different regions.
%
%   MSBACKADJ(...,'STEPSIZE',SS) sets the steps for the shifting window.
%   The default is 200 s.u., which means a background point is estimated
%   for windows at every 200 s.u. SS may also be a function handle. The
%   referenced function is evaluated at each X value to compute the
%   distance between adjacent windows.
%
%   MSBACKADJ(...,'REGRESSIONMETHOD',RM) sets the method used to regress
%   the window estimated points to a soft curve. The default is 'pchip';
%   i.e., shape-preserving piecewise cubic interpolation. Other options are
%   'linear' and 'spline' interpolation.
%
%   MSBACKADJ(...,'ESTIMATIONMETHOD',EM) sets the method used to find the
%   likely background value at every window. Default is 'quantile', in
%   which the quantile value is set to 10%. An alternative method is 'em',
%   which assumes a doubly stochastic model; i.e., every sample is the
%   i.i.d. draw of any of two normal distributed classes (background or
%   peaks). Because the class label is hidden, the distributions are
%   estimated with an expectation-maximization algorithm. The ultimate
%   background value is the mean of the background class.
%
%   MSBACKADJ(...,'SMOOTHMETHOD',SM) sets the method used to smooth the
%   curve of estimated points, useful to eliminate the effect of possible
%   outliers. Options are 'none' (default), 'lowess' (linear fit), 'loess'
%   (quadratic fit), or 'rlowess' and 'rloess' (robust linear and quadratic
%   fit).
%
%   MSBACKADJ(...,'QUANTILEVALUE',QV) changes the default quantile value.
%   The default is 0.10.
%
%   MSBACKADJ(...,'PRESERVEHEIGHTS',true) sets the baseline subtraction
%   mode to preserve the height of the tallest peak in the signal when
%   subtracting the baseline. By default heights are not preserved.
%
%   MSBACKADJ(...,'SHOWPLOT',SP) plots the background estimated points, the
%   regressed baseline, and the original signal. When SP is TRUE, the first
%   signal in Y is used. If MSBACKADJ is called without output arguments, a
%   plot will be shown unless SP is FALSE. SP can also contain an index to
%   one of the signals in Y.
%
%   Examples: 
%
%      % Correct the baseline of SELDI-TOF mass-spectrograms:
%      load sample_lo_res
%
%      % Adjust the baseline of a group of spectrograms.
%      YB = msbackadj(MZ_lo_res,Y_lo_res);
%
%      % Plot the third spectrogram in Y_lo_res and its estimated baseline.
%      msbackadj(MZ_lo_res,Y_lo_res,'SHOWPLOT',3);
%
%      % Plot the estimated baseline for the fourth spectrogram in Y_lo_res
%      % using an anonymous function to describe a MZ dependent parameter.
%      wf = @(mz) 200 + .001 .* mz;
%      msbackadj(MZ_lo_res,Y_lo_res(:,4),'STEPSIZE',wf);
%
%
%   See also MSALIGN, MSHEATMAP, MSLOWESS, MSNORM, MSPREPRODEMO,
%   MSRESAMPLE, MSSGOLAY, MSVIEWER. 

%   Copyright 2003-2008 The MathWorks, Inc.


% References: 
% [1] Lucio Andrade and Elias Manolakos, "Signal Background Estimation and
%     Baseline Correction Algorithms for Accurate DNA Sequencing" Journal
%     of VLSI, special issue on Bioinformatics 35:3 pp 229-243 (2003)

% check inputs
bioinfochecknargin(nargin,2,mfilename);
% set defaults
stepSize = 200;
windowSize = 200;
regressionMethod = 'pchip';
estimationMethod = 'quantile';
smoothMethod = 'none';
quantileValue = 0.1;
preserveHeights = false;
maxNumWindows = 1000;
if nargout == 1
    plotId = 0; 
else
    plotId = 1;
end

if  nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:msbackadj:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'stepsize','windowsize','regressionmethod',...
              'estimationmethod','quantilevalue','preserveheights',...
              'smoothmethod','showplot'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:msbackadj:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:msbackadj:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % step size
                    stepSize = pval;
                    if ~isscalar(stepSize) && ~isa(stepSize,'function_handle')
                        error(message('bioinfo:msbackadj:NotValidStepSize'))
                    end
                case 2 % window size
                    windowSize = pval;
                    if ~isscalar(windowSize) && class(windowSize,'function_handler')
                        error(message('bioinfo:msbackadj:NotValidWindowSize'))
                    end
                case 3 % regression method
                    regressionMethods = {'cubic','pchip','spline','linear'};
                    regressionMethod = strmatch(lower(pval),regressionMethods); 
                    if isempty(regressionMethod) 
                        error(message('bioinfo:msbackadj:NotValidRegressionMethod'))
                    end
                    regressionMethod = regressionMethods{max(2,regressionMethod)};
                case 4 % estimation method
                    estimationMethods = {'quantile','em'};
                    estimationMethod = strmatch(lower(pval),estimationMethods); 
                    if isempty(estimationMethods) 
                        error(message('bioinfo:msbackadj:NotValidEstimationMethod'))
                    end
                    estimationMethod = estimationMethods{estimationMethod};
                case 5 % quantile value
                    quantileValue = pval;
                case 6 % preserve heights
                    preserveHeights = bioinfoprivate.opttf(pval);
                case 7 % smoothing method
                    smoothMethods = {'none','lowess','loess','rlowess','rloess'};
                    smoothMethod = strmatch(lower(pval),smoothMethods); 
                    if isempty(smoothMethod) 
                        error(message('bioinfo:msbackadj:NotValidSmoothMethod'))
                    elseif length(smoothMethod)>1
                        error(message('bioinfo:msbackadj:AmbiguousSmoothMethod', pname));
                    end
                    smoothMethod = smoothMethods{smoothMethod};
                case 8 % show
                    if bioinfoprivate.opttf(pval) 
                        if isnumeric(pval)
                            if isscalar(pval)
                                plotId = double(pval); 
                            else
                                plotId = 1;
                                warning(message('bioinfo:msbackadj:SPNoScalar'))
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

% validate X and Y

if ~isnumeric(Y) || ~isreal(Y)
   error(message('bioinfo:msbackadj:IntensityNotNumericAndReal')) 
end

if ~isnumeric(X) || ~isreal(X)
   error(message('bioinfo:msbackadj:XNotNumericAndReal')) 
end

if size(X,1) ~= size(Y,1)
   error(message('bioinfo:msbackadj:NotEqualNumberOfSamples'))
end
 
numSignals = size(Y,2);

if (plotId~=0) && ~any((1:numSignals)==plotId)
    warning(message('bioinfo:msbackadj:InvalidPlotIndex'))
end

multiple_X = false;
if size(X,2)>1
    multiple_X = true;
    if size(X,2) ~= numSignals
        error(message('bioinfo:msbackadj:NotEqualNumberOfXScales'))
    end
end

% change scalars to function handlers
if isnumeric(stepSize)   
    stepSize   = @(x) repmat(stepSize,size(x));   
end
if isnumeric(windowSize) 
    windowSize = @(x) repmat(windowSize,size(x)); 
end

% allocate space for Xp and WE
Xp = zeros(maxNumWindows,1);
WE = nan(maxNumWindows,1);

% calculates the location of the windows (when it is the same for all the
% signals)
if ~multiple_X
    Xpid = max(0,X(1));
    Xend = X(end);
    id = 1;
    while Xpid <= Xend
        Xp(id) = Xpid;
        Xpid = Xpid + stepSize(Xpid);
        id = id + 1;
        if id > maxNumWindows
            error(message('bioinfo:msbackadj:maxNumWindowsExceeded'))
        end
    end
    numWindows = id-1;
end

% iterate for every signal
for ns = 1:numSignals 
if nargout>0 || (ns == plotId)
    % find the location of the windows (when it is different for every
    % signal, otherwise this was done out of the loop)
    if multiple_X
       Xpid = max(0,X(1,ns));
       Xend = X(end,ns);
       id = 1; Xp = zeros(1,maxNumWindows);
       while Xpid <= Xend
           Xp(id) = Xpid;
           Xpid = Xpid + stepSize(Xpid);
           id = id + 1;
           if id > maxNumWindows
                error(message('bioinfo:msbackadj:maxNumWindowsExceeded'))
           end
       end
       Xp(id:end)=[];
       numWindows = id-1;
       nnss = ns;
    else
       nnss = 1;
    end    
    Xpt = Xp(1:numWindows); 
    Xw = windowSize(Xpt);
    
    % find the estimated baseline for every window
    for nw = 1:numWindows
        subw = Y(X(:,nnss)>=Xpt(nw) & X(:,nnss)<= (Xpt(nw)+Xw(nw)),ns);
        switch estimationMethod
            case 'quantile'
                WE(nw) = quantile(subw,quantileValue);
            case 'em'
                WE(nw) = em2c1d(subw);
        end
    end % for nw = 1:numWindows
    
    % smooth the estimated points
    if ~isequal('none',smoothMethod)
        WE(1:numWindows) = ...
            bioinfoprivate.masmooth(Xpt+Xw/2,WE(1:numWindows),10,smoothMethod,2);
    end
            
    % regress the estimated points
    b = interp1(Xpt+Xw/2,WE(1:numWindows),X(:,nnss),regressionMethod);
    
    if (ns == plotId)
       figure
       plot(X(:,nnss),Y(:,ns))
       hold on
       plot(X(:,nnss),b,'r','linewidth',2)
       plot(Xpt+Xw/2,WE(1:numWindows),'kx')
       title(sprintf('Signal ID: %d',ns));
       xlabel('Separation Units')
       ylabel('Relative Intensity')
       legend('Original Signal','Regressed baseline','Estimated baseline points')
       axis([min(X(:,nnss)) max(X(:,nnss)) min(Y(:,ns)) max(Y(:,ns))])
       grid on
       hold off
       setAllowAxesRotate(rotate3d(gcf),gca,false)
    end
    
    % apply the correction
    if preserveHeights
        K = 1 - b/max(Y(:,ns));
        Y(:,ns) = (Y(:,ns) - b) ./ K;
        %[YMax,locMax] = max(Y(:,ns));
        %K = 1 - b(locMax)/YMax;
        %Y(:,ns) = (Y(:,ns) - b) / K;
    else
        Y(:,ns) = (Y(:,ns) - b);
    end 

end % if nargout>0 || (ns == plotId)    
end % for ns = 1:numSignals 

if nargout == 0 
    clear Y
end



    
    
