function [XOUT,YOUT,YF,filtLength] = msresample(X,Y,N,varargin)
%MSRESAMPLE resamples a signal with peaks
%
%   [XOUT,YOUT] = MSRESAMPLE(X,Y,N) resamples a signal with peaks to N
%   equally or linearly spaced points. 
%
%   X and Y are column vectors where paired values represent points in the
%   signal. Y can be a matrix with several signals, all sharing the same X
%   scale. Units in the X scale (separation units or s.u.) may quantify
%   wavelength, frequency, distance, time or m/z depending on the type of
%   instrument that generates the signal. 
%
%   MSRESAMPLE analyzes the X vector and automatically detects if the
%   distances (in X) between samples are constant or increase/decrease
%   linearly, independently of missing samples and accuracy errors. The
%   same detected relation is used to set N points in the output signal
%   (XOUT) between the range [min(X) max(X)].
%
%   When input arguments are set such that down-sampling takes place,
%   MSRESAMPLE first applies a low pass filter to minimize aliasing. For
%   the antialias filter a linear-phase FIR filter is designed using
%   least-squares error minimization. The cutoff frequency is automatically 
%   set by the largest down-sampling ratio when comparing the same regions
%   of the X and XOUT vectors.
%
%   [XOUT,YOUT] = MSRESAMPLE(X,Y,XOUT) resamples a signal with peaks to the
%   points specified in XOUT. Use this option to resample a new signal to
%   the same X scale of a previously sampled signal. When XOUT is given
%   MSRESAMPLE does not analyze the scale in the X vector, MSRESAMPLE only
%   applies the antialias filter and interpolates the output signal YOUT.
%
%   MSRESAMPLE(...,'UNIFORM',TRUE) forces the spaces between samples in the
%   output signal XOUT to be uniform. Defaults to FALSE, MSRESAMPLE
%   automatically selects the spacing between samples.
%
%   MSRESAMPLE(...,'RANGE',R) uses a 1-by-2 vector with the range for the
%   desired output signal. R defaults to [min(X) max(X)]. When R values go
%   beyond the values in X, MSRESAMPLE extrapolates the signal with zeros
%   and returns a warning message. 
%
%   MSRESAMPLE(...,'RANGEWARNOFF',TRUE) omits the RANGE warning message.
%   Defaults to FALSE.  
%
%   MSRESAMPLE(...,'MISSING',TRUE) analyzes the X vector for dropped
%   samples. Defaults to FALSE. Note: checking for dropped samples can be
%   time consuming and might not be worthwhile if the down-sample factor is
%   large. Dropped samples can only be recovered if the original X values
%   follow a linear or a quadratic function of the sample index.
%
%   MSRESAMPLE(...,'WINDOW',W) sets the window used during the filter
%   design process. W defaults to 'Flattop'. Other options are 'Blackman',
%   'Hamming', and 'Hanning'.
%
%   MSRESAMPLE(...,'CUTOFF',F) sets the cutoff frequency. F can be a scalar
%   value between 0 and 1. 1 corresponds to the Nyquist frequency or half
%   the sampling frequency. By default, MSRESAMPLE automatically estimates 
%   F by inspecting X and XOUT; however, F can be underestimated if X
%   presents anomalies. 
%
%   MSRESAMPLE(...,'SHOWPLOT',SP) plots the original and the resampled
%   signals. When SP is TRUE, the first signal in Y is used. If MSRESAMPLE
%   is called without output arguments, a plot will be shown unless SP is
%   false. SP may also contain an index to one of the signals in Y.
%  
%   Example:
%
%     load sample_hi_res
%     MZ = MZ_hi_res;
%     Y = Y_hi_res;
%    
%     % Resample the spectrogram by a 1/2 factor.
%     N = numel(MZ)/2  % find out the number of required samples
%     [mz1,y1] = msresample(MZ,Y,N);
%
%     % Resample the spectrogram to have 10000 samples between 2000 and
%     % 11000 Da. and plot the original and the resampled signals.
%     R = [2000 11000]; % set the range
%     [mz2,y2] = msresample(MZ,Y,10000,'RANGE',R,'SHOWPLOT',true);
%
%   See also MSALIGN, MSBACKADJ, MSHEATMAP, MSLOWESS, MSNORM, MSPREPRODEMO,
%   MSSGOLAY, MSVIEWER. 

%   Copyright 2003-2012 The MathWorks, Inc.


% check inputs
bioinfochecknargin(nargin,3,mfilename);

%set defaults
checkMissing = false;
evenlySpaced = false;
windowType = 'flattop';
cutoffGiven = false;
rangewarn = true;

if nargout == 0
    plotId = 1; 
else
    plotId = 0;
end

% validate required inputs 
if size(X,1)~=size(Y,1)
    error(message('bioinfo:msresample:differentSize'))
end
if ~isnumeric(Y) || ~isreal(Y) 
   error(message('bioinfo:msresample:IntensityNotNumericAndReal')) 
end
if ~isnumeric(X) || ~isreal(X) || ~isvector(X)
   error(message('bioinfo:msresample:XNotNumericAndReal')) 
end

if isscalar(N) && rem(N,1)==0
    xoutgiven = false;
    R = [max(0,min(X)) max(X)];
elseif isnumeric(X) && isreal(X) && isvector(X)
    xoutgiven = true;
    XOUT = N;
    N = numel(XOUT);
    R = [max(0,min(XOUT)) max(XOUT)];
else
    error(message('bioinfo:msresample:badNValue'))
end

if  nargin > 3
    if rem(nargin,2) == 0
        error(message('bioinfo:msresample:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'range','uniform','missing','window','showplot','cutoff','rangewarnoff'};
    for j=1:2:nargin-3
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isequal(k,[1 7])
            k = 1; % special case because 'range' is a partial completion 
                   % of 'rangewarnoff' and we need to be backwards compatible 
                   % to previous cases of incomplete use of 'range'
        end
        if isempty(k)
            error(message('bioinfo:msresample:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:msresample:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % range
                    if xoutgiven 
                        warning(message('bioinfo:msresample:ignoredRange'))
                    else
                        R = pval;
                        if numel(R)~=2 || diff(R)<=0  
                            error(message('bioinfo:msresample:badRange'))
                        end
                        R(1) = max(R(1),0);
                    end
                case 2 % uniform
                    if xoutgiven 
                        warning(message('bioinfo:msresample:ignoredUniform'))
                    else
                        evenlySpaced = bioinfoprivate.opttf(pval,pname,mfilename);
                    end
                case 3 % missing
                    checkMissing = bioinfoprivate.opttf(pval,pname,mfilename);
                case 4 % window
                    windowTypes = {'flattop','blackman','hamming','hanning','kaiser'};
                    h = find(strncmpi(pval,windowTypes,numel(pval)));
                    if isempty(h)
                        error(message('bioinfo:msresample:NotValidWindow'))
                    end
                    windowType = windowTypes{h};
                 case 5 % show
                    if bioinfoprivate.opttf(pval) 
                        if isnumeric(pval)
                            if isscalar(pval)
                                plotId = double(pval); 
                            else
                                plotId = 1;
                                warning(message('bioinfo:msresample:SPNoScalar'))
                            end
                        else
                            plotId = 1;
                        end
                    else
                        plotId = 0;
                    end        
                 case 6 % cutoff
                    F = pval(1);
                    cutoffGiven = true;
                    if F<=0 || F>1
                        error(message('bioinfo:msresample:InvalidCutoffFrequency'))
                    end
                 case 7 % norangewarn
                    rangewarn = ~bioinfoprivate.opttf(pval,pname,mfilename);
            end
        end
    end
end

if rangewarn && (R(1)<min(X) || R(2)>max(X))
    warning(message('bioinfo:msresample:extrapolation', sprintf( '%f', min( X ) ), sprintf( '%f', max( X ) )))
end

% make sure there are no dropped samples 
if checkMissing
    [X,Y] = msdroppedsamples(X,Y);
end

% figure out the new X vector
if ~xoutgiven
if evenlySpaced
    XOUT = R(1)+(0:(N-1))'*diff(R)/N;
else % robust polyfit of observed spaces between samples
    nX = cast(numel(X),class(X));
    xx = (2-nX:2:nX-2)./(nX-2); % centered and scaled domain to improve the 
                                % numerical properties of the fitting 
    dx = diff(X)';              % observed spaces to regress
    h = true(size(dx));         % logical vector to exclude outliers
    in = inf;
    count = 0;
    while (sum(h)~=in) && (count<10)
        count = count + 1;
        in = sum(h);
        P = polyfit(xx(h),dx(h),1);
        res = (polyval(P,xx)-dx).^2; % residuals
        h = res<=2*mean(res);  % outliers are inside two dev stds
    end
    % find the intersample distance after the robust fit at the start and
    % the end of the range,
    SE = polyval(P,xx([1 end]));
    
    if abs(mean(SE)./diff(SE)) > 1000
        XOUT = R(1)+(0:(N-1))'*diff(R)/N;
    else %spaces increase/decrease linearly
        K1 = diff(R).*2./sum(SE)./(N-1);
        K2 = diff(SE).*K1./(N-2); 
        XOUT = [R(1) + (0:N-2).*(K1.*SE(1)) + ((-1:N-3).*(0:N-2)./2).*K2 , R(2)]';
    end
end
end

if cutoffGiven
    DF = F;
else
    % estimate variable sampling factor (this loop is slightly faster than
    % histc)
    nX = numel(X);
    df = zeros(1,N);
    i = 1; j = 1;
    while j <= N
        count=0;
        while (i <= nX) && (X(i)<XOUT(j))
            i = i + 1;
            count = count+1;
        end
        df(j) = count;
        j = j+1;
    end
    df(1) = df(2); % the first bucket contains all samples to the left of the 
                   % range 
    
    % 'df' now contains the observed variable sampling ratio at every X
    % value, to minimize the effect of outliers two steps are taken:
    %    1) eliminate low X (Xmin = 50) values because df can easily be
    %    overestimated due to the small increments when X follows a
    %    quadratic law. 
    %    2) fit df to a low order polynomial (polOrder = 10), so we can
    %    figure out the maximum downsampling rate without taking a likely
    %    outlier. 
    
    Xmin = 50;
    polOrder = 10;
    startAt = floor(min(N/5,find(XOUT>Xmin,1,'first'))); 
    DF = 1/max(polyval(polyfit((startAt:N)/N,df(startAt:N),polOrder),(startAt:N)/N));
    
    % this plot is only for debugging purposes (not documented)
    plotFlag = false;
    if plotFlag 
        figure %#ok<UNRCH>
        plot(XOUT,df,'r'); 
        hold on
        plot(XOUT,polyval(polyfit((startAt:N)/N,df(startAt:N),polOrder),(1:N)/N))
        title('Variable downsampling ratio')
        xlabel('Separation Units')
        hold off
    end
end

if DF >= 1  % it is UPSAMPLE, then just interpolate
    YOUT = interp1(X,Y,XOUT,'spline',0);  
else % DF < 1  it is DOWNSAMPLE, do antialias filter before interpolate
    filtLength = max(21,round(1/DF)*4+1);
    YF = filter(msfir(filtLength,DF,windowType),1,[Y;zeros(filtLength-1,size(Y,2))]);
    YF = YF(1+(filtLength-1)/2:end-(filtLength-1)/2,:);
    YOUT = interp1(X,YF,XOUT,'linear',0);
end

numSignals = size(Y,2);
if (plotId~=0) && ~any((1:numSignals)==plotId)
    warning(message('bioinfo:msresample:InvalidPlotIndex'))
elseif plotId > 0
    figure
    plot(X,Y(:,plotId),'.')
    hold on
    plot(XOUT,YOUT(:,plotId),'r.-')
    title(sprintf('Signal ID: %d  Cutoff Freq: %f ',plotId,DF))
    xlabel('Separation Units')
    ylabel('Relative Intensity')
    legend('Original samples','Up/down-sampled signal')
    axis([min(X) max(X) min(Y(:,plotId)) max(Y(:,plotId))])
    setAllowAxesRotate(rotate3d(gcf),gca,false)
    grid on
    hold off
end

if nargout == 0
    clear XOUT
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = msfir(L,F,windowType)
% Designs a low pass anti-aliasing filter. For more advance filter design
% options see the Signal Processing Toolbox.
% Length (L) is the size of filter.
% Cut-off frequency (F) is a value between [0 1] using normalized
% frequencies. 

% Reference for window coefficients: Bruel & Kjaer, Windows to FFT Analysis
%                                   (Part I), Technical Review, No. 3, 1987 

N = L-1; % filter order

% f = firls(N,[0;F;F;1],[1;1;0;0]);
if rem(L,2)==0 % Type II linear phase filter (even)
    x = pi/2 * F * (1:2:N);
    f = F * sin(x)./x;
    f = [fliplr(f) f];
else % Type I linear phase filter (odd)
    x = pi * F * (1:N/2);
    f = F * sin(x)./x;    
    f = [fliplr(f) F f];
end

x = (2*pi/N)*(0:N);
switch windowType
    case 'hanning' % Hanning window
        w = 0.5 - 0.5*cos(x);
    case 'hamming'  % Hamming window
        w = 0.54 - 0.46*cos(x);
    case 'blackman' % Blackman window 
        w = 0.42 - 0.5*cos(x) + 0.08*cos(2*x);
    case 'flattop' % Flattop window
        w = 0.2156 - 0.4160*cos(x) + 0.2781*cos(2*x) - 0.0836*cos(3*x) +...
            0.0069*cos(4*x);
    case 'kaiser'  % Kaiser window (undocumented)
        w = kaiser(N+1,7.8562)'; % only valid when SPT is present 
end

f = f.* w; % Apply windowing
f = f/sum(f);

    
    
 
    
    
