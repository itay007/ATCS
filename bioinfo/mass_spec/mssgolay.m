function Y = mssgolay(X,Y,varargin)
%MSSGOLAY provides least-squares polynomial smoothing of signals with peaks
%
%   YOUT = MSSGOLAY(X,Y) smoothes the signal Y using the least-squares
%   digital polynomial filters (Savitzky and Golay filters). The default
%   span (or frame) is 15 samples.
%
%   X and Y are column vectors where paired values represent points in the
%   signal. Y can be a matrix with several signals, all sharing the same X
%   scale. Units in the X scale (separation units or s.u.) may quantify
%   wavelength, frequency, distance, time or m/z depending on the type of
%   instrument that generates the signal. 
%
%   MSSGOLAY(...,'SPAN',S) allows you to modify the frame size for the
%   smoothing function. If S is greater than 1, then the window is of size
%   SP (in samples, i.e., independently of the X vector). Higher values
%   will smooth the signal more, at the expense of computation time. If S
%   is less than 1, then the window size is taken to be a fraction of the
%   number of points in the data; e.g., for S = 0.05, the window size equal
%   to 5% of the number of points in X. 
%
%   Notes:
%   1) The original Savitzky and Golay algorithm assumes a uniformly
%   spaced X vector; MSSGOLAY also allows one that is not uniformly
%   spaced. Therefore, the sliding frame for smoothing is centered using
%   the closest samples in terms of the X value and not in terms of the X
%   index.
%   2) When the X vector does not have repeated values or NaNs, the
%   algorithm will be approximately two times faster.
%   3) When the X vector is evenly spaced, the least-squares fitting is
%   performed once and the signal is filtered with the same coefficients,
%   speeding up the algorithm considerably.
%   4) If the X vector is evenly spaced and S is even, S is automatically
%   incremented by 1 to include both edge samples in the current frame. 
%
%   MSSGOLAY(...,'DEGREE',D) sets the D of the polynomial to be fitted to
%   all the points in the moving frame. The default DEGREE is 2. D must be
%   smaller than S. 
%
%   MSSGOLAY(...,'SHOWPLOT',SP) plots the smoothed signal over the
%   original. When SP is TRUE, the first signal in Y is used. The
%   default is FALSE. No plot is generated unless MSSGOLAY is called
%   without output arguments. SP can also contain an index to one of the
%   signals in Y.  
%
%   Example: 
%
%       load sample_lo_res
%
%       % Smooth a group of spectrograms.
%       YS = mssgolay(MZ_lo_res,Y_lo_res);
%
%       % Plot the third spectrogram in Y_lo_res and its smoothed signal.
%       mssgolay(MZ_lo_res,Y_lo_res,'SHOWPLOT',3);
%
%   See also MSALIGN, MSBACKADJ, MSHEATMAP, MSLOWESS, MSNORM, MSPREPRODEMO,
%   MSRESAMPLE, MSVIEWER.

% Copyright 2003-2008 The MathWorks, Inc.



% set default parameters
span = 15;
degree = 2; 
if nargout == 1
    plotId = 0; 
else
    plotId = 1;
end

% deal with the various inputs
if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:mssgolay:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'span','degree','showplot'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:mssgolay:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:mssgolay:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % span
                    if ~isnumeric(pval) || ~isscalar(pval)
                        error(message('bioinfo:mssgolay:LowessSpanNumeric'));
                    end
                    span = pval;
                    if span <= 0
                        error(message('bioinfo:mssgolay:LowessSpanPositive'));
                    end
                    if span >= numel(X)
                        error(message('bioinfo:mssgolay:LowessSpanLessThanX'));
                    end
                case 2 % degree
                    if ~isnumeric(pval) || ~isscalar(pval)
                        error(message('bioinfo:mssgolay:DegreeNumeric'));
                    end
                    degree = pval;
                    if degree<0
                        error(message('bioinfo:mssgolay:LowessBadNumericDegree'));
                    end
                case 3 % show
                    if bioinfoprivate.opttf(pval) 
                        if isnumeric(pval)
                            if isscalar(pval)
                                plotId = double(pval); 
                            else
                                plotId = 1;
                                warning(message('bioinfo:mssgolay:SPNoScalar'))
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
if size(X,2) ~= 1
   error(message('bioinfo:mssgolay:OnlyOneXScale'))
end
if ~isnumeric(Y) || ~isreal(Y)
   error(message('bioinfo:mssgolay:IntensityNotNumericAndReal')) 
end
if ~isnumeric(X) || ~isreal(X)
   error(message('bioinfo:mssgolay:XNotNumericAndReal')) 
end
if size(X,1) ~= size(Y,1)
   error(message('bioinfo:mssgolay:NotEqualNumberOfSamples'))
end

numSignals = size(Y,2);

if (plotId~=0) && ~any((1:numSignals)==plotId)
    warning(message('bioinfo:mssgolay:InvalidPlotIndex'))
end

% iterate for every signal
for ns = 1:numSignals
if nargout>0 || (ns == plotId)
    
    if (ns == plotId)
        yo = Y(:,ns);
    end
    
    Y(:,ns) = bioinfoprivate.masmooth(X,Y(:,ns),span,'sgolay',degree);
        
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

