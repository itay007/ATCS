function [Y,P] = msnorm(X,Y,varargin)
%MSNORM normalize a set of signals with peaks
%
%   YOUT = MSNORM(X,Y) normalizes a group of signals with peaks by
%   standardizing the area under the curve (AUC) to the group median. 
%
%   X and Y are column vectors where paired values represent points in the
%   signal. Y can be a matrix with several signals, all sharing the same X
%   scale. Units in the X scale (separation units or s.u.) may quantify
%   wavelength, frequency, distance, time or m/z depending on the type of
%   instrument that generates the signal. 
%
%   MSNORM(...,'QUANTILE',QV) sets a 1-by-2 vector with the quantile
%   limits to select a reduced set of X positions. Default is [0 1], i.e.
%   it uses the whole AUC. For example, when QV = [0.9 1], only the largest
%   10% of intensities in every signal are used to compute the AUC. When QV
%   is a scalar it represents the lower quantile limit and the upper
%   quantile limit is set to 1.
%
%   MSNORM(...,'LIMITS',LIM) sets a 1-by-2 vector with the valid X range
%   to pick points for normalization. This parameter is useful to eliminate
%   noise from the AUC calculation, for example the matrix noise that
%   appears in the low-mass region of SELDI mass spectrometers.
%
%   MSNORM(...,'CONSENSUS',CON) sets a consensus rule; to include a X
%   position into the AUC its intensity must be within the quantile limits
%   of at least a CON proportion of the signals in Y. The same X positions
%   are used to normalize all the signals. CON is a scalar between 0 and 1.
%
%   MSNORM(...,'METHOD',METHOD) sets the target to which the AUC of every
%   signal is normalized. METHOD can be 'Median' (default) or 'Mean'. 
%
%   MSNORM(...,'MAX',MAX). After individually normalizing every signal,
%   they are further scaled to adjust the overall maximum intensity to MAX.
%   MAX is a scalar. If it is omitted, no postscaling is performed. 
%
%   MSNORM(X,Y,P) employs the information in P to normalize a new set of
%   signals using the same parameters to select the X positions and output
%   scale as a previous normalization. P is a structure created by MSNORM.
%   If a consensus proportion was given in the previous normalization, no
%   new X positions are selected; normalization is performed using the
%   the same X positions. 
%
%   [YOUT,P] = MSNORM(...) returns a structure containing the necessary
%   parameters to normalize another set of signals.
%
%   Examples: 
%   
%      load sample_lo_res
%      Y = Y_lo_res(:,[1 2 5 6]);
%      MZ = MZ_lo_res;
%
%      % Normalize the AUC of every spectrogram to the median.
%      Y1 = msnorm(MZ,Y);
%
%      % Normalize the AUC of every spectrogram to the median, eliminating
%      % low-mass (m/z < 1,000) noise, and postscaling such that the
%      % maximum intensity is 50.
%
%      Y2 = msnorm(MZ,Y,'LIMITS',[1000 inf],'MAX',50);
%
%      % Normalize the ion intensities of all spectra to the maximum
%      % intensity of the single highest peak from any of the spectra in 
%      % the range above 1000 MZ.
%
%      Y3 = msnorm(MZ,Y,'QUANTILE',[1 1],'LIMITS',[1000 inf]);
%
%      % Select MZ regions whose intensities are within the third quartile
%      % in at least 90% of the spectrograms.
%
%      [Y4,S] = msnorm(MZ,Y,'QUANTILE',[0.5 0.75],'CONSENSUS',0.9);
%
%      % Use the same MZ regions to normalize another set of spectrograms.
%
%      Y5 = msnorm(MZ,Y,S);
%
%   See also MSALIGN, MSBACKADJ, MSHEATMAP, MSLOWESS, MSPREPRODEMO,
%   MSRESAMPLE, MSSGOLAY, MSVIEWER. 

%   Copyright 2003-2012 The MathWorks, Inc.


% References:
% [1] M. Wagner, D. Nalk, and A. Pothen, "Protocols for disease
%     classification from mass spectrometry data" Proteomics 3,pp.1692-1698
%     (2003).
% [2] G. Satten, et al in "Standardization and denoising algorithms for
%     mass spectra to classify whole-organism bacterial specimens"
%     Bioinformatics 2004 Jun 24 (Epub ahead of print)
% [3] Leping Li, David M. Umbach, Paul Terry and Jack A. Taylor,
%     "Application of the GA/KNN method to SELDI proteomics data" PNAS 2003.
% [4] Lilien et al in "Probabilistic Disease Classification of
%     Expression-Dependent Proteomic Data from Mass Spectrometry of Human
%     Serum", Journal of Computational Biology, 10(6) 2003, pp. 925-946.

% check inputs
bioinfochecknargin(nargin,2,mfilename);
% set defaults
parametersGiven = false;
outputScaleGiven = false;
consensusRateGiven = false;
methodIsMedian = true;
quantileValues = [0 1];
consensusRate = 0.75;
XLimits = [0,inf];

% validate X and Y
if ~isnumeric(Y) || ~isreal(Y)
   error(message('bioinfo:msnorm:IntensityNotNumericAndReal')) 
end

if ~isnumeric(X) || ~isreal(X) || ~isvector(X)
   error(message('bioinfo:msnorm:XNotNumericAndReal')) 
end

if size(X,1) ~= size(Y,1)
   error(message('bioinfo:msnorm:NotEqualNumberOfSamples'))
end

[numSamples, numSignals] = size(Y);

nvarargin =  numel(varargin); 
% check if third input is the structure with parameters
if nvarargin == 1 && isstruct(varargin{1})
    P = varargin{1};
    if ~isfield(P,'targetedSumI') || ~isfield(P,'Xh')
            error(message('bioinfo:msnorm:NoMsnormStruct'));
    end
    if ~isnan(P.Xh) % it is a normalization with consensus
        if numel(P.Xh)~=numSamples || ~islogical(P.Xh)
            error(message('bioinfo:msnorm:InconsistentMsnormStruct'));
        end
    end
    parametersGiven = true;
else
    if nvarargin
        if rem(nvarargin,2)
            error(message('bioinfo:msnorm:IncorrectNumberOfArguments', mfilename));
        end
        okargs = {'quantile','consensusrate','limits','max','method'};
        for j=1:2:nvarargin
            pname = varargin{j};
            pval = varargin{j+1};
            if ~ischar(pname)
                error(message('bioinfo:msnorm:ParameterNameNochar'));
            end
            k = find(strncmpi(pname, okargs,length(pname)));
            if isempty(k)
                error(message('bioinfo:msnorm:UnknownParameterName', pname));
            elseif length(k)>1
                error(message('bioinfo:msnorm:AmbiguousParameterName', pname));
            else
                switch(k)
                    case 1 % 'quantile'
                        if isscalar(pval)
                            if pval<0 || pval>1
                                error(message('bioinfo:msnorm:badScalarQuantile'))
                            else
                                quantileValues = [pval 1];
                            end
                        else
                            if numel(pval)~=2 || diff(pval)<0
                                error(message('bioinfo:msnorm:badVectorQuantile'))
                            else
                                quantileValues = [pval(1) pval(2)];
                            end
                        end
                    case 2 % 'consensusrate'
                        if ~isscalar(pval) || pval<0 || pval>1
                            error(message('bioinfo:msnorm:badConsensus'))
                        end
                        consensusRate = pval;
                        consensusRateGiven = true;
                    case 3 % 'limits'
                        if numel(pval)~=2 || diff(pval)<=0
                            error(message('bioinfo:msnorm:badRange'))
                        end
                        XLimits = [max(0,max(pval(1),X(1))),min(pval(2),X(end))];
                    case 4 % 'max'
                        if ~isscalar(pval)
                            error(message('bioinfo:msnorm:badScale'))
                        end
                        outputScaleGiven = true;
                        outputScale = pval;
                    case 5 % 'method'
                        if ischar(pval)
                            okmethods = {'mean','median'};
                            nm = strmatch(lower(pval), okmethods);
                            if isempty(nm)
                                error(message('bioinfo:msnorm:MethodNameNotValid'));
                            elseif length(nm)>1
                                error(message('bioinfo:msnorm:AmbiguousMethodName', pval));
                            else
                                methodIsMedian = nm == 2;
                            end
                        else
                             error(message('bioinfo:msnorm:MethodNameNotValid'));
                        end
                end
            end
        end
    end
end

if ~parametersGiven

    % search for the selected X values
    Xl = X>=XLimits(1) & X<=XLimits(2);
    Q = zeros(numSignals,2);
    % loop so we avoid more copies of the huge Y
    for ind = 1:numSignals
        Q(ind,:) = quantile(Y(Xl,ind),quantileValues);
    end
    Q = Q';
    if consensusRateGiven
        msCount = zeros(numSamples,1);
        for k = 1:numSignals
            t = Y(:,k)>=Q(1,k) & Y(:,k)<=Q(2,k) & Xl;
            if ~sum(t)
                [dummy,h] = min(abs(Y(:,k)-Q(1,k)) + inf.*~Xl);  
                msCount(h) = msCount(h)+1;
            else
                msCount = msCount+t;
            end
        end % for k = ...
        Xh = msCount/numSignals >= consensusRate;
        if ~sum(Xh)
            error(message('bioinfo:msnorm:SelectedXisEmpty'))
        end
        % compute the sum of the intensities of selected points
        sumI = sum(Y(Xh,:));
    else %~consensusRateGiven
        sumI = zeros(numSignals,1);
        for k = 1:numSignals
            t = Y(:,k)>=Q(1,k) & Y(:,k)<=Q(2,k) & Xl;
            if ~sum(t)
                [dummy,h] = min(abs(Y(:,k)-Q(1,k)) + inf.*~Xl);  
                sumI(k) = Y(h,k);
            else
                sumI(k) = sum(Y(t,k));
            end
        end
        Xh = NaN; %will indicate it is a normalization does not use consensus
    end

    if methodIsMedian
        medSumI = median(sumI);
    else
        medSumI = mean(sumI);
    end

    % scale every signal
    for k = 1:numSignals
        Y(:,k) = Y(:,k) * (medSumI/sumI(k));
    end
 
    if outputScaleGiven
        K = outputScale/max(max(Y(Xl,:)));
        Y = Y * K;
        targetedSumI = K * medSumI;
    else
        targetedSumI = medSumI;
    end
    if nargout > 1
        P.Xh = Xh;
        P.targetedSumI = targetedSumI;
        if ~consensusRateGiven
            P.XLimits = XLimits;
            P.quantileValues = quantileValues;
        end
    end
else % parametersGiven
    if isnan(P.Xh) % it is a normalization without consensus
        Xl = X>=P.XLimits(1) & X<=P.XLimits(2);
        Q = zeros(numSignals,2);
        for ind = 1:numSignals
            Q(ind,:) = quantile(Y(Xl,ind),P.quantileValues);
        end
        Q = Q';
        sumI = zeros(numSignals,1);
        for k = 1:numSignals
            t = Y(:,k)>=Q(1,k) & Y(:,k)<=Q(2,k) & Xl;
            if ~sum(t)
                [dummy,h] = min(abs(Y(:,k)-Q(1,k)) + inf.*~Xl);  
                sumI(k) = Y(h,k);
            else
                sumI(k) = sum(Y(t,k));
            end
        end
    else % it is a normalization with consensus, X positions already chosen
        sumI = sum(Y(P.Xh,:));
    end
    % scale every signal
    for k = 1:numSignals
        Y(:,k) = Y(:,k) * (P.targetedSumI/sumI(k));
    end
end

if nargout == 0 
    clear Y
end
