function [YD,YN] = waveletdenoise(Y,waveletBasis,waveletLevels,noiseEstimator,thresholdMultiplier)
% WAVELETDENOISE Implements undecimated wavelet denoise with the Daubechies filters.
%
%  [YD,YN] = WAVELETDENOISE(Y,WAVELETBASIS,WAVELETLEVELS,NOISEESTIMATOR,THRESHOLDMULTIPLIER)

%   Copyright 2006-2007 The MathWorks, Inc.


%[1] Donoho, D.L. and Johnstone, I.M. "Adapting to unknown smoothness via
%    wavelet shrinkage," J. Am. Statist. Asso. 90 pp. 1200-1224 (1995). 
%[2] Gilbert Strang and Truong Nguyen, "Wavelets and Filter Banks"
%    Wellesley Cambridge Press (1996). 

% set defaults
if nargin<5
    thresholdMultiplier = 1;
    if nargin<4
        noiseEstimator = 'mad'; % std, none
        if nargin<3
            waveletLevels = 10;              
            if nargin<2
                waveletBasis = 4;
            end
        end
    end
end
            
persistent h0 h1 f0 f1 perwb perwl
if isempty(h0) || perwb~=waveletBasis || perwl~=waveletLevels
    perwb = waveletBasis;
    perwl = waveletLevels;
    [h0 h1 f0 f1] = mydau(waveletBasis,2.^(waveletLevels-1));
end

if ~isnumeric(Y) || ~isreal(Y)
   error(message('bioinfo:waveletdenoise:IntensityNotNumericAndReal')) 
end

if any(isnan(Y))
    replacenans = true;
    originalnans = isnan(Y);
    Y(originalnans) = 0;
else
    replacenans = false;
end


[n,numSpectrograms] = size(Y);
YD = zeros(n,numSpectrograms); % allocating output

% defining some constants a needed variables
m = waveletBasis*2-1;        % basis length
scs = 2.^(0:waveletLevels);  % scales
bll = ceil((n+cumsum([0 scs(1:waveletLevels)].*m))./scs) + m; % block length with lags
ys = ceil((sum(scs(1:waveletLevels)*m)+n)./scs(waveletLevels+1)).*fliplr(scs);

yh = cell(waveletLevels,1);   
yl = zeros(1,ys(1));

for j = 1:numSpectrograms

    % wavelet decomposition without decimation
    yl(:) = 0;
    yl(1:n) = Y(:,j);
    for i = 1:waveletLevels-1
        yh{i} = filter(h1,1,yl(:,1:bll(i)),[],2);
        yl(:,1:bll(i)) = filter(h0,1,yl(:,1:bll(i)),[],2);
        yl = reshape(yl,scs(i+1),ys(i+1));
    end
    yh{waveletLevels} = filter(h1,1,yl(:,1:bll(waveletLevels)),[],2);
    yl(:,1:bll(waveletLevels)) = filter(h0,1,yl(:,1:bll(waveletLevels)),[],2);

    switch noiseEstimator
        case 'mad' %VisuShrink
            th = thresholdMultiplier * sqrt(2*log(n)) * mad(yh{1}(1:n),1)/.6745;
        case 'std'
            th = thresholdMultiplier * std(yh{1}(1:n));
        otherwise % 'none'
            th = thresholdMultiplier;
    end

    % soft thresholding  (Strang & Nguyen 1996, pp387)
    yh = cellfun(@(w) max(0,(abs(w)-th)).*sign(w),yh,'Uniform',false);
    yl = max(0,(abs(yl)-th)).*sign(yl);

    % hard thresholding  (Strang & Nguyen 1996, pp387)
    % yh = (abs(yh)>th).*yh;
    % yl = (abs(yl)>th).*yl;

    %wavelet reconstruction without upsampling
    for i = waveletLevels:-1:2
        vl = filter(f0,1,yl(:,1:bll(i)),[],2);
        vh = filter(f1,1,yh{i},[],2);
        yl(:,1:bll(i)-m) = (vl(:,m+1:end)+vh(:,m+1:end))/2;     
        yl = reshape(yl,scs(i-1),ys(i-1));
        
    end
           
    vl = filter(f0,1,yl(:,1:bll(1)),[],2);
    vh = filter(f1,1,yh{1},[],2);
    yl(:,1:bll(1)-m) = (vl(:,m+1:end)+vh(:,m+1:end))/2; 
    
    YD(:,j) = yl(1:n);

    %plot([YD(:,j) Y(:,j)]) % verify perfect reconstruction
    %pause
    
end

if replacenans 
    Y(originalnans) = NaN;
end

if nargout>1 %noise signal
    YN = Y - YD;
end

function [h0 h1 f0 f1] = mydau(k,u)
% Calculates Daubechies coefficients for an undecimated transform
a = bsxfun(@rdivide,0.5-k:k,toeplitz([inf 1-(2:k)],[inf 1:k+k-1]));
P(1:2:k+k) = prod(~a+a,2);
R = roots([P,1,fliplr(P)]); 
R = R(abs(R)<1);
[dummy,i] = sort(abs(R+1),'descend');
h = real(poly([R(i(1:k-1));-ones(k,1)]));
h = h.*(sqrt(2)/sum(h)); 
f0 = h;
h1 = h.*((-1).^(1:k+k));
h0 = fliplr(f0);
f1 = fliplr(h1);
