function [nmz,y] = msdroppedsamples(mz,y)
%MSDROPPEDSAMPLES analyzes a MZ vector for dropped samples
%
%   NMZ = MSDROPPEDSAMPLES(MZ) Analyzes a MZ vector to check for dropped
%   samples, and inserts them if needed. This function assumes that the
%   sampling period increases linearly or is constant, which means the MZ
%   vector is either a quadratic or linear function of the sample index,
%   respectively. Recovering missing samples is important; if dropped
%   samples are not recovered, high mass/charge peaks can show a
%   substantial shift compared to their real position.
%
%   [NMZ,NY] = MSDROPPEDSAMPLES(MZ,Y) If an ion intensity vector is
%   provided, the missing values are interpolated. 
%
%   Example:
%
%      load sample_hi_res
%      plot(diff(mz),'.')
%      mz =  msdroppedsamples(mz);
%      hold on
%      plot(diff(mz),'.r')

% Copyright 2003-2006 The MathWorks, Inc.


if (nargout == 2) && (nargin == 1)
    error(message('bioinfo:msdroppedsamples:noIntensityVectorSupplied'))
end

if (nargin > 1) && (size(mz,1)~=size(y,1))
    error(message('bioinfo:msdroppedsamples:differentSize'))
end

if ~isvector(mz)
    error(message('bioinfo:msdroppedsamples:NotAVector'))
end

transposeMZ = false;
if size(mz,1)<size(mz,2)
    transposeMZ = true;
    mz = mz';
end

n = numel(mz);

% try first a uniform sampling period
Bo = diff(mz);
Be = quantile(Bo,0.5);        % first estimation assuming roughly less 
                              % than 50% samples are lost
x = Bo(Bo < Be * 1.5);
Be = mean(x);
if sum(abs((x/Be-1))>.05)>n*.05    
    % mz was not evenly spaced then assuming a quadratic period
    Bo = diff(mz)./sqrt(mz(2:n)); % B observed
    Be = quantile(Bo,0.5);    % first estimation assuming roughly less 
                              % than 50% samples are lost   
    x = Bo(Bo < Be * 1.5);                              
    Be = mean(x);             % refining the estimation
    if sum(abs((x/Be-1))>.05)>n*.05   
        % could not recover law governing mz
        error(message('bioinfo:msdroppedsamples:BadMZvector'))
    end
end
  
h = Bo >= Be * 1.5;           % find sites where samples are missing

if sum(h) > n/2
    error(message('bioinfo:msdroppedsamples:TooManyDrops'))
end

m = max(0,round(Bo/Be));      % multiplicity of samples
N = sum(m)+1;                 % size of new MZ vector
nmz = nan(N,1,class(mz));     % allocate space for the output
g = cumsum([1;m(:)]);         % indices of MZ values in new MZ vector
nmz(g) = mz;                  % 
j = find(isnan(nmz));         % indices of missing samples in new MZ vector
nmz(j) = interp1q(g,mz,j);    % interpolate the MZ missing values

if transposeMZ
    nmz = nmz';
end

if nargout>1
    ny = interp1(mz,y,nmz(j)); % interpolate missing intensities
    y(g,:) = y;                 % re-arrange the old intensity vector
    y(j,:) = ny;                % and put in the interpolated values
end

