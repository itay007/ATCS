function [MZ,Y] = msbatchprocessing(repository,files)
% MSBATCHPROCESSING Example function for BIODISTCOMPDEMO
%
% [MZ,Y] = MSBATCHPROCESSING(REPOSITORY,FILES) Preprocesses the
% spectrogram in files FILES and returns the mass/charge (MZ) and ion
% intensities (Y) vectors.
%
% Hard-coded parameters in the preprocessing steps have been adjusted to
% deal with the high-resolution spectrograms of the example.

% Copyright 2004-2012 The MathWorks, Inc.

K = numel(files);
Y = zeros(15000,K); % need to preset the size of Y for memory performance
MZ = zeros(15000,1);

parfor k = 1:K

    file = [repository files{k}]; 
    
    % read the two-column text file with mass-charge and intensity values
    fid = fopen(file,'r');
    ftext = textscan(fid, '%f%f');
    fclose(fid);
    signal = ftext{1};
    intensity = ftext{2};

    % resample the signal to 15000 points between 2000 and 11900
    mzout = (sqrt(2000)+(0:(15000-1))'*diff(sqrt([2000,11900]))/15000).^2;
    [mz,YR] = msresample(signal,intensity,mzout);
    
    % align the spectrograms to two good reference peaks
    P = [3883.766 7766.166];
    YA = msalign(mz,YR,P,'WIDTH',2);

    % estimate and adjust the background
    YB = msbackadj(mz,YA,'STEP',50,'WINDOW',50);

    % reduce the noise using a nonparametric filter
    Y(:,k) = mslowess(mz,YB,'SPAN',5);
    
    % the mass/charge vector is the same for all spectra after the resample
    if k==1
        MZ(:,k) = mz;
    end

end
