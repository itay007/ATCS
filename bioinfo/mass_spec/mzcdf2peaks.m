function [P,T] = mzcdf2peaks(S)
% MZCDF2PEAKS converts a mzCDF struct to a list of peaks.
%
%   [P,T] = MZCDF2PEAKS(MZCDFSTRUCT) extracts peak information from an
%   mzCDF structure. P is a cell with peak lists, a peak list is a 2 column 
%   matrix with the mass/charge (MZ) and ion intensity for every peak. T is
%   a vector with the retention time of every scan.
%
%   Example:
%
%       % Extract the peak information from the file results.cdf and create
%       % a dot plot:
%       mzcdf_struct = mzcdfread('results.cdf');
%       [peaks,time] = mzcdf2peaks(mzcdf_struct);
%       msdotplot(peaks,time)
%
%   Note that the file results.cdf is not provided.
%
%   See also MSPALIGN, MSDOTPLOT, MSPPRESAMPLE, MZCDFREAD.

%   Copyright 2007 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);

if ~isfield(S,'intensity_values') || ~isfield(S,'mass_values') || ~isfield(S,'scan_index')
    error(message('bioinfo:mzcdf2peaks:invalidFieldnames'))
end

I = S.scan_index+1;
J = [S.scan_index(2:end);numel(S.mass_values)];

n = numel(I);

P = cell(n,1);

for i =1:n
    P{i} = [single(S.mass_values(I(i):J(i)))  single(S.intensity_values(I(i):J(i)))];
end

if nargout>1
    if ~isfield(S,'scan_acquisition_time') 
      error(message('bioinfo:mzcdf2peaks:invalidFieldnames'))
    end
    T = S.scan_acquisition_time(:);
end
