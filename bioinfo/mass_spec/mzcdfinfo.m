function out = mzcdfinfo(filename)
%MZCDFINFO Information about netCDF file.
%   INFO = MZCDFINFO(FILENAME) returns a structure whose fields contain
%   information about a netCDF file with mass spectrometric data.  FILENAME
%   is a string containing a file name, or a path and a file name, of a
%   netCDF file that conforms to the ANDI/MS or the ASTM E2077-00 (2005)
%   standard specification or earlier specifications. INFO is a structure
%   with the following fields: 
%
%       Filename - Name of the file
%
%       FileTimeStamp - NetCDF file date time stamp
%
%       FileSize - Size of the file in bytes
%
%       NumberOfScans - Number of scans in the file
%
%       StartTime - Run start time
%
%       EndTime - Run end time
%
%       TimeUnits - Units for time
%
%       GlobalMassMin - Minimum M/Z value in all scans
%
%       GlobalMassMax - Maximum M/Z value in all scans
%
%       GlobalIntensityMin - Minimum intensity in all scans
%
%       GlobalIntensityMax - Maximum intensity in all scans
%
%       ExperimentType - Indicates if raw or centroided data
%
%       NOTE: if optional attributes are not in the netCDF file the field
%       value will be set to 'N/A' or NaN.
%
%   Example:
%
%       % Get netCDF file information.
%       info = mzcdfinfo('results.cdf');
%
%   Note that the file results.cdf is not provided. 
%
%   See also JCAMPREAD, MZCDF2PEAKS, MZCDFREAD, MZXML2PEAKS, MZXMLINFO,
%   MZXMLREAD, TGSPCINFO, TGSPCREAD.

% Copyright 2008 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);

if( ~ischar( filename ) )
    error(message('bioinfo:mzcdfinfo:invalidUsage'));
end
% Check if file exists
if ~exist(filename,'file')
     error(message('bioinfo:mzcdfinfo:FileNotFound', filename));
end

% opens the netcdf for reading
nc_global = netcdf.getConstant('NC_GLOBAL');
cdfid = netcdf.open(filename,'NOWRITE');
cdfidc = onCleanup(@()netcdf.close(cdfid)); 

out.Filename = filename;

[ndims, nvars, natts] = netcdf.inq(cdfid);

attName = cell(natts,1);
attValue = cell(natts,1);
% read global attributes
for i = 1:natts
    attName{i} = netcdf.inqAttName(cdfid,nc_global,i-1);
    attValue{i} = netcdf.getAtt(cdfid,nc_global,attName{i});
end

h = strmatch('netcdf_file_date_time_stamp',attName,'exact');
if isempty(h)
    out.FileTimeStamp = 'N/A';
else
    out.FileTimeStamp = attValue{h};
end

tmp = dir(filename);
out.FileSize = tmp.bytes;

% read variable names:
varName = cell(nvars,1);
for i = 1:nvars
    varName{i} = netcdf.inqVar(cdfid,i-1);
end
h = strmatch('scan_acquisition_time',varName,'exact')-1;
if isempty(h)
    out.NumberOfScans = NaN;
    out.StartTime = NaN;
    out.EndTime = NaN;
else
    T =  netcdf.getVar(cdfid,h,'double');
    out.NumberOfScans = numel(T);
    out.StartTime = min(T);
    out.EndTime = max(T);
end

h = strmatch('units',attName,'exact');
if isempty(h)
    out.TimeUnits = 'N/A';
else
    out.TimeUnits = attValue{h};
end

h = strmatch('global_mass_min',attName,'exact');
if isempty(h)
    out.GlobalMassMin = NaN;
else
    out.GlobalMassMin = attValue{h};
end

h = strmatch('global_mass_max',attName,'exact');
if isempty(h)
    out.GlobalMassMax = NaN;
else
    out.GlobalMassMax = attValue{h};
end

h = strmatch('global_intensity_min',attName,'exact');
if isempty(h)
    out.GlobalIntensityMin = NaN;
else
    out.GlobalIntensityMin = attValue{h};
end

h = strmatch('global_intensity_max',attName,'exact');
if isempty(h)
    out.GlobalIntensityMax = NaN;
else
    out.GlobalIntensityMax = attValue{h};
end

h = strmatch('experiment_type',attName,'exact');
if isempty(h)
    out.ExperimentType = 'N/A';
else
    out.ExperimentType = attValue{h};
end


