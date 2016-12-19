function out = mzcdfread(filename,varargin)
% MZCDFREAD Read a netCDF file.
%
%   OUT = MZCDFREAD(FILENAME) reads a netCDF file with mass spectrometric
%   data into MATLAB as a structure. FILENAME is a string containing a file
%   name, or a path and a file name, of a netCDF file that conforms to the
%   ANDI/MS or the ASTM E2077-00 (2005) standard specification or earlier
%   specifications. Fields in the output structure OUT correspond to
%   variables and global attributes in the netCDF file. When a netCDF
%   variable contains local attributes an additional field is created, 
%   the name of the field is composed from the concatenation of the
%   variable name and the string '_attributes'. 
%
%   MZCDFREAD(...,'TIMERANGE', RANGE) reads only the subset of the scans
%   that occurred within a specified time range in the chromatographic (or
%   any other) separation experiment. RANGE is a two-element numeric array.
%   Time units are usually indicated within the netCDF global attributes
%   (see MZCDFINFO). Default is to extract all scans ([0,INF]).
%
%   MZCDFREAD(...,'SCANINDICES',INDICES) reads a subset of the scans
%   specified with a vector of indices. Default is to extract all scans
%   ([1:NumberOfScans]). To indicate a range of indices use
%   [START_IND:END_IND]. To extract a single scan use only the scan index. 
%
%   MZCDFREAD(...,'VERBOSE',T/F) show progress of reading the netCDF file.
%   Default is true.
%
%   Example:
%
%       out = mzcdfread('results.cdf');
%       % view the second scan in the netCDF file:
%       idx1 = out.scan_index(2)+1
%       idx2 = out.scan_index(3)
%       y = out.intensity_values(idx1:idx2);
%       z = out.mass_values(idx1:idx2);
%       stem(z,y,'marker','none')
%
%       title(sprintf('Time: %f',out.scan_acquisition_time(2)))
%       xlabel(out.mass_values_attributes.units)
%       ylabel(out.intensity_values_attributes.units)
%
%   Note that the file results.cdf is not provided.
%
%   See also JCAMPREAD, MZCDFINFO, MZCDF2PEAKS, MZXMLREAD, TGSPCINFO,
%   TGSPCREAD.

% Copyright 2008 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);

if( ~ischar( filename ) )
    error(message('bioinfo:mzcdfread:invalidUsage'));
end
% Check if file exists
if ~exist(filename,'file')
     error(message('bioinfo:mzcdfread:FileNotFound', filename));
end

%Check parameter value pairs
numvaragin = numel(varargin);

%Set default values
tmrange = [];
indices = [];
verbose = true;
getSubset = false;
blcksize = 2^21;

if rem(numvaragin,2)
    error(message('bioinfo:mzcdfread:IncorrectNumberOfArguments', mfilename));
end

% The allowed inputs
okargs = {'timerange','scanindices','verbose'};
for j=1:2:numvaragin
    % Lookup the pair
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % timerange
            if ~isnumeric(pval) || numel(pval)~=2 || diff(pval)<=0
                error(message('bioinfo:mzcdfread:badRange'))
            end
            tmrange = pval(:)';
            getSubset = true;
        case 2  % scanindices
            if ~isnumeric(pval) || ~isvector(pval) || any(pval<=0) || any(rem(pval,1))
                error(message('bioinfo:mzcdfread:InvalidScanIndices'))
            end
            indices = pval(:);
            getSubset = true;
        case 3 %verbose
            verbose = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end

end

% opens the netcdf for reading
nc_global = netcdf.getConstant('NC_GLOBAL');
cdfid = netcdf.open(filename,'NOWRITE');
cdfidc = onCleanup(@()netcdf.close(cdfid)); 

% inquires the netcdf file about the number of dimensions, the number of
% variables, the number of global attributes and the index of the unlimited
% dimension (if any)
[ndims, nvars, natts, recdim] = netcdf.inq(cdfid);
if verbose
 disp([...
  sprintf('\nReading filename: %s  \n',filename)...
  sprintf('Number of dimensions: %d   Unlimited Dimension:  %d\n',ndims,recdim)...
  sprintf('Number of variables:  %d   Number of attributes: %d\n',nvars,natts)])
end

% read global attributes
for i = 1:natts
    attName = netcdf.inqAttName(cdfid,nc_global,i-1);
    attValue = netcdf.getAtt(cdfid,nc_global,attName);
    out.(attName) = attValue;
end

if verbose
    disp('Exploring global netCDF attributes:')
    fnames = fieldnames(out);
    for i = 1:natts
        if isnumeric(out.(fnames{i}))
            disp(sprintf('     %s -> %f',fnames{i},out.(fnames{i})))
        else
            disp(sprintf('     %s -> %s',fnames{i},out.(fnames{i})))
        end
    end
    disp(' ')
end

% explore dimensions
dimName = cell(ndims,1);
dimLength = zeros(ndims,1);
for i = 1:ndims
    [dimName{i},dimLength(i)] = netcdf.inqDim(cdfid,i-1);
end

if verbose
    u = [dimName arrayfun(@(x) {sprintf('%d',x)},dimLength) num2cell(1:ndims)']';
    disp('Exploring netCDF dimensions:')
    disp([sprintf('                          Name  Length  Dim_ID \n') sprintf('%30s %7s    %2d\n',u{:})])
end

pdim = strmatch('point_number',dimName);
if isempty(pdim)
     error(message('bioinfo:mzcdfread:point_NumberMissing'))
end
numP = dimLength(pdim);

% read variable names
varName = cell(nvars,1);
varType = zeros(nvars,1);
varDim = cell(nvars,1);
varNatts = zeros(nvars,1);
for i = 1:nvars
    [varName{i},varType(i),varDim{i},varNatts(i)] = netcdf.inqVar(cdfid,i-1);
end

% read all scans (no subsets)
if ~getSubset
    if verbose
        disp('Reading netCDF variables:')
        disp('                          Name    Type  Dim_ID  NumAtt')
    end
    for i = 1:nvars
        if verbose
            str = sprintf('%d ',varDim{i}+1);
            disp(sprintf('%30s %6d  %6s %6d',varName{i},varType(i),str,varNatts(i)))
        end
        if any(varDim{i}+1 == pdim) % when dimension is point_number read by block
                                    % to make it more memory efficient
            out.(varName{i}) = zeros(numP,1,class(netcdf.getVar(cdfid,i-1,0)));
            for j = 0:blcksize:numP-1
                k = min(j+blcksize,numP);
                out.(varName{i})(j+1:k) = netcdf.getVar(cdfid,i-1,j,k-j);
            end
        else % read the whole variable at once
            out.(varName{i}) = netcdf.getVar(cdfid,i-1);
        end
        if ischar(out.(varName{i})) && isvector(out.(varName{i}))
            out.(varName{i}) = strtrim(strrep(out.(varName{i})(:)',char(0),''));
        end
        for j = 1:varNatts(i)
            attName = netcdf.inqAttName(cdfid,i-1,j-1);
            attValue = netcdf.getAtt(cdfid,i-1,attName);
            if verbose
                if isnumeric(attValue)
                    disp(sprintf('%30s -> %f',attName,attValue))
                else
                    disp(sprintf('%30s -> %s',attName,attValue))
                end
            end
            out.([varName{i} '_attributes']).(attName) = attValue;
        end
    end
    return  
end

sdim = strmatch('scan_number',dimName);
if isempty(sdim)
     error(message('bioinfo:mzcdfread:scan_NumberMissing'))
end
numS = dimLength(sdim);

% extracting info for reading a subset of scans, if these variables can not be
% found the function errors
vidsi = strmatch('scan_index',varName,'exact')-1;
vidsat = strmatch('scan_acquisition_time',varName,'exact')-1;

if isempty(vidsi)
    error(message('bioinfo:mzcdfread:scan_indexMissing'))
end
I =  netcdf.getVar(cdfid,vidsi,'double');

if ~isempty(tmrange) 
    if isempty(vidsat)
        error(message('bioinfo:mzcdfread:scan_acquisition_timeMissing'))
    else
        T =  netcdf.getVar(cdfid,vidsat,'double');
    end
end

% Figure out the subset of scans to read
J = [I(2:end);numP];   % find the end index for every scan:
h = true(numS,1);      % start with all scans
if ~isempty(indices)   % filter only those that are in SCANINDICES
    h = h & ismember((1:numS)',indices);
end
if ~isempty(tmrange)   % filter only those that are in TIMERANGE
    h = h & ((T>=tmrange(1))&(T<=tmrange(2)));
end

if ~any(h)
    warning(message('bioinfo:mzcdfread:emptySubset'))
end

% read the subset of scans
if verbose
    disp('Reading netCDF variables:')
    disp('                          Name    Type  Dim_ID  NumAtt')
end
for i = 1:nvars
    if verbose
        str = sprintf('%d ',varDim{i}+1);
        disp(sprintf('%30s %6d  %6s %6d',varName{i},varType(i),str,varNatts(i)))
    end
    if any(varDim{i}+1 == pdim) % when dimension is point_number read by block
        II = I(h);
        JJ = J(h);
        hh = find(II(2:end)~=JJ(1:end-1));
        if ~isempty(hh)
          II = II([1;hh+1]);
          JJ = JJ([hh;end]);
        end
        out.(varName{i}) = zeros(sum(JJ-II),1,class(netcdf.getVar(cdfid,i-1,0,1)));
        newScanIdx = zeros(numel(II),1,class(netcdf.getVar(cdfid,vidsi,0,1)));
        r = 1;
        for l = 1:numel(II)
            newScanIdx(l) = r;
            for j = II(l):blcksize:JJ(l)-1
                k = min(j+blcksize,JJ(l));
                out.(varName{i})(r:r+k-j-1) = netcdf.getVar(cdfid,i-1,j,k-j);
                r = r + k - j;
            end
        end
    elseif any(varDim{i}+1 == sdim) % when dimension is scan_number read the elements of subsets only
        tmp = netcdf.getVar(cdfid,i-1);
        out.(varName{i}) = tmp(h);
    else % read the whole variable at once
        out.(varName{i}) = netcdf.getVar(cdfid,i-1);
    end
    if ischar(out.(varName{i})) && isvector(out.(varName{i}))
        out.(varName{i}) = strtrim(strrep(out.(varName{i})(:)',char(0),''));
    end
    for j = 1:varNatts(i)
        attName = netcdf.inqAttName(cdfid,i-1,j-1);
        attValue = netcdf.getAtt(cdfid,i-1,attName);
        if verbose
            if isnumeric(attValue)
                disp(sprintf('%30s -> %f',attName,attValue))
            else
                disp(sprintf('%30s -> %s',attName,attValue))
            end
        end
        out.([varName{i} '_attributes']).(attName) = attValue;
    end
end
% scan_index is updated because it is used to extract scans from the
% structure, zero indexing convention is preserved as the netcdf file.
out.scan_index = newScanIdx-1;
