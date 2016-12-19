function data = tgspcread(filename,varargin)
% TGSPCREAD Reads a Thermo-Galactic SPC format file into MATLAB as a structure.
%
%   OUT = TGSPCREAD(FILENAME) reads a Thermo Scientific SPC format file
%   into MATLAB as a structure. FILENAME is a string containing a file
%   name, or a path and a file name, of a file that conforms to the SPC
%   specification. OUT is a MATLAB structure with the following fields:
%   'Header','X','Y','Z'. The fields 'X','Y', and 'Z', store the X, Y and Z
%   data from the file. The Header field is a structure with the following
%   fields: 'Filename','FileSize', 'ExperimentType', 'NumDataPoints', 'XFirst',
%   'XLast','NumScans', 'XLabel', 'YLabel', 'ZLabel', 'CollectionTime',
%   'CollectionTimeDatenum', 'Resolution', 'SourceInstrument',
%   'InterferogramPeakPointNumber', 'Comment','CustomAxisUnitLabel',
%   'SubScanHeaders','ZValues'. 
%
%   TGSPCREAD(...,'ZRANGE', RANGE) specifies which range of Z values to
%   extract scans from the SPC file.  RANGE is a two-element numeric
%   array, [START END].  START and END must be scalar values and START
%   must be less then END.  The default is to extract all scans.  If
%   SCANINDICES option is used, then ZRANGE option can not be used.
%
%   TGSPCREAD(...,'SCANINDICES',INDICES) specifies an index, a vector of
%   indices or range of scans indices with which to extract scans. INDICES
%   is a positive integer or a vector of integers.  To indicate a range of
%   indices, use [START_IND:END_IND]. START_IND is an integer indicating
%   the scan number with which to start extracting scan information. 
%   END_IND is an integer indicating the scan number with which to stop
%   extracting scan information. START_IND must be less then END_IND.  The
%   default is to extract all scans.  If ZRANGE option is used, then
%   SCANINDICES option can not be used.
%
%   TGSPCREAD(...,'VERBOSE',T/F) show progress of reading the SPC file.
%   Default is false.
%
%   Example:
%
%       out = tgspcread('results.spc');
%       % plot the first scan in the SPC file:
%       plot(out.X,out.Y(:,1));
%
%   Note that the file results.spc is not provided.
%   See also JCAMPREAD, MZCDF2PEAKS, MZCDFINFO, MZCDFREAD, MZXML2PEAKS,
%   MZXMLINFO, MZXMLREAD, TGSPCINFO.

% Copyright 2009 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);
[zrange,scanindices,verbose,headeronly] = parse_inputs(varargin{:});

[fid, theMessage] = fopen(filename,'rb','l');
if fid == -1
    error(message('bioinfo:tgspcread:CouldNotOpenFile', theMessage));
end
dirInfo = dir(filename);
c = onCleanup(@()fclose(fid));
if(verbose)
    fprintf('Reading header for file: %s\n',filename);
end
data = [];
% BYTE ftflgs; /* Flag bits defined below */
header.ftflgs = fread(fid,1,'*uint8');
header.flags = lookupFtflgs(header.ftflgs);
% BYTE fversn; /* 0x4B=> new LSB 1st, 0x4C=> new MSB 1st, 0x4D=> old
header.fversn = fread(fid,1,'*uint8');
if header.fversn == 77 % old format
    
    %     BYTE oftflgs;
    % BYTE oversn; /* 0x4D rather than 0x4C or 0x4B */
    % short oexp; /* Word rather than byte */
    % float onpts; /* Floating number of points */
    % float ofirst; /* Floating X coordinate of first pnt (SP rather than DP) */
    % float olast; /* Floating X coordinate of last point (SP rather than DP) */
    % BYTE oxtype; /* Type of X units */
    % BYTE oytype; /* Type of Y units */
    % WORD oyear; /* Year collected (0=no date/time) - MSB 4 bits are Z type */
    % BYTE omonth; /* Month collected (1=Jan) */
    % BYTE oday; /* Day of month (1=1st) */
    % BYTE ohour; /* Hour of day (13=1PM) */
    % BYTE ominute; /* Minute of hour */
    % char ores[8]; /* Resolution text (null terminated unless 8 bytes used) */
    % WORD opeakpt;
    % WORD onscans;
    % float ospare[7];
    % char ocmnt[130];
    % char ocatxt[30];
    % char osubh1[32]; /* Header for first (or main) subfile included in main header */
    
    
    header.fexp = fread(fid,1,'*uint16');
    header.fexperStr = lookupFexper(0);
    header.fnpts = fread(fid,1,'single');
    header.ffirst = fread(fid,1,'*single');
    header.flast = fread(fid,1,'*single');
    header.fxtype = fread(fid,1,'*uint8');
    header.fxtypeStr = lookupFxzwtype(header.fxtype);
    header.fytype = fread(fid,1,'*uint8');
    header.fytypeStr = lookupFytype(header.fytype);
    header.fyear = fread(fid,1,'*uint16');
    header.fztype = bitshift(header.fyear,-12);
    header.fztypeStr = lookupFxzwtype(header.fztype);
    header.fyear = double(bitshift(bitshift(header.fyear,4),-4));
    if header.fyear < 100
        header.fyear = header.fyear + 1900;
    end
    header.fmonth = fread(fid,1,'uint8');
    header.day= fread(fid,1,'uint8');
    header.hour = fread(fid,1,'uint8');
    header.fminute = fread(fid,1,'uint8');
    
    try
        header.fdatenum = datenum(header.fyear,header.fmonth,header.day,header.hour,header.fminute,0);
    catch
        header.fdatenum = 0;
    end
    header.fdatestr = datestr(header.fdatenum);
    header.fres= fread(fid,8,'uint8=>char')';
    header.fpeakpt = fread(fid,1,'*uint16');
    header.fnscans = fread(fid,1,'*uint16');
    header.fspare = fread(fid,7,'*single');
    header.fsource = '';
    header.fcmnt= fread(fid,130,'uint8=>char')';
    header.fcatxt= fread(fid,30,'uint8=>char')';
    header.fnsub = 1;
    if(verbose)
        fprintf('File contains %d scans\n',header.fnsub);
    end
    %header.osubh1= fread(fid,32,'uint8=>char')';
    % subheader
    % BYTE subflgs; /* Flags as defined above */
    % char subexp; /* Exponent for sub-file's Y values (80h=>float) */
    % WORD subindx; /* Integer index number of trace subfile (0=first) */
    % float subtime; /* Floating time for trace (Z axis coordinate) */
    % float subnext; /* Floating time for next trace (May be same as beg) */
    % float subnois; /* Floating peak pick noise level if high byte nonzero */
    % DWORD subnpts; /* Integer number of subfile points for TXYXYS type */
    % DWORD subscan; /* Integer number of co-added scans or 0 (for collect) */
    % float subwlevel; /* Floating W axis value (if fwplanes non-zero) */
    % char subresv[4]; /* Reserved area (must be set to zero) */
    
    subheader.subflgs= fread(fid,1,'*uint8');
    subheader.subexp= fread(fid,1,'uint8=>char');
    subheader.subindx= fread(fid,1,'*uint16');
    subheader.subtime= fread(fid,1,'*single');
    subheader.subnext= fread(fid,1,'*single');
    subheader.subnois= fread(fid,1,'*single');
    subheader.subnpts= fread(fid,1,'*uint32');
    subheader.subscan= fread(fid,1,'*uint32');
    subheader.subwlevel= fread(fid,1,'*single');
    subheader.subresv= fread(fid,4,'uint8=>char')';
    
    if isnan(scanindices)
        reportVal = true;
        scanindices = 1;
    end
    valWords = fread(fid,header.fnpts,'*uint32');
    % need to switch the words
    valWords = bitor(bitshift(valWords,16),bitshift(valWords,-16));
    signBit = bitget(valWords,32);
    vals = 2^double(header.fexp)*(double(valWords) - 2^32 * double(signBit))/2^32;
    %vals = 2^double(header.oexp)*(fread(fid,header.onpts,'int32'))/2^32;
    xdata = linspace(header.ffirst,header.flast,header.fnpts)';
    
    % header */
else % new format
    % need to check endianness
    if header.fversn == 76
        
        fclose(fid);
        fid = fopen(filename,'rb','b');
        header.ftflgs = fread(fid,1,'*uint8');
        header.fversn = fread(fid,1,'*uint8');
    end
    if header.fversn <75 || header.fversn>76
         error(message('bioinfo:tgspcread:BadFormat', filename));
    end
    % new format
    
    
    % BYTE fexper; /* Instrument technique code (see below) */
    header.fexper = fread(fid,1,'*uint8');
    header.fexperStr = lookupFexper(header.fexper);
    % char fexp; /* Fraction scaling exponent integer (80h=>float) */
    header.fexp = double(fread(fid,1,'uint8=>char'));
    % DWORD fnpts; /* Integer number of points (or TXYXYS directory position) */
    header.fnpts = fread(fid,1,'uint32');
    % double ffirst; /* Floating X coordinate of first point */
    header.ffirst = fread(fid,1,'*double');
    % double flast; /* Floating X coordinate of last point */
    header.flast = fread(fid,1,'*double');
    % DWORD fnsub; /* Integer number of subfiles (1 if not TMULTI) */
    header.fnsub = fread(fid,1,'*uint32');
    % BYTE fxtype; /* Type of X axis units (see definitions below) */
    header.fxtype = fread(fid,1,'*uint8');
    header.fxtypeStr = lookupFxzwtype(header.fxtype);
    % BYTE fytype; /* Type of Y axis units (see definitions below) */
    header.fytype = fread(fid,1,'*uint8');
    header.fytypeStr = lookupFytype(header.fytype);
    % BYTE fztype; /* Type of Z axis units (see definitions below) */
    header.fztype = fread(fid,1,'*uint8');
    header.fztypeStr = lookupFxzwtype(header.fztype);
    % BYTE fpost; /* Posting disposition (see GRAMSDDE.H) */
    header.fpost = fread(fid,1,'*uint8');
    % DWORD fdate; /* Date/Time LSB: min=6b,hour=5b,day=5b,month=4b,year=12b */
    header.fdate = fread(fid,1,'*uint32');
    header.fdatenum = unpackFdate(header.fdate);
    header.fdatestr = datestr(header.fdatenum);
    % char fres[9]; /* Resolution description text (null terminated) */
    header.fres = deblank(fread(fid,9,'uint8=>char')');
    % char fsource[9]; /* Source instrument description text (null terminated) */
    header.fsource = deblank(fread(fid,9,'uint8=>char')');
    % WORD fpeakpt; /* Peak point number for interferograms (0=not known) */
    header.fpeakpt = fread(fid,1,'*uint16');
    % float fspare[8]; /* Used for Array Basic storage */
    header.fspare = fread(fid,8,'*single');
    % char fcmnt[130]; /* Null terminated comment ASCII text string */
    header.fcmnt = deblank(fread(fid,130,'uint8=>char')');
    % char fcatxt[30]; /* X,Y,Z axis label strings if ftflgs=TALABS */
    header.fcatxt = fread(fid,30,'uint8=>char')';
    % DWORD flogoff; /* File offset to log block or 0 (see above) */
    header.flogoff = fread(fid,1,'*uint32');
    % DWORD fmods; /* File Modification Flags (see below: 1=A,2=B,4=C,8=D..) */
    header.fmods = fread(fid,1,'*uint32');
    % BYTE fprocs; /* Processing code (see GRAMSDDE.H) */
    header.fprocs = fread(fid,1,'*uint8');
    % BYTE flevel; /* Calibration level plus one (1 = not calibration data) */
    header.flevel = fread(fid,1,'*uint8');
    % WORD fsampin; /* Sub-method sample injection number (1 = first or only ) */
    header.fsampin = fread(fid,1,'*uint16');
    % float ffactor; /* Floating data multiplier concentration factor (IEEE-32) */
    header.ffactor = fread(fid,1,'*single');
    % char fmethod[48]; /* Method/program/data filename w/extensions comma list */
    header.fmethod = deblank(fread(fid,48,'uint8=>char')');
    % float fzinc; /* Z subfile increment (0 = use 1st subnext-subfirst) */
    header.fzinc = fread(fid,1,'*single');
    % DWORD fwplanes; /* Number of planes for 4D with W dimension (0=normal) */
    header.fwplanes = fread(fid,1,'*uint32');
    % float fwinc; /* W plane increment (only if fwplanes is not 0) */
    header.fwinc = fread(fid,1,'*single');
    % BYTE fwtype; /* Type of W axis units (see definitions below) */
    header.fwtype = fread(fid,1,'*uint8');
    header.fwtypeStr = lookupFxzwtype(header.fwtype);
    % char freserv[187]; /* Reserved (must be set to zero) */
    header.freserv = deblank(fread(fid,187,'uint8=>char')');
    
    % we may need to read xdata if header.flags.xyxys
    if header.flags.xvals &&  ~header.flags.xyxys
        % Read xdata
        xdata = fread(fid,header.fnpts,'*single');
    elseif header.flags.xyxys
        xdata = cell(1,header.fnsub);
    else
        xdata = linspace(header.ffirst,header.flast,header.fnpts)';
    end
    subheader(header.fnsub).subflgs = 0;
    vals = cell(1,header.fnsub);
    reportVal = true(1,header.fnsub);
    if(verbose)
        fprintf('File contains %d scans\n',header.fnsub);
    end
    if isnan(scanindices)
        scanindices = 1:header.fnsub;
    end
    for traceCount = 1:header.fnsub
        % subheader
        
        % BYTE subflgs; /* Flags as defined above */
        subheader(traceCount).subflgs = fread(fid,1,'*uint8');
        % char subexp; /* Exponent for sub-file's Y values (80h=>float) */
        subheader(traceCount).subexp = fread(fid,1,'uint8=>char')';
        % WORD subindx; /* Integer index number of trace subfile (0=first) */
        subheader(traceCount).subindx = fread(fid,1,'*uint16');
        % float subtime; /* Floating time for trace (Z axis coordinate) */
        subheader(traceCount).subtime = fread(fid,1,'*single');
        % float subnext; /* Floating time for next trace (May be same as beg) */
        subheader(traceCount).subnext = fread(fid,1,'*single');
        % float subnois; /* Floating peak pick noise level if high byte nonzero */
        subheader(traceCount).subnois = fread(fid,1,'*single');
        % DWORD subnpts; /* Integer number of subfile points for TXYXYS type */
        subheader(traceCount).subnpts = fread(fid,1,'*uint32');
        % DWORD subscan; /* Integer number of co-added scans or 0 (for collect) */
        subheader(traceCount).subscan = fread(fid,1,'*uint32');
        % float subwlevel; /* Floating W axis value (if fwplanes non-zero) */
        subheader(traceCount).subwlevel = fread(fid,1,'*single');
        % char subresv[4]; /* Reserved area (must be set to zero) */
        subheader(traceCount).subresv = deblank(fread(fid,4,'uint8=>char')');
        
        % now read the data
        % we may need to read xdata if header.flags.xyxys
        if header.flags.xyxys
            numPoints = double(subheader(traceCount).subnpts);
            % Read xdata
            xdata{traceCount} = fread(fid,numPoints,'*single');
        else
            numPoints = header.fnpts;
        end
        % fexp == 0x80 when we are using IEEE floats
        if header.fexp == 128
            vals{traceCount} = fread(fid,numPoints,'*single');
        else
            if header.flags.sprec
                valWords = fread(fid,numPoints,'*uint16');
                % need to switch the words
                signBit = bitget(valWords,16);
                vals{traceCount} = 2^double(header.fexp)*(double(valWords) - 2^32 * double(signBit))/2^32;
            else
                valWords = fread(fid,numPoints,'*uint32');
                % need to switch the words
                signBit = bitget(valWords,32);
                vals{traceCount} = 2^double(header.fexp)*(double(valWords) - 2^32 * double(signBit))/2^32;
            end
        end
        if(headeronly) || subheader(traceCount).subtime < zrange(1) || subheader(traceCount).subtime > zrange(2) || ~ismember(traceCount,scanindices)
            vals(traceCount) = [];
            reportVal(traceCount) = false;
            if header.flags.xyxys
                xdata{traceCount} = [];
            end
        end
    end
end
outputHeader.Filename = filename;
outputHeader.FileSize = dirInfo.bytes;
outputHeader.ExperimentType = header.fexperStr;
outputHeader.NumDataPoints = header.fnpts;
outputHeader.XFirst = header.ffirst;
outputHeader.XLast = header.flast;
outputHeader.NumScans = header.fnsub;
outputHeader.XLabel = header.fxtypeStr;
outputHeader.YLabel = header.fytypeStr;
outputHeader.ZLabel = header.fztypeStr;
outputHeader.CollectionTime = header.fdatestr;
outputHeader.CollectionTimeDatenum = header.fdatenum;
outputHeader.Resolution = header.fres;
outputHeader.SourceInstrument = strtrim(deblank(header.fsource));
outputHeader.InterferogramPeakPointNumber = header.fpeakpt;
outputHeader.Comment = strtrim(deblank(header.fcmnt));
outputHeader.CustomAxisUnitLabel = strtrim(deblank(header.fcatxt));


data.Header = outputHeader;
data.Header.SubScanHeaders = subheader;
allZVals = [subheader.subtime];
data.Header.ZValues = allZVals;
if headeronly
    data = data.Header;
else
    if iscell(xdata) && numel(xdata) == 1
        xdata = xdata{1};
    end
    data.X = xdata;
    if iscell(vals) && ~iscell(xdata)
        data.Y = cell2mat(vals(reportVal));
    else
        data.Y = vals(:,reportVal);
    end
    data.Z = allZVals(reportVal);
end


function str = lookupFexper(fexper)
% Instrument Technique fexper settings */
% In older software, the TCGRAM in ftflgs must be set if fexper is non-zero. */
% A general chromatogram is specified by a zero fexper when TCGRAM is set. */
switch(fexper)
    case 0
        str = 'General SPC';
        % format has this as 
        % str = 'General SPC (could be anything)';
    case 1
        str = 'Gas Chromatogram';
    case 2
        str = 'General Chromatogram';
        % str = 'General Chromatogram (same as SPCGEN with TCGRAM)';
    case 3
        str = 'HPLC Chromatogram';
    case 4
        str = 'FT-IR, FT-NIR, FT-Raman Spectrum or Igram';
        %str = 'FT-IR, FT-NIR, FT-Raman Spectrum or Igram (Can also be used for scanning IR.)';
    case 5
        str = 'NIR Spectrum';      
        % str = 'NIR Spectrum (Usually multi-spectral data sets for calibration.)';
    case 7
        str = 'UV-VIS Spectrum';
        % str = 'UV-VIS Spectrum (Can be used for single scanning UV-VIS-NIR)';
    case 8
        str = 'X-ray Diffraction Spectrum';
    case 9
        str = 'Mass Spectrum';
        % str = 'Mass Spectrum (Can be single, GC-MS, Continuum, Centroid or TOF.)';
    case 10
        str = 'NMR Spectrum or FID';
    case 11
        str = 'Raman Spectrum';
        % str = 'Raman Spectrum (Usually Diode Array, CCD, etc. use SPCFTIR for FT-Raman.)';
    case 12
        str = 'Fluorescence Spectrum';
    case 13
        str = 'Atomic Spectrum';
    case 14
        str = 'Chromatography Diode Array Spectra';
end

function str = lookupFxzwtype(val)
% Possible settings for fxtype, fztype, fwtype. */
% XEV and XDIODE - XMETERS are new and not supported by all software. */
switch val
    case 0
        str = 'Arbitrary';
    case 1
        str = 'Wavenumber (cm-1)';
    case 2
        str = 'Micrometers (um)';
    case 3
        str = 'Nanometers (nm)';
    case 4
        str = 'Seconds';
    case 5
        str = 'Minutes';
    case 6
        str = 'Hertz (Hz)';
    case 7
        str = 'Kilohertz (KHz)';
    case 8
        str = 'Megahertz (MHz)';
    case 9
        str = 'Mass (M/z)';
    case 10
        str = 'Parts per million (PPM)';
    case 11
        str = 'Days';
    case 12
        str = 'Years';
    case 13
        str = 'Raman Shift (cm-1)';
    case 14
        str = 'eV';
    case 15
        str = 'XYZ text labels in fcatxt';
    case 16
        str = 'Diode Number';
    case 17
        str = 'Channel';
    case 18
        str = 'Degrees';
    case 19
        str = 'Temperature (F)';
    case 20
        str = 'Temperature (C)';
    case 21
        str = 'Temperature (K)';
    case 22
        str = 'Data Points';
    case 23
        str = 'Milliseconds (mSec)';
    case 24
        str = 'Microseconds (uSec)';
    case 25
        str = 'Nanoseconds (nSec)';
    case 26
        str = 'Gigahertz (GHz)';
    case 27
        str = 'Centimeters (cm)';
    case 28
        str = 'Meters (m)';
    case 29
        str = 'Millimeters (mm)';
    case 30
        str = 'Hours';
    case 255
        str = 'Double interferogram';
end
function [str, mbPositive] = lookupFytype(val)
% Possible settings for fytype. (The first 127 have positive peaks.) */
% YINTENS - YDEGRK and YEMISN are new and not supported by all software. */
switch val
    case 0,
        str = 'Arbitrary Intensity';
    case 1
        str = 'Interferogram';
    case 2
        str = 'Absorbance';
    case 3
        str = 'Kubelka-Monk';
    case 4
        str = 'Counts';
    case 5
        str = 'Volts';
    case 6
        str = 'Degrees';
    case 7
        str = 'Milliamps';
    case 8
        str = 'Millimeters';
    case 9
        str = 'Millivolts';
    case 10
        str = 'Log(1/R)';
    case 11
        str = 'Percent';
    case 12
        str = 'Intensity';
    case 13
        str = 'Relative Intensity';
    case 14
        str = 'Energy';
    case 16
        str = 'Decibel';
    case 19
        str = 'Temperature (F)';
    case 20
        str = 'Temperature (C)';
    case 21
        str = 'Temperature (K)';
    case 22
        str = 'Index of Refraction [N]';
    case 23
        str = 'Extinction Coeff. [K]';
    case 24
        str = 'Real';
    case 25
        str = 'Imaginary';
    case 26
        str = 'Complex';
    case 128
        str = 'Transmission (ALL HIGHER MUST HAVE VALLEYS!)';
    case 129
        str = 'Reflectance';
    case 130
        str = 'Arbitrary or Single Beam with Valley Peaks';
    case 131
        str = 'Emission';
end
mbPositive = val<=127;

function flags = lookupFtflgs(val)
% Possible bit FTFLGS flag byte settings. */
% Note that TRANDM and TORDRD are mutually exclusive. */
% Code depends on TXVALS being the sign bit. TXYXYS must be 0 if TXVALS=0. */
% In old software without the fexper code, TCGRAM specifies a chromatogram. */

% TSPREC 1 /* Single precision (16 bit) Y data if set. */
flags.sprec = logical(bitand(val,uint8(1)));
% TCGRAM 2 /* Enables fexper in older software (CGM if fexper=0) */
flags.cgram = logical(bitand(val,uint8(2)));
% TMULTI 4 /* Multiple traces format (set if more than one subfile) */
flags.multi = logical(bitand(val,uint8(4)));
% TRANDM 8 /* If TMULTI and TRANDM=1 then arbitrary time (Z) values */
flags.randm = logical(bitand(val,uint8(8)));
% TORDRD 16 /* If TMULTI abd TORDRD=1 then ordered but uneven subtimes */
flags.ordrd = logical(bitand(val,uint8(16)));
% TALABS 32 /* Set if should use fcatxt axis labels, not fxtype etc. */
flags.alabs = logical(bitand(val,uint8(32)));
% TXYXYS 64 /* If TXVALS and multifile, then each subfile has own X's */
flags.xyxys = logical(bitand(val,uint8(64)));
% TXVALS 128 /* Floating X value array preceeds Y's (New format only) */
flags.xvals = logical(bitand(val,uint8(128)));


% /* FMODS spectral modifications flag setting conventions: */
% /* "A" (2^01) = Averaging (from multiple source traces) */
% /* "B" (2^02) = Baseline correction or offset functions */
% /* "C" (2^03) = Interferogram to spectrum Computation */
% /* "D" (2^04) = Derivative (or integrate) functions */
% /* "E" (2^06) = Resolution Enhancement functions (such as deconvolution) */
% /* "I" (2^09) = Interpolation functions */
% /* "N" (2^14) = Noise reduction smoothing */
% /* "O" (2^15) = Other functions (add, subtract, noise, etc.) */
% /* "S" (2^19) = Spectral Subtraction */
% /* "T" (2^20) = Truncation (only a portion of original X axis remains) */
% /* "W" (2^23) = When collected (date and time information) has been modified */
% /* "X" (2^24) = X units conversions or X shifting */
% /* "Y" (2^25) = Y units conversions (transmission->absorbance, etc.) */
% /* "Z" (2^26) = Zap functions (features removed or modified) */


function dNum = unpackFdate(fdate)
% DWORD fdate; /* Date/Time LSB: min=6b,hour=5b,day=5b,month=4b,year=12b */
[year, fdate] =  getPackedBits(fdate,12);
[month, fdate] = getPackedBits(fdate,4);
[day, fdate] = getPackedBits(fdate,5);
[hour, fdate] = getPackedBits(fdate,5);
minutes = getPackedBits(fdate,6);
dNum = datenum(year, month, day, hour, minutes, 0);

function [val, theRem] = getPackedBits(theInt,numBits)
val = double(bitshift(theInt,numBits-32));
theRem = bitshift(theInt,numBits);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zrange,scanindices,verbose,headeronly] = parse_inputs(varargin)
% Parse the varargin parameter/value inputs

% Check that we have the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:tgspcread:IncorrectNumberOfArguments', mfilename));
end

% The allowed inputs
okargs = {'zrange','scanindices','verbose','headeronly'};

% Set default values
verbose = true;
zrange = [];
scanindices = NaN;
headeronly = false;

% Loop over the values
for j=1:2:nargin
    % Lookup the pair
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        
        case 1  % zrange
            if ~isnumeric(pval) || ~isvector(pval) || numel(pval)<2||numel(pval)>2
                error(message('bioinfo:tgspcread:ZRangeNonNumeric'))
            elseif any(scanindices>0)
                error(message('bioinfo:tgspcread:ConflictingPVPairsLevelsOrIndices'))
                
            end
            zrange = pval;
            
        case 2  % scanindices
            if ~isnumeric(pval) || ~isvector(pval) || any(pval<0)||numel(pval)<2
                error(message('bioinfo:tgspcread:IndicesNonNumeric'))
            elseif  ~isempty(zrange)
                error(message('bioinfo:tgspcread:ConflictingPVPairsLevelsOrTimeRange'))
            end
            scanindices = pval;
        case 3  % verbose
            verbose = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 4  % headeronly
            headeronly = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end
    
end
if isempty(zrange)
    zrange = [-inf inf];
end
