function info = tgspcinfo(filename)
% TGSPCINFO Reports information about an SPC format file.
%   INFO = TGSPCINFO(FILENAME) returns a structure whose fields contain
%   information about an SPC format file.  FILENAME is a string containing a file
%   name, or a path and a file name, of an SPC format file.
%   INFO is a structure with the following fields:
%
%   Filename - The name of the file.
%
%   FileSize - The size of the file in bytes.
%
%	ExperimentType - The instrumental technique that created the original
%    data. 
% 
%	NumDataPoints - The number of data points (Y values)
% 
%	XFirst - The first X value in the file.
% 
%   XLast -  The last X value in the file.
% 
%   NumScans - The number of scans or sub files in the file.
% 
%   XLabel - The label for X values.
% 
%   YLabel - The label for Y values.
% 
%   ZLabel -The label for Z values.
% 
%   CollectionTime - Date and time information for when the scans were
%   collected.
% 
%   CollectionTimeDatenum - CollectionTime as a MATLAB datenum (see HELP
%   DATENUM).
% 
%   Resolution - The instrument resolution.
%
%   SourceInstrument - The source instrument
% 
%   InterferogramPeakPointNumber -  TODO
% 
%   Comment - A free text comment.
% 
%   CustomAxisUnitLabel - TODO
% 
%   SubScanHeaders - Header information for sub scans.
% 
%   ZValues - A vector containing the Z values of the scans in the file.
%
%   Example:
%
%   % Get information about a file called sample.spc. 
%
%          info = tgspcinfo('sample.spc');
%
%   Note that the file sample.spc is not provided. 
%
%   See also JCAMPREAD, MZCDF2PEAKS, MZCDFINFO, MZCDFREAD, MZXML2PEAKS,
%   MZXMLINFO, MZXMLREAD, TGSPCREAD.

% Copyright 2009 The MathWorks, Inc.



bioinfochecknargin(nargin,1,mfilename);

info = tgspcread(filename,'headeronly',true);
