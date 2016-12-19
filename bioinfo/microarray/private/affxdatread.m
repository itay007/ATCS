function out = affxdatread(filename)
%AFFXDATREAD reads Affymetrix DAT files.
%
%   DatStruct = AFFXDATREAD(DATFILE) reads an Affymetrix image file DATFILE
%   and creates a structure DatStruct.
%
%   See also AFFYREAD, AGFEREAD, CELINTENSITYREAD, GPRREAD,
%   PROBELIBRARYINFO, PROBESETLINK, PROBESETLOOKUP, PROBESETPLOT,
%   PROBESETVALUES, SPTREAD.
%
%   GeneChip and Affymetrix are registered trademarks of Affymetrix, Inc.

% Copyright 2007 The MathWorks, Inc.


fid = fopen(filename,'rb','ieee-le');
if fid == -1
    error(message('bioinfo:affxdatread:badDatFile', filename));
end
% Read the header.
% Some information is not returned. We could change this if necessary but
% the old reader did not return this.

fileType = fread(fid,1,'uchar'); %#ok
fclose(fid);
% check that magic numbers
if fileType == hex2dec('FC')
    out = readVersion1(filename);
elseif fileType == 59
    out = readCommandConsoleFormat(filename);
else
    error(message('bioinfo:affxdatread:invalidDatFile', filename))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = readVersion1(filename)
fid = fopen(filename,'rb','ieee-le');
fileType = fread(fid,1,'uchar'); %#ok
try
    NumPixelsPerRow = fread(fid,1,'uint16');
    NumRows = fread(fid,1,'uint16');
    totalPixels = fread(fid,1,'uint32');%#ok
    MinData = fread(fid,1,'uint32');
    MaxData = fread(fid,1,'uint32');
    meanPixel = fread(fid,1,'double');%#ok
    stdPixel = fread(fid,1,'double');%#ok
    numPerRowChar = fread(fid,9,'*char')';%#ok
    numRowsChar = fread(fid,9,'*char')';%#ok
    pixelWidthChar = fread(fid,7,'*char')';
    pixelHeightChar = fread(fid,7,'*char')';%#ok
    scanSpeedChar = fread(fid,6,'*char')';
    temperatureChar = fread(fid,7,'*char')';%#ok
    laserPowerChar = fread(fid,4,'*char')';%#ok
    scanDate = strtrim(fread(fid,18,'*char')');
    chipData = fread(fid,220,'*char')';
    dcOffset = fread(fid,1,'double');%#ok
    stdDCOffset = fread(fid,1,'double');%#ok
    dcOffsetSamples = fread(fid,1,'uint32');%#ok
    upperLeft = fread(fid,2,'uint16');
    upperRight = fread(fid,2,'uint16');
    lowerRight = fread(fid,2,'uint16');
    lowerLeft = fread(fid,2,'uint16');

    cellMargin = fread(fid,1,'uint16');
    experimentName = fread(fid,154,'*char')';%#ok

    % read the image
    imageData = fread(fid,NumPixelsPerRow*NumRows,'*uint16');
catch
    fclose(fid);
    error(message('bioinfo:affxdatread:invalidDatFile', filename))
end
fclose(fid);

% get the full path information
[thePath,thefile,theExt] = fileparts(filename);
if isempty(thePath)
    thePath = pwd;
end
dirInfo = dir(filename);


chipDataFields = strread(chipData,'%s','delimiter',char(hex2dec('14')));

% populate the output structure

out.Name = [thefile,theExt];
out.DataPath = thePath;
out.LibPath = thePath;
out.FullPathName = [thePath,filesep,thefile,theExt];
out.ChipType = strtok(chipDataFields{3},'.');
if isfield(dirInfo,'datenum')
    datestr(dirInfo.datenum,'mmm dd yyyy HH:MM:SS');
else
    out.Data = dirInfo.date;
end
out.NumPixelsPerRow = NumPixelsPerRow;
out.NumRows = NumRows;
out.MinData = MinData;
out.MaxData = MaxData;
out.PixelSize = sscanf(pixelWidthChar,'XIN=%d');
out.CellMargin = cellMargin;
out.ScanSpeed = sscanf(scanSpeedChar,'VE=%d');
out.ScanDate = datestr(datenum(scanDate,'mm/dd/yy HH:MM:SS'),...
    'dd-mmm-yyyy HH:MM:SS');
out.ScannerID = chipDataFields{1};
out.UpperLeftX = upperLeft(1);
out.UpperLeftY = upperLeft(2);
out.UpperRightX = upperRight(1);
out.UpperRightY = upperRight(2);
out.LowerLeftX = lowerLeft(1);
out.LowerLeftY = lowerLeft(2);
out.LowerRightX = lowerRight(1);
out.LowerRightY = lowerRight(2);
out.ServerName = '';
out.Image = reshape(imageData,NumPixelsPerRow,NumRows)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = readCommandConsoleFormat(filename)

genStruct = affxgenericread(filename);

% get the full path information
[thePath,thefile,theExt] = fileparts(filename);
if isempty(thePath)
    thePath = pwd;
end

% populate the output structure

out.Name = [thefile,theExt];
out.DataPath = thePath;
out.LibPath = thePath;
out.FullPathName = [thePath,thefile,theExt];

out.ChipType = genStruct.params(strcmp('affymetrix-array-type',{genStruct.params.name})).value;

out.Date =  genStruct.params(strcmp('affymetrix-scan-date',{genStruct.params.name})).value;


out.NumPixelsPerRow = genStruct.params(strcmp('affymetrix-pixel-cols',{genStruct.params.name})).value;
out.NumRows = genStruct.params(strcmp('affymetrix-pixel-rows',{genStruct.params.name})).value;
out.MinData = genStruct.params(strcmp('affymetrix-min-pixel-intensity',{genStruct.params.name})).value;
out.MaxData = genStruct.params(strcmp('affymetrix-max-pixel-intensity',{genStruct.params.name})).value;
out.PixelSize = genStruct.params(strcmp('affymetrix-pixel-size' ,{genStruct.params.name})).value;
out.CellMargin = 0;
out.ScanSpeed = 0;
out.ScanDate = genStruct.params(strcmp('affymetrix-scan-date' ,{genStruct.params.name})).value;
out.ScannerID = genStruct.params(strcmp('affymetrix-scanner-id',{genStruct.params.name})).value;

pixelDataNdx = strcmpi({genStruct.dataGroups.dataSets.name},'Pixel');
gridDataNdx = find(strcmpi({genStruct.dataGroups.dataSets.name},'GlobalGrid'));

colNames = {genStruct.dataGroups.dataSets(gridDataNdx).cols.name};
theCol = strcmpi(colNames,'Upper left x');
out.UpperLeftX = genStruct.dataGroups.dataSets(gridDataNdx).data(theCol);
theCol = strcmpi(colNames,'Upper left y');
out.UpperLeftY = genStruct.dataGroups.dataSets(gridDataNdx).data(theCol);
theCol = strcmpi(colNames,'Upper right x');
out.UpperRightX = genStruct.dataGroups.dataSets(gridDataNdx).data(theCol);
theCol = strcmpi(colNames,'Upper right y');
out.UpperRightY = genStruct.dataGroups.dataSets(gridDataNdx).data(theCol);
theCol = strcmpi(colNames,'Lower left x');
out.LowerLeftX = genStruct.dataGroups.dataSets(gridDataNdx).data(theCol);
theCol = strcmpi(colNames,'Lower left y');
out.LowerLeftY = genStruct.dataGroups.dataSets(gridDataNdx).data(theCol);
theCol = strcmpi(colNames,'Lower right x');
out.LowerRightX = genStruct.dataGroups.dataSets(gridDataNdx).data(theCol);
theCol = strcmpi(colNames,'Lower right y');
out.LowerRightY = genStruct.dataGroups.dataSets(gridDataNdx).data(theCol);

out.ServerName = '';
out.Image = reshape(genStruct.dataGroups.dataSets(pixelDataNdx).data,out.NumPixelsPerRow,out.NumRows)';
