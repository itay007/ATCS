function out = affxexpread(fileName)
%AFFXEXPREAD reads Affymetrix EXP files.
%
%   ExpStruct = AFFXEXPREAD(FILE) reads an Affymetrix experiment file FILE
%   and creates a structure ExpStruct.
%
%   See also AFFYREAD, AGFEREAD, CELINTENSITYREAD, GPRREAD,
%   PROBELIBRARYINFO, PROBESETLINK, PROBESETLOOKUP, PROBESETPLOT,
%   PROBESETVALUES, SPTREAD.
%
%   GeneChip and Affymetrix are registered trademarks of Affymetrix, Inc.

% Copyright 2007 The MathWorks, Inc.


% get the parts of the file

[thePath,theName,theExt] = fileparts(fileName);

% set up the output structure
out.Name = [theName,theExt];
if isempty(thePath)
    fullFileName = which(out.Name);
    thePath = fileparts(fullFileName);
end

out.DataPath = thePath;
% Fill this in later once we know the chip type
out.LibPath = '';
out.FullPathName = [thePath,filesep,theName,theExt];

% we need the date stamp of the file.
dirInfo = dir(fileName);

% open the file
fid = fopen(fileName,'r');
if fid == -1
    error(message('bioinfo:affxexpread:badExpFile', fileName));
end
% Extract the lines -- note that we use bufsize as some sample files had
% very long empty lines that were longer than the default buffer size of
% textscan.
allLines = textscan(fid,'%s','delimiter','\n');
fclose(fid);

% pull the cell array apart
lines = allLines{1};
if numel(lines) >=2
    header = lines{1};
    headerVer = lines{2};
else
    header = '';
    headerVer = '';
end
% Check that the file header is really an EXP file. Warn if not.
if ~strcmpi(strtrim(header),'Affymetrix GeneChip Experiment Information') ||...
        ~strncmpi(strtrim(headerVer),'Version',7)
    warning(message('bioinfo:affxexpread:BadExpHeader'));
end

% now parse out the fields of interest

out.ChipType = getValFromTag(lines,'Chip Type');
out.LibPath = fileparts(which([out.ChipType '.CDF']));
if isfield(dirInfo,'datenum')
    out.Date = datestr(dirInfo.datenum,'mmm dd yyyy HH:MM:SS');
else
    out.Data = dirInfo.date;
end
out.ChipLot = getValFromTag(lines,'Chip Lot');
out.Operator = getValFromTag(lines,'Operator');
out.SampleType = getValFromTag(lines,'Sample Type');
out.SampleDesc = getValFromTag(lines,'Description');
out.Project = getValFromTag(lines,'Project');
out.Comments = getValFromTag(lines,'Comments');
out.Reagents = getValFromTag(lines,'Solution Type');
out.ReagentLot = getValFromTag(lines,'Solution Lot');
% The protocol details are stored between these two lines
[out.Protocol, protocolLine] = getValFromTag(lines,'Protocol');
[out.Station, stationLine] = getValFromTag(lines,'Station');
out.Module = getValFromTag(lines,'Module');
out.HybridizeDate = getValFromTag(lines,'Hybridize Date');
out.ScanPixelSize = getValFromTag(lines,'Pixel Size');
out.ScanFilter = getValFromTag(lines,'Filter');
out.ScanDate = getValFromTag(lines,'Scan Date');
out.ScannerID = getValFromTag(lines,'Scanner ID');
out.NumberOfScans = str2double(getValFromTag(lines,'Number of Scans'));
out.ScannerType = getValFromTag(lines,'Scanner Type');

% pull out the protocol steps
numProtocolSteps = stationLine-protocolLine-1;
if numProtocolSteps < 0
    warning(message('bioinfo:expread:UnexpectedProtocolOrder'));
end
out.NumProtocolSteps = numProtocolSteps;
out.ProtocolSteps = lines(protocolLine+1:stationLine-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value, lineNum] = getValFromTag(lines,tag)
% helper function that extracts the value from the file given a tag

lineNum = find(strncmpi(lines,tag,numel(tag)));
if isempty(lineNum)
    value = '';
    return
end
if ~isscalar(lineNum)
    warning(message('bioinfo:expread:badExpHeader'));
    lineNum = lineNum(1);
end
value = strtrim(lines{lineNum}(numel(tag)+1:end));

