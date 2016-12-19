function output = galread(filename)
%GALREAD reads GenePix Array List format data files.
%
%   GALDATA = GALREAD(FILE) reads in a GenePix Array List format
%   from FILE and creates a structure GALDATA, containing these fields:
%       Header 
%       BlockData
%       IDs
%       Names
%
%   The BlockData field is an Nx3 array. The columns of this array are the
%   Block data, the Column data, and the Row data, respectively.
%
%   For more information on the GAL format, see
%   http://www.moleculardevices.com/pages/software/gn_genepix_file_formats.html 
%
%   Examples:
%
%       % Copy the sample GAL file from
%       % http://www.moleculardevices.com/pages/software/gn_genepix_file_formats.html#example 
%       % to a local file called Demo.gal.
%
%       % Then use GALREAD to load the data into a MATLAB structure.
%       galread('Demo.gal')
%
%   See also AFFYREAD, AGFEREAD, GEOSOFTREAD, GPRREAD, IMAGENEREAD, SPTREAD.
%
%   GenePix is a registered trademark of Molecular Devices Corporation 

% Copyright 2002-2006 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename); 
try
    fopenMessage = '';
    [fid, fopenMessage] = fopen(filename,'rt');
catch theException %#ok<NASGU>
    fid = -1;
end

if fid == -1
    error(message('bioinfo:galread:CannotOpenGALFile', filename, fopenMessage));
end

% first line should be ATF versions
checkHeader = fgetl(fid);

if ~strncmp(checkHeader,'ATF',3)
    error(message('bioinfo:galread:BadGALFile'))
end

fileSize = sscanf(fgetl(fid),'%d');

for count = 1:fileSize(1)
    line = strrep(fgetl(fid),'"','');
    [field,val] = strtok(line,'='); 
    val = deblank(val(2:end));
    v = str2num(val); %#ok
    if ~isempty(v)
        blockNum = sscanf(field,'Block%d');
        if isempty(blockNum)
            header.(field) = v;
        else
            header.Block(blockNum,:) = v;
        end
        
    else
        header.(field) = val;
    end
end

% now deal with the data
colNames = strread(fgetl(fid),'%s','delimiter','\t');
colNames = strrep(colNames,'"','');

% suck the data into memory
currPos = ftell(fid);
fseek(fid,0,1);
endPos = ftell(fid);
fseek(fid,currPos,-1);

data = strread(fread(fid,endPos-currPos,'uchar=>char'),'%s','delimiter','\t','emptyvalue',NaN);

fclose(fid);

% pull this apart

numRows = length(data)/fileSize(2);
if floor(numRows) ~= numRows
    error(message('bioinfo:galread:ProblemsReadingGAL', filename));
end
data = reshape(data,fileSize(2),numRows)';

% find the column ordering -- the format says that this is flexible
blockCol = find(strcmpi(colNames,'Block'));
columnCol = find(strcmpi(colNames,'Column'));
rowCol = find(strcmpi(colNames,'Row'));
idCol = find(strcmpi(colNames,'ID'));
nameCol = find(strcmpi(colNames,'Name'));
udCol = find(strncmpi(colNames,'User',4));

% let str2num deal with the unknown shape
% if this is slow then change it to sscanf
blocks = str2num(char(data(:,blockCol))); %#ok
columns = str2num(char(data(:,columnCol))); %#ok
rows = str2num(char(data(:,rowCol))); %#ok

blockdata = [blocks,columns,rows];

% deal with empty vals
if ~isempty(idCol)
    IDs = data(:,idCol);
else 
    IDs = [];
end
if ~isempty(nameCol)
    names = data(:,nameCol);
else
    names = [];
end

if ~isempty(udCol)
    userDefined = data(:,udCol);
else
    userDefined = [];
end

% create the output structure
output.Header = header;
output.BlockData = blockdata;
output.IDs =IDs;
output.Names = names;
if ~isempty(userDefined)
    output.UserDefined = userDefined;
end
