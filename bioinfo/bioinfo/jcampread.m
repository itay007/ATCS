function output = jcampread(filename,warnFlag)
%JCAMPREAD reads JCAMP-DX (DX) files
%
%   JCAMPDATA = JCAMPREAD(FILE) reads in JCAMP-DX format data from FILE
%   and creates a structure JCAMPDATA, containing these fields:
%           Title
%           DataType
%           Origin
%           Owner
%           Blocks
%           Notes
%
%   The Blocks field of the structure is an array of structures
%   corresponding to each set of data in the file. These structures have
%   the following fields:
%           XData
%           YData
%           XUnits
%           YUnits
%           Notes
%
%   The function supports Versions 4.24 and 5 of the JCAMP-DX
%   format for infrared and mass spectrometry.
%
%   For more details on the JCAMP-DX format, see
%   http://www.jcamp-dx.org/index.html.
%
%   Sample data is available from http://www.jcamp-dx.org/testdata.html.
%
%   Example:
%
%       % Read in a sample JCAMP-DX file and plot the frequency spectrum.
%       jcampStruct = jcampread('brukaffn.dx')
%       data = jcampStruct.Blocks(1);
%       plot(data.XData,data.YData);
%       title(jcampStruct.Title);
%       xlabel(data.XUnits);
%       ylabel(data.YUnits);
%
%   See also MSLOWESS, MSSGOLAY, MSVIEWER, MZXMLREAD.

%  Note: http://www.jcamp-dx.org
%  Copyright 1997-2007  International Union of Pure and Applied Chemistry

% Copyright 2003-2012 The MathWorks, Inc.



% Check first if the input is a file identifier (used for the compound
% structure) or a filename string.
openedFile = 0;
if nargin < 2
    warnFlag = true;
end
try
    pos = ftell(filename);
catch allExceptions %#ok<NASGU>
    pos = -1;
end
if (pos < 0)
    try
        fid = fopen(filename,'rt');
        openedFile = 1;
    catch allExceptions %#ok<NASGU>
        fid = -1;
    end
    if fid == -1
        error(message('bioinfo:jcampread:CannotOpenJCAMPFile', filename));
    end
else
    fid = filename;
end
output = [];

% Read in and verify the first 5/6 lines to make sure the format is correct
headerField = {'##TITLE','##JCAMP-DX','##DATA TYPE','##DATA CLASS','##ORIGIN','##OWNER','##NUM DIM'};
headerField2 = {'##TITLE','##JCAMPDX','##DATATYPE','##DATACLASS','##ORIGIN','##OWNER','##NUMDIM'};
allHeaders = union(headerField,headerField2);
headers = cell(size(headerField));
count = 1;
while count <= 6
    % Retrieve the Labeled Data Record
    pos = ftell(fid);
    line=getnext(fid);
    [field,rem]=strtok(line,'='); 
    field = strtrim(field);
    if strcmp(field,'##')||~strncmp(field,'##',2)
        continue;
    elseif ~ismember(field,allHeaders)
        if warnFlag
            warning(message('bioinfo:jcampread:BadJCAMPHeader', headerField{ count }, field))
            warnFlag = false;
        end
        fseek(fid,pos,-1);
        break;
    end
    [found position] = ismember(field,headerField);
    if ~found
        [~, position] = ismember(field,headerField2);
    end
    % Remove whitespaces and comments from  the value
    value = strtok(rem,'=');
    headers{position} = strtrim(strtok(value,'$'));

    % Skip the Data Class field if the version is less than 5
    % use str2double to deal with [].
    if (position==3) && (str2double(headers{2})<5)
        count=4;
    end
    count = count+1;
end

% Check if using a supported version of JCAMP file format
if ~isempty(headers{2})
    version = sscanf(headers{2},'%f');
    if (version < 5)&&(version ~= 4.24)
        warning(message('bioinfo:jcampread:BadJCAMPVersion', sprintf( '%0.2f', version )));
    end
end

numBlocks = 0;
DataSpec = struct('XFactor',1,'YFactor',1,'XUnits','','YUnits','',...
    'NPoints',-1,'FirstX',-1,'LastX',-1,'FirstY',-1);
nCounter = 1; % Keeps track of the notes index
notes = cell(1,2);
inNote = false;
while ~feof(fid)
    pos = ftell(fid);
    line = getnext(fid);

    % Make sure it is a vald line and extract field/value pairs
    try
        [field,value]=strread(line,'##%s %s','delimiter','=');
        inNote = false;
    catch allExceptions
        % some example files have multi-line notes without the ## delimiter
        if inNote
            notes{nCounter,1}=line;
            nCounter = nCounter+1;
        end
        continue;
    end
    try
        field = upper(strtrim(field{1}));
    catch allExceptions
        continue;
    end
    try
        value = strtrim(value{1});
    catch allExceptions %#ok<NASGU>
        value = '';
    end

    if strcmp(field,'XUNITS')
        DataSpec.XUnits = value;
    elseif strcmp(field,'YUNITS')
        DataSpec.YUnits = value;
	elseif strcmp(field,'ORIGIN')
        if isempty(headers{5})
            headers{5} = value;
        end
    elseif strcmp(field,'OWNER')
        if isempty(headers{6})
            headers{6} = value;
        end
    elseif strcmp(field,'XFACTOR')
        DataSpec.XFactor = sscanf(value,'%f');
    elseif strcmp(field,'YFACTOR')
        DataSpec.YFactor = sscanf(value,'%f');
    elseif strcmp(field,'FIRSTX')
        DataSpec.FirstX = sscanf(value,'%f');
    elseif strcmp(field,'LASTX')
        DataSpec.LastX = sscanf(value,'%f');
    elseif strcmp(field,'NPOINTS')
        DataSpec.NPoints = sscanf(value,'%d');
    elseif strcmp(field,'FIRSTY')
        DataSpec.FirstY = sscanf(value,'%f');
    elseif strcmp(field,'BLOCKS')
        numBlocks = sscanf(value,'%f');
    elseif strcmp(field,'TITLE')
        % Run subroutine for compound data
        fseek(fid, pos, 'bof'); % Move back one line to start at TITLE
        dataBlocks = compoundRead(fid,numBlocks, warnFlag);
        break;
    elseif strcmp(field,'XYDATA')
        % Run subroutine for xydata
        dataBlocks(1) = xyRead(fid,DataSpec);
        break;
    elseif strcmp(field,'NTUPLES')
        % Run subroutine for ntuples
        numBlocks = 1; %#ok Multiple dataBlocks, so expect an additional END field
        dataBlocks = ntupleRead(fid);
        break;
    elseif strcmp(field,'PEAK TABLE')
        % Run subroutine for peak table
        dataBlocks(1) = peakRead(fid,DataSpec);
        break;
    elseif strcmp(field,'END')
        % If you find an END without a multiple dataBlock structure, there
        % is an error in the file format.
        if numBlocks < 1
            error(message('bioinfo:jcampread:BadJCAMPDataType', headers{ 3 }))
        end
    else
        % Store field/value as a note
        notes{nCounter,1}=field;
        notes{nCounter,2}=value;
        nCounter = nCounter+1;
         % some example files have multi-line notes without the ##
         % delimiter
        inNote = true;
    end
end

if(openedFile)
    fclose(fid);
end

if ~exist('dataBlocks','var')
    error(message('bioinfo:jcampread:BadJCAMPDataType', headers{ 3 }))
end

output.Title = headers{1};
output.DataType = headers{3};
if ~isempty(headers{4})
    output.DataClass = headers{4};
end
output.Origin = headers{5};
output.Owner = headers{6};
output.Blocks = dataBlocks;
output.Notes = notes;

% -----------------------------------------------------------------
% PEAKREAD
% Handles uncompressed data in the (XY..XY) format
% -----------------------------------------------------------------
function dataBlock = peakRead(fid,DataSpec)
if (DataSpec.NPoints>0)
    x = zeros(1,DataSpec.NPoints);
    y = zeros(1,DataSpec.NPoints);
else
    x = [];
    y = [];
end

j = 1;
data = getnext(fid);
% some files have delimiter ',' so we need to be robust to that.
formatString = '%f, %f';
delimiter = ';';
try % try the default format of commas and semi-colons.
    [x0,y0] = strread(data,formatString,'delimiter',delimiter); %#ok
catch allExceptions %#ok<NASGU>
    try % try using comma as delimiter for everything
        formatString = '%f %f';
        delimiter = ',';
        [x0,y0] = strread(data,formatString,'delimiter',delimiter); %#ok
    catch anotherException%#ok<NASGU> % go back to default
        warning(message('bioinfo:jcampread:BadDelimiter'));
        formatString = '%f, %f';
        delimiter = ';';
    end
end
while ~strncmp(data,'##',2)
    [xi yi] = strread(data,formatString,'delimiter',delimiter);
    k = j+length(xi)-1;
    x(j:k) = xi;
    y(j:k) = yi;
    j = k+1;
    data = getnext(fid);
end

fseek(fid,-(length(data)+getNewLineLength(fid)+1),'cof');

dataBlock.XData = x*DataSpec.XFactor;
dataBlock.YData = y*DataSpec.YFactor;
dataBlock.XUnits = DataSpec.XUnits;
dataBlock.YUnits = DataSpec.YUnits;

% -----------------------------------------------------------------
% XYREAD
% Handles data in the (X++(Y..Y)) form, including SQZ,DIF, and DUP
% compressed data.
% -----------------------------------------------------------------
function dataBlock = xyRead(fid,DataSpec)
tol = 1; % Set the tolerance for the checkpoints
xTolWarnFlag = true;
yTolWarnFlag = true;

SQZ = 'ihgfedcba@ABCDEFGHI';
DIF = 'rqponmlkj%JKLMNOPQR';
DUP = 'STUVWXYZs';

if (DataSpec.NPoints>0)
    x = linspace(DataSpec.FirstX/DataSpec.XFactor,...
        DataSpec.LastX/DataSpec.XFactor,DataSpec.NPoints);
    y = zeros(1,DataSpec.NPoints);
else
    error(message('bioinfo:jcampread:BadJCAMPDataPoints'));
end

j=1;
mode = 1;

data = getnext(fid);
while ~strncmp(data,'##',2)
  values = regexp(data,'[+-]*[a-s_A-Z_%@]?\d*\.?\d*([Ee][+-]\d+)?','match');
    xval = sscanf(values{1},'%f');
    if abs(xval-round(x(j)))>tol
        if xTolWarnFlag
            warning(message('bioinfo:jcampread:XOutofTol'))
            xTolWarnFlag = false;
        end
    end
    for i = 2:length(values)
        num = values{i};
        if ~isempty(strfind(SQZ,num(1)))
            dig = find(SQZ==num(1))-10;
            num = sscanf([sprintf('%d',dig) num(2:end)],'%f');
            if (i==2)&&(j>1)&&(mode==2)
                if abs(num-y(j))>tol
                    if yTolWarnFlag
                        warning(message('bioinfo:jcampread:YOutofTol'))
                        yTolWarnFlag = false;
                    end
                    y(j) = num;
                end
            else
                y(j) = num;
            end
            j = j+1;
        elseif ~isempty(strfind(DIF,num(1)))
            difDig = find(DIF==num(1))-10;
            difValue = sscanf([sprintf('%d',difDig) num(2:end)],'%f');
            y(j) = y(j-1)+difValue;
            j = j+1;
            mode = 2;
        elseif ~isempty(strfind(DUP,num(1)))
            dupDig = find(DUP==num(1));
            if length(num)==1
                dupValue = dupDig;
            else
                dupValue = sscanf([sprintf('%d',dupDig) num(2:end)],'%f');
            end
            for k = 1:dupValue-1
                if mode==1
                    y(j) = y(j-1);
                else
                    y(j) = y(j-1)+difValue;
                end
                j = j+1;
            end
        else
            if (i==2)&&(j>1)&&(mode==2)
                if abs(sscanf(num,'%f')-y(j))>tol
                    warning(message('bioinfo:jcampread:YOutofTol'))
                end
            else
                y(j) = sscanf(num,'%f');
            end
            j = j+1;
        end
    end
    if (mode==2)
        j = j-1; % rewind one value for the checkpoint
    end
    data=getnext(fid);
end
fseek(fid,-(length(data)+getNewLineLength(fid)+1),'cof');

dataBlock.XData = x*DataSpec.XFactor;
dataBlock.YData = y*DataSpec.YFactor;
dataBlock.XUnits = DataSpec.XUnits;
dataBlock.YUnits = DataSpec.YUnits;

% -----------------------------------------------------------------
% COMPOUNDREAD
% This subfunction recursively call the main function to read in each
% subblock of data.
% -----------------------------------------------------------------
function dataBlocks = compoundRead(fid,numBlocks,warnFlag)
dataBlocks = [];
if nargin < 3
    warnFlag = true;
end
for i = 1:numBlocks
    try
        DataStr = jcampread(fid,warnFlag);
        try
            dataBlocks(i) = DataStr.Blocks(1); %#ok<*AGROW>
        catch allExceptions %#ok<NASGU>
            dataBlocks(i).XData = DataStr.Blocks(1).XData;
            dataBlocks(i).YData = DataStr.Blocks(1).YData;
            dataBlocks(i).XUnits = DataStr.Blocks(1).XUnits;
            dataBlocks(i).YUnits = DataStr.Blocks(1).YUnits;
        end
        dataBlocks(i).Title = DataStr.Title;
        dataBlocks(i).DataType = DataStr.DataType;
        dataBlocks(i).Owner = DataStr.Owner;
        dataBlocks(i).Origin = DataStr.Origin;
        dataBlocks(i).Notes = DataStr.Notes;
    catch theException
        msgstr = theException.message;
        if warnFlag
            warning(message('bioinfo:jcampread:CompoundStructureError', msgstr));
            warnFlag = false;
        end
        while (~feof(fid))&&(~strncmp(fgetl(fid),'##END',5))
        end
    end
    if ~exist('dataBlocks','var')
        dataBlocks = struct;
    end
end

% -----------------------------------------------------------------
% NTUPLEREAD
% Handles Spectral Series Data
% -----------------------------------------------------------------
function dataBlocks = ntupleRead(fid)
blockNum = 1;
DataSpec = struct('XFactor',1,'YFactor',1,'XUnits','','YUnits','',...
    'NPoints',-1,'FirstX',-1,'LastX',-1,'FirstY',-1);
data = getnext(fid);
while ~strncmp(data,'##END',5)
    try
        t=regexp(data,'##([\s\w-/]*)\s*=\s*(.*)','tokens');
        field = t{1}{1};
        value = t{1}{2};
    catch allExceptions
        data = getnext(fid);
        continue;
    end

    if strncmp(field,'$',1)
        continue; % Ignore custom fields designated by ##$...
    elseif strcmp(field,'VAR_NAME')
        name = strread(value,'%s','delimiter',',');
    elseif strcmp(field,'SYMBOL')
        cellsymbol = strread(value,'%s','delimiter',',');
    elseif strcmp(field,'VAR_DIM')
        vardim = strread(value,'%d','delimiter',',');
    elseif strcmp(field,'UNITS')
        units = strread(value,'%s','delimiter',',');
    elseif strcmp(field,'FIRST')
        first = strread(value,'%f','delimiter',',');
    elseif strcmp(field,'LAST')
        last = strread(value,'%f','delimiter',',');
    elseif strcmp(field,'FACTOR')
        fctrs = strread(value,'%f','delimiter',',');
    elseif strcmp(field,'PAGE')
        [zvar zval] = strread(value,'%s %f','delimiter','=');
        zvar = strtrim(zvar{1}); % Extract one character
    elseif strcmp(field,'NPOINTS')
        DataSpec.NPoints = sscanf(value,'%d');
    elseif strcmp(field,'DATA TABLE')

        % Determine the type of data
        form = strread(value,'%s','delimiter',',');
        mode = 0;
        vars = regexp(form{1},'\((\w)(\w)..\1\2\)','tokens');
        if isempty(vars)
            mode = 1;
            vars = regexp(form{1},'\((\w.*)\+\+\((\w)..\2\)','tokens');
            if isempty(vars)
                error(message('bioinfo:jcampread:BadNtupleFormat', form{ 1 }));
            end
        end
        if ~isempty(cellsymbol)
            Xidx = find(strcmp(cellsymbol,vars{1}{1}));
            Yidx = find(strcmp(cellsymbol,vars{1}{2}));
            Zidx = find(strcmp(cellsymbol,zvar));
        else
            error(message('bioinfo:jcampread:BadJCAMPntuple'));
        end
        if exist('fctrs','var')&&~isempty(fctrs)
            DataSpec.XFactor = fctrs(Xidx);
            DataSpec.YFactor = fctrs(Yidx);
        end
        if exist('units','var')&&~isempty(units)
            DataSpec.XUnits = units{Xidx};
            DataSpec.YUnits = units{Yidx};
        end
        if exist('first','var')&&~isempty(units)
            DataSpec.FirstX = first(Xidx);
            DataSpec.FirstY = first(Yidx);
        end
        if exist('last','var')&&~isempty(last)
            DataSpec.LastX = last(Xidx);
        end
        if DataSpec.NPoints < 0
            try
                DataSpec.NPoints = vardim(Xidx);
            catch allExceptions
                error(message('bioinfo:jcampread:BadJCAMPDataPoints'));
            end
        end

        % Run subroutine based on table format
        if mode
            dataBlock = xyRead(fid,DataSpec);
        else
            dataBlock = peakRead(fid,DataSpec);
        end

        % Restructure result and add extra values
        if exist('name','var')&&~isempty(name)
            dataBlock.XName = name{Xidx};
            dataBlock.YName = name{Yidx};
            dataBlock.ZName = name{Zidx};
        end

        try
            dataBlock.ZUnits = units{Zidx};
        catch allExceptions %#ok<NASGU>
            dataBlock.ZUnits = '';
        end
        dataBlock.ZData = zval;
        dataBlocks(blockNum) = dataBlock; 
        blockNum = blockNum+1;

        % Move back up to the end of the table so next line is Page or END
%         fseek(fid,endTable,'bof');
    end
    data = getnext(fid);
end

% -----------------------------------------------------------------
% GETNEXT - gets the next valid line from the file,
% ignores comments and custom fields.
% -----------------------------------------------------------------
function data = getnext(fid)
data = fgetl(fid);
while strncmp(data,'$$',2)||strncmp(data,'##$',3)||isempty(data)
    data = fgetl(fid);
end
data = strtrim(data);
if data<0
    error(message('bioinfo:jcampread:NoJCAMPEnd'))
else
    if ~isempty(strfind(data,'$'))
        data = strtrim(strread(data,'%s',1,'delimiter','$'));
        data = data{1};
    end
end
% -----------------------------------------------------------------
% GETNEWLINELENGTH - gets the length of a newline in this file.
% -----------------------------------------------------------------
function newLineLength = getNewLineLength(fid)

currentPos = ftell(fid);
fgetlLength = length(fgetl(fid));
fseek(fid,currentPos,'bof');
fgetsLength = length(fgets(fid));
newLineLength = fgetsLength-fgetlLength;
