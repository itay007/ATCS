function out = affxgenericread(filename)
%AFFXGENERICREAD reads Affymetrix CC Generic files.
%
%   genStruct = AFFXGENERICREAD(FILE) reads an Affymetrix CC Generic file
%   FILE and creates a structure genStruct.
%
%   See also AFFYREAD, AGFEREAD, CELINTENSITYREAD, GPRREAD,
%   PROBELIBRARYINFO, PROBESETLINK, PROBESETLOOKUP, PROBESETPLOT,
%   PROBESETVALUES, SPTREAD.
%
%   GeneChip and Affymetrix are registered trademarks of Affymetrix, Inc.

% Copyright 2007-2012 The MathWorks, Inc.


fid = fopen(filename,'rb','ieee-be', 'US-ASCII');
if fid == -1
    error(message('bioinfo:affxgenericread:badFile', filename));
end
fileType = fread(fid,1,'uchar'); 
% check that magic numbers
if  fileType ~= 59
    error(message('bioinfo:affxgenericread:badFile', filename));
end

paramValStruct = struct('name','','value',0,'type',0);
dataColStruct = struct('name','','type',0,'size',0);
dataSetStruct = struct('name','','nextDataSet',0,'firstDataElem',0,...
    'numParams',0,'params',[],'numCols',0,'cols',[],'numRows',0,'data',[]);
dataGroupStruct = struct('name','','nextGroupPos',0,'firstDataSet',0,...
    'numDataSets',0,'dataSets',[]);

dataTypes = {'*char','*uchar','*int16','*uint16','*int32','*uint32','*float','double','*char','*char'};

% read the file header and generic data header
out.name = filename;
out.version = fread(fid,1,'uchar'); 
out.numDataGroups = fread(fid,1,'int32'); 
out.firstGroupPos = fread(fid,1,'uint32');
out.dataType = readStr(fid); 
out.fileID = readStr(fid); 
out.dateTime = readWStr(fid); 
out.locale = readWStr(fid); 
out.numParams = fread(fid,1,'int32');
headerParams = repmat(paramValStruct,out.numParams,1);

for count = 1:out.numParams
    name = readWStr(fid);
    headerParams(count).name = name;
    theStartPos = ftell(fid);
    value = readStr(fid);
    headerParams(count).value = value;
    type = readWStr(fid);
    headerParams(count).type = type;
    theEndPos = ftell(fid);
    if ~isempty(strfind(type,'text/x-calvin'))
        type = strrep(type,'text/x-calvin-','');
        type = strrep(type,'unsigned-','u');
        % treating all ints as int32 seems to be the way to go
        type = regexprep(type,'integer-.*$','int');
        fseek(fid,theStartPos,-1);
        fread(fid,1,'int32');
        value = fread(fid,1,type);
        headerParams(count).value = value;
        fseek(fid,theEndPos,-1);
    end
end
out.params = headerParams;

% skip the parent headers
% numParentHeaders = fread(fid,1,'int32');
fseek(fid,out.firstGroupPos,-1);

dataGroup = repmat(dataGroupStruct,out.numDataGroups,1);
for groupCount = 1:out.numDataGroups
    %fscanf(fid,'%s',1)
    dataGroup(groupCount).nextGroupPos = fread(fid,1,'uint32');
    dataGroup(groupCount).firstDataSet = fread(fid,1,'uint32');
    dataGroup(groupCount).numDataSets = fread(fid,1,'int32');
    dataGroup(groupCount).name = readWStr(fid);

    dataSet = repmat(dataSetStruct,dataGroup(groupCount).numDataSets,1);
    for setCount = 1:dataGroup(groupCount).numDataSets
        dataSet(setCount).firstDataElem =  fread(fid,1,'uint32');
        dataSet(setCount).nextDataSet =  fread(fid,1,'uint32');
        dataSet(setCount).name = readWStr(fid);
        dataSet(setCount).numParams = fread(fid,1,'int32');

        if dataSet(setCount).numParams > 0
            paramVals = repmat(paramValStruct,dataSet(setCount).numParams,1);
            for paramCount = 1:dataSet(setCount).numParams
                paramVals(paramCount).name = readWStr(fid);
                paramVals(paramCount).value = readStr(fid);
                paramVals(paramCount).type = readWStr(fid);
            end
            dataSet(setCount).params = paramVals;
        end
        dataSet(setCount).numCols = fread(fid,1,'uint32');

        if dataSet(setCount).numCols > 0
            dataCol = repmat(dataColStruct,dataSet(setCount).numCols,1);
            for colCount = 1:dataSet(setCount).numCols
                dataCol(colCount).name = readWStr(fid);
                dataCol(colCount).type = fread(fid,1,'char');
                dataCol(colCount).size = fread(fid,1,'int32');
            end
            dataSet(setCount).cols = dataCol;

            dataSet(setCount).numRows = fread(fid,1,'uint32');
            % This assumes that all columns are the same size
            if range([dataCol.size]) == 0 && range([dataCol.type]) == 0

                [dataSet(setCount).data, countCheck] = fread(fid,...
                    (dataSet(setCount).numRows*dataSet(setCount).numCols),...
                    dataTypes{dataCol(1).type+1});  % only tried with uint16
                
                if countCheck ~= (dataSet(setCount).numRows*dataSet(setCount).numCols)
                    warning(message('bioinfo:affxgenericread:readWarning', filename));
                end
            else
                dataSet(setCount).data = zeros(dataSet(setCount).numRows,dataSet(setCount).numCols);
                for rowCount = 1:dataSet(setCount).numRows
                    for colCount = 1:dataSet(setCount).numCols
                        dataSet(setCount).data(rowCount,colCount) = double(fread(fid,1,dataTypes{dataCol(colCount).type+1}));
                    end
                end
            end
        end
        fStatus = fseek(fid,dataSet(setCount).nextDataSet,-1);
        if fStatus == -1
           error(message('bioinfo:affxgenericread:readError', filename));
        end
    end
    dataGroup(groupCount).dataSets = dataSet;
    if dataGroup(groupCount).nextGroupPos > 0
        fseek(fid,dataGroup(groupCount).nextGroupPos,-1);
    end
end
fclose(fid);
out.dataGroups = dataGroup;
%%%%%%%%%%%%%%%%%%%%%%

function str = readStr(fid)
strLength = fread(fid,1,'int32');
str=fread(fid,strLength,'int8=>char')';
str(str ==0) = '';
%%%%%%%%%%%%%%%%%%%%%%

function str = readWStr(fid)
strLength = fread(fid,1,'int32');
str= fread(fid,strLength,'uint16=>char')';
