function outputStruct=geoseriesread(geotext)
% GEOSERIESREAD reads in Gene Expression Omnibus GSE format data.
%
%   GEODATA = GEOSERIESREAD(FILE) reads in Gene Expression Omnibus (GEO)
%   GSE format Series data from FILE and creates a structure GEODATA,
%   containing the following fields:
%       Header
%       Data
%
%   Where Data is a DataMatrix Object.   
%
%   FILE can also be a URL or a MATLAB character array that contains the
%   text of a GEO Series format file.
%
%   Example:
%
%       % Get a file from GEO and save it to a file.
%       geodata = getgeodata('GSE11287','ToFile','GSE11287.txt')
%
%       % In subsequent MATLAB sessions you can use geoseriesread to access
%       % the local copy from disk instead of accessing it from the GEO web 
%       % site.
%       geodata = geoseriesread('GSE11287.txt')
%
%       % You can access the sample IDs using the colnames property of the
%       % DataMatrix Object.
%       sampleIDs = geodata.Data.colnames
%
%   See also AFFYREAD, AGFEREAD, GALREAD, GEOSOFTREAD, GETGEODATA, GPRREAD,
%   SPTREAD.

% Copyright 2008-2012 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);

if ~ischar(geotext) && ~iscellstr(geotext)
    error(message('bioinfo:geoseriesread:InvalidStringInput'));
end

% If the input is a string of GenBank data then first character should be ^
isAstring = ischar(geotext) && ~isempty(geotext) && geotext(1) == '^';

if ~isAstring
    if ~iscellstr(geotext)
        if (strfind(geotext, '://'))
            if (~usejava('jvm'))
                error(message('bioinfo:geoseriesread:NoJava'))
            end
            % must be a URL
            geotext = urlread(geotext);
            % clean up "&amp;"
            geotext = strrep(geotext, '&amp;', '&');
        elseif exist(geotext,'file')
            % it is a file
            geotext = fopen(geotext, 'rt');
            % Close file on function exit
            cleanup = onCleanup(@() fclose(geotext));
        else
            error(message('bioinfo:geoseriesread:FileNotExist', geotext))
        end
    else
        error(message('bioinfo:geoseriesread:InvalidSOFTFile'))
    end
end

geotext = textscan(geotext,'%s','delimiter','\n','Whitespace','','ReturnOnError',false);
geotext = geotext{1};

emptyLines = cellfun('isempty',geotext);
geotext(emptyLines) = [];
numLines = size(geotext,1);
numSamples = 0;
% format as follows
% ^ start line with ID
% ! comment lines
% # column descriptions
% data values


% outputStruct.Header
outputStruct.Header.Series = [];
outputStruct.Header.Samples = [];
sampleIDs = {};

lineCount = 1;
while lineCount < numLines
    theLine = geotext{lineCount};
    if strncmpi(theLine,'!Series',7)
        theField = strtrim(regexp(theLine,'Series_([^"]*)"','tokens','once'));
        theValue = regexp(theLine,'"([^"]*)"','tokens','once');
        try
            if~isempty(theField)
                if  ~isfield(outputStruct.Header.Series,(genvarname(theField{1})))
                    outputStruct.Header.Series.(genvarname(theField{1})) = theValue{1};
                else
                    outputStruct.Header.Series.(genvarname(theField{1})) = sprintf('%s\n%s',outputStruct.Header.Series.(genvarname(theField{1})),theValue{1});
                end
            end
        catch allExceptions %#ok<NASGU>
        end
    elseif strncmpi(theLine,'!Sample',7)
        theField = strtrim(regexp(theLine,'Sample_([^"]*)"','tokens','once'));
        theValue = regexp(theLine,'"([^"]*)"','tokens');
        for cellCount = 1:numel(theValue)
            if isscalar(theValue{cellCount})
                theValue{cellCount} = char(theValue{cellCount});
            end
        end
        try
            if~isempty(theField)
                if  ~isfield(outputStruct.Header.Samples,(genvarname(theField{1})))
                    outputStruct.Header.Samples.(genvarname(theField{1})) = theValue;
                else
                    outputStruct.Header.Samples.(genvarname(theField{1})) =...
                        [outputStruct.Header.Samples.(genvarname(theField{1}));theValue];
                end
            end        
        catch allExceptions %#ok<NASGU>
        end
    elseif strncmpi(theLine,'"ID_REF"',7)
        theValue = regexp(theLine,'"([^"]*)"','tokens');
        
        sampleIDs = [theValue{2:end}];
        numSamples = numel(sampleIDs);
        break
    elseif strncmp(theLine,'^SAMPLE=',8) || strncmp(theLine,'^database',9)
          error(message('bioinfo:geoseriesread:SoftFile'))
    end
    lineCount = lineCount+1;
end

fullData = sprintf('%s\n',geotext{lineCount+1:end-1});
% handle case with the string 'null'
if ~isempty(regexpi(fullData,'null','once'))
    fullData = strrep(fullData,'null','NaN ');
    fullData = strrep(fullData,'NULL','NaN ');
    fullData = strrep(fullData,'Null','NaN ');
end
if ~isempty(regexp(geotext{lineCount+1},'"','once'))
    firstValFormat = '%q';
else
    firstValFormat = '%f';
end
formatString = [firstValFormat repmat('%f',1,numSamples)];
allData = textscan(fullData,formatString,'delimiter','\t','ReturnOnError',false);
outputStruct.Data = bioma.data.DataMatrix(cell2mat(allData(2:end)),'Rownames',allData{1},'Colnames',sampleIDs);
