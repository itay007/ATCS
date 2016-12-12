function output = ilmnbslookup(filename,theID,varargin)
%ILMNBSLOOKUP looks up Illumina BeadStudio Probe Sequence and Annotation information.
%
%   ANNOTATION = ILMNBSLOOKUP(ANNOTATIONFILE,ID) looks up Illumina
%   BeadStudio probe sequence and annotation information for target ID from
%   annotation file, ANNOTATIONFILE and creates a structure ANNOTATION,
%   containing the information from the annotation file.
%
%   ID can be a single ID or a cell array of IDs.
%
%   Annotation files, for example HumanRef-8_V3_0_R0_11282963_A.txt, can be
%   downloaded from gene expression - probe sequences and annotation area
%   of the download site on www.illumina.com. Sample data files are also
%   available for download from www.ilumina.com.
%
%   ILMNBSLOOKUP(...,'LOOKUPFIELD',LOOKUP) looks for ID in field LOOKUP.
%   The default lookup field is Search_Key.
%
%   Example:
%
%       % Read in a sample BeadStudio file.
%       % Note that this file is not provided with Bioinformatics Toolbox.
%       ilmnStruct = ilmnbsread('TumorAdjacent-probe-raw.txt')
%
%       % Lookup annotation information for the 10th entry in annotation
%       % file 'HumanRef-8_V3_0_R0_11282963_A.txt'.
%       srchCol = find(strcmpi('Search_Key',ilmnStruct.TextColumnNames));
%
%       annotation = ilmnbslookup('HumanRef-8_V3_0_R0_11282963_A.txt',...
%                                        ilmnStruct.TextData{10,srchCol})
%
%   See also ILMNBSREAD.
%
%   Illumina is a registered trademark of Illumina, Inc.

% Copyright 2007-2010 The MathWorks, Inc.

%   Undocummented:
%   ILMNBSLOOKUP(...,'DELIMITER',delimiter) lets you set the delimiter for
%   the file. Default is '\t'.
%
%   Previously undocumented but no longer needed (R2013a):
%   ILMNBSLOOKUP(...,'BUFFERSIZE',bufSize) lets you set the buffersize used
%   by textscan. Default is 4095.

% Handle inputs
bioinfochecknargin(nargin, 2, mfilename)
[lookupField,delimiter] = parse_inputs(varargin{:});

% These are the names that we expect to see in the files.
expectedCols = {'Search_Key','Species','Source','Transcript','ILMN_Gene',...
    'Source_Reference_ID','RefSeq_ID','Unigene_ID','Entrez_Gene_ID',...
    'GI','Accession','Symbol','Protein_Product','Probe_Id','Array_Address_Id',...
    'Probe_Type','Probe_Start','Probe_Sequence','Chromosome','Probe_Chr_Orientation',...
    'Probe_Coordinates','Cytoband','Definition','Ontology_Component',...
    'Ontology_Process','Ontology_Function','Synonyms','Obsolete_Probe_Id'};
expectedControlCols = {'Probe_Id','Array_Address_Id','Probe_Sequence',...
    'Reporter_Group_Name','Reporter_Group_id','Reporter_Composite_map'};

% These are the expected names for miRNA files.
expectedmiRNACols = {'miRNA','ILMN_Gene','Search_Key','Probe_Sequence',...
    'Array_Address_Id','Illumicode_Seq','Chromosome','Ploidy','Species',...
    'Probe_Coordinates','Mature_miRNA_seq','Probe_Chr_Orientation'};
expectedmiRNAControlCols = {'Probe_Id','Array_Address_Id',...
    'Reporter_Group_Name','Reporter_Group_id','Color'};


% Handle old CSV format files
[~,~,ext] = fileparts(filename);
if strcmpi(ext,'.csv')
    output = ilmnbslookupcsv(filename,theID,lookupField);
    return
end
cleanUpBGX = false;
if strcmpi(ext,'.bgx')
    try
        cleanUpBGX = true;
        filename = gunzip(filename,tempname);
    catch mExcp
        if strcmp(mExcp.identifier,'MATLAB:checkfilename:invalidFilename')
            error(message('bioinfo:ilmnbslookup:BadBGXFile', filename));
        else
            rethrow(mExcp)
        end
        
    end
    filename = filename{1};
end

try
    fopenMessage = '';
    [fid, fopenMessage] = fopen(filename,'rt');
catch theException %#ok<NASGU>
    fid = -1;
end

if fid == -1
    bgxCleanup(filename,cleanUpBGX);
    error(message('bioinfo:ilmnbslookup:CannotOpenFile', filename, fopenMessage));
end


% Read the header information

headerlines = 0;
maxHeaderSize = 30;
numProbes = Inf;
numControls = 0;
while headerlines < maxHeaderSize
    line = fgetl(fid);
    headerlines = headerlines + 1;
    probeVal = sscanf(lower(line),'number of probes %d');
    controlVal = sscanf(lower(line),'number of controls %d');
    miRna = strcmpi(line,sprintf('Assay Format\tmiRNA'));
    if miRna
        expectedCols = expectedmiRNACols;
        expectedControlCols = expectedmiRNAControlCols;
    end
    if ~isempty(probeVal)
        numProbes = probeVal;
    end
    if ~isempty(controlVal)
        numControls = controlVal;
    end
    if strcmpi(line,'[Probes]')
        break
    end
end
if headerlines >= maxHeaderSize
    fclose(fid);
    bgxCleanup(filename,cleanUpBGX);
    error(message('bioinfo:ilmnbslookup:UnknownFormat', filename));
end
% Read the column names
colNames = strread(fgetl(fid),'%s','delimiter',delimiter);
% set up some variables
numCols = numel(colNames);
formatStr = repmat('%q',1,numCols);
numTries = 0;
maxTries = 10;
controlData = [];
while numTries < maxTries;
    if isinf(numProbes)
        probeData = textscan(fid,formatStr,'delimiter',delimiter);
    else
        probeData = textscan(fid,formatStr,numProbes,'delimiter',delimiter);
    end
    if numControls > 0
        line = fgetl(fid);
        if ~strcmpi(line,'[Controls]')
            line = fgetl(fid);
            if ~strcmpi(line,'[Controls]')
                numControls = 0;
                break
            end
        end
        controlColNames = strread(fgetl(fid),'%s','delimiter',delimiter);
        numControlCols = numel(controlColNames);
        
        controlFormatStr = repmat('%q',1,numControlCols);
        controlData = textscan(fid,controlFormatStr,numControls,'delimiter',delimiter);
    end
    break
end

fclose(fid);
bgxCleanup(filename,cleanUpBGX);
% Create output struct
allCols = union(expectedCols,expectedControlCols);
output = cell2struct(repmat({''},1,numel(allCols)),allCols,2);
if numControls > 0
    output = doLookup(output,controlData,theID,expectedControlCols,controlColNames,lookupField,filename,true);
end
output = doLookup(output,probeData,theID,expectedCols,colNames,lookupField,filename,false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bgxCleanup(filename,tf)
% cleanup the bgx temp file
if tf
    if iscell(filename)
        filename = char(filename);
    end
    pathstr = fileparts(filename);
    delete(filename);
    rmdir(pathstr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = doLookup(output,theData,theID,expectedCols,colNames,lookupField,filename,controlFlag)
% Do the lookup based for either probe or control information

lookupCol = find(strcmpi(lookupField,expectedCols));
if isempty(lookupCol)
    if ~controlFlag
        error(message('bioinfo:ilmnbslookup:LookupFieldNotValid', lookupField));
    else
        lookupCol = 1;
    end
end
[colFound, colMap] = ismember(lower(expectedCols),lower(strrep(colNames,' ','_')));

% warn if don't find the expected columns
missingFields = find(~colFound);
if ~isempty(missingFields)
    for count = missingFields
        warning(message('bioinfo:ilmnbslookup:ExpectedFieldNotFound', expectedCols{ count }, filename));
    end
end

% now loop over the IDs

if ischar(theID)
    theID = {theID};
end
% Make an array if necessary
numIDs = numel(theID);
if numel(output) ~= numIDs
    output = repmat(output,numIDs,1);
end

for IDCount = 1:numel(theID)
    thisID = theID{IDCount};
    row = find(strcmpi(thisID,[theData{:,colMap(lookupCol)}]));
    
    if isempty(row)
        if ~controlFlag && isempty(output(IDCount).Probe_Id)
            warning(message('bioinfo:ilmnbslookup:IDNotFound', thisID, filename));
        end
        continue
    end
    
    if numel(row)>1
        warning(message('bioinfo:ilmnbslookup:IDNotUnique', thisID, filename));
    end
    
    % Put data into structure for output
    for count = 1:max(colMap)
        if colFound(count)
            % if multiple matches are found return as a cell array.
            % Otherwise return as strings.
            if numel(row) > 1
                data = {theData{colMap(count)}{row}};
            else
                data = theData{colMap(count)}{row};
            end
            output(IDCount).(expectedCols{count}) = data;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = ilmnbslookupcsv(filename,theID,lookupField)
% Handle old CSV format files

% try to open the file
try
    fopenMessage = '';
    [fid, fopenMessage] = fopen(filename,'rt');
catch theException %#ok<NASGU>
    fid = -1;
end

if fid == -1
    error(message('bioinfo:ilmnbslookup:CannotOpenFile', filename, fopenMessage));
end

% Read the first line and extract the column names

colNames = strread(fgetl(fid),'%s','delimiter',',');

% set up some variables
numCols = numel(colNames);

formatStr = repmat('%q',1,numCols);

allData = textscan(fid,formatStr,'delimiter',',');

fclose(fid);

output = struct('Search_key','',...
    'Target','',...
    'ProbeId','',...
    'Gid','',...
    'Transcript','',...
    'Accession','',...
    'Symbol','',...
    'Type','',...
    'Start','',...
    'Probe_Sequence','',...
    'Definition','',...
    'Ontology','',...
    'Synonym','');

expectedCols = fieldnames(output);
lookupCol = find(strcmpi(lookupField,expectedCols));
if isempty(lookupCol)
    error(message('bioinfo:ilmnbslookup:LookupFieldNotValid', lookupField));
end
[colFound, colMap] = ismember(expectedCols,colNames);

% warn if don't find the expected columns
missingFields = find(~colFound);
if ~isempty(missingFields)
    for count = missingFields
        warning(message('bioinfo:ilmnbslookup:ExpectedFieldNotFound', expectedCols{ count }, filename));
    end
end

% now loop over the IDs

if ischar(theID)
    theID = {theID};
end
% Make an array if necessary
numIDs = numel(theID);
output(numIDs) = output;

for IDCount = 1:numel(theID)
    thisID = theID{IDCount};
    row = find(strcmpi(thisID,[allData{:,colMap(lookupCol)}]));
    
    if isempty(row)
        warning(message('bioinfo:ilmnbslookup:IDNotFound', thisID, filename));
        return
    end
    
    if numel(row)>1
        warning(message('bioinfo:ilmnbslookup:IDNotUnique', thisID, filename));
    end
    
    % Put data into structure for output
    for count = 1:max(colMap)
        if colFound(count)
            % if multiple matches are found return as a cell array.
            % Otherwise return as strings.
            if numel(row) > 1
                data = {allData{colMap(count)}{row}};
            else
                data = allData{colMap(count)}{row};
            end
            output(IDCount).(expectedCols{count}) = data;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lookupfield,delimiter] = parse_inputs(varargin)
% Parse the varargin parameter/value inputs

% Check that we have the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:ilmnbslookup:IncorrectNumberOfArguments', mfilename));
end

% The allowed inputs
okargs = {'lookupfield','delimiter'};
lookupfield = 'Search_key';
delimiter = '\t';
% deal with the various inputs
for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % columns
            if ischar(pval)
                lookupfield = pval;
            else
                error(message('bioinfo:ilmnbsread:LookupFieldNotChar'));
            end
        case 2  % delimiter
            if ischar(pval)
                delimiter = pval;
            else
                error(message('bioinfo:ilmnbsread:DelimiterFieldNotChar'));
            end
    end
end

