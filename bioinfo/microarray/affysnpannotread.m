function output = affysnpannotread(filename, snpID, varargin)
%AFFYSNPANNOTREAD reads an annotation CSV file for Affymetrix Mapping arrays.
%
%   ADATA = AFFYSNPANNOTREAD(FILENAME, PID) reads the annotation CSV file
%   for Mapping 10K arrays; Mapping 100K array set and Mapping 500K array
%   set, and returns a structure containing these fields:
%           ProbeSetIDs
%           Chromosome
%           ChromPosition
%           Cytoband
%           Sequence
%           AlleleA
%           AlleleB
%           Accession
%           FragmentLength       
%   PID can be a single probe set ID or a cell array of the probe set IDs.
%   The data in the structure should be ordered by PID. FILENAME should be
%   the full path to the annotation file.
% 
%   The CSV Annotation files, for example,
%   Mapping50K_Xba240.na25.annot.csv, can be downloaded from
%   http://www.affymetrix.com/support/technical/annotationfilesmain.affx.
%
%   AFFYSNPANNOTREAD(...,'LOOKUPFIELD',LOOKUP) returns annotation for PID
%   from the column specified by LOOKUP. There are more columns in the
%   annotation file than specified in fields of the default output
%   structure. If LOOKUPFIELD is specified, the function returns a
%   structure with the field of LOOKUP and ordered by PID. LOOKUP should be
%   the same as the column header to be looked up.
% 
%   Example:
% 
%       % Read the Mapping50K_Xba240.na25.annot.csv file for Mapping50K Xba SNP array
%       cdf = affyread('C:\AffyLibFiles\Mapping50K_Xba240.CDF');
%       probesetIDs = {cdf.ProbeSets.Name}';
%       snpInfo = affysnpannotread('Mapping50K_Xba240.na25.annot.csv',probesetIDs)
%
%   See also AFFYSNPINTENSITYSPLIT, AFFYREAD.
%
%   Affymetrix is a registered trademark of Affymetrix, Inc.

% Copyright 2008-2011 The MathWorks, Inc.

% Undocummented:
%   PIDS = AFFYSNPANNOTREAD(FILENAME,'') returns a cell aray with all the
%   probe set IDs available in the file.
%   COLH = AFFYSNPANNOTREAD(...,'LOOKUPFIELD','?') returns a cell array with
%   all the column headers available in the file.

% Handle inputs
bioinfochecknargin(nargin, 2, mfilename)
lookupField = parse_inputs(varargin{:});

if strcmp(lookupField,'?')
    readOnlyColumnHeaders = true;
elseif isempty(snpID)
    readOnlyColumnHeaders = false;
    readOnlyProbeSetIds = true;
    lookupField = {'Probe Set ID'};
else
    readOnlyColumnHeaders = false;
    readOnlyProbeSetIds = false;
    %== Filter SNP IDs
    if ischar(snpID)
       snpID = {snpID};
    end
    if ~iscellstr(snpID)
       error(message('bioinfo:affysnpannotread:InvalidIDType'));
    end
end

%== Try to open the file
try
    fopenMessage = '';
    [fid, fopenMessage] = fopen(filename,'rt');
catch ME %#ok<NASGU>
    fid = -1;
end
if fid == -1
    error(message('bioinfo:affysnpannotread:CannotOpenFile', filename, fopenMessage));
end

%== Check for Probe Set ID at the beginning of the first line
tline = fgetl(fid);

% First line(s) do not contain columns names.
while isempty(regexp(tline, 'Probe Set ID', 'once'))
    tline = fgetl(fid);
    
    %== Error if reach the end-of file and still not find 'Probe Set ID'
    if tline == -1
        fclose(fid);
         error(message('bioinfo:affysnpannotread:WrongAnnotFileFormat', filename));
    end
end

% Read the second line and extract the column names
colNames = strread(tline,'%q','delimiter',',');

if readOnlyColumnHeaders
    output = colNames;
    return
end

% Number of columns
numCols = numel(colNames);

 
defaultFlag = true;
%==Check if lookup field match one of the column names
if ~isempty(lookupField)
    lookupCol = find(strcmpi(lookupField,colNames));
    if isempty(lookupCol)
        fclose(fid);
        error(message('bioinfo:affysnpannotread:LookupFieldNotValid', lookupField));
    end
    if lookupCol == 1
        output = struct('ProbeSetID', '');
    else
        lookupCol = [1 lookupCol];
        lookupField = {'Probe Set ID', lookupField};
        output = struct('ProbeSetIDs', '',...
            strrep(lookupField{end}, ' ', '_'), '');
    end
    defaultFlag = false;
else % default outputs
    lookupField = {'Probe Set ID',...
        'Chromosome',...
        'Physical Position',...
        'Cytoband',...
        'Flank',...
        'Allele A',...
        'Allele B',...
        'Associated Gene',...
        'Fragment Enzyme'}';
    % Fragment Enzyme Type Length Start Stop - n25, Fragment Enzyme
    % Length Start Stop - n24, check Fragment Enzyme only.
    numDefaultFields = length(lookupField);
    lookupCol = ones(numDefaultFields, 1);
    for i = 1:numDefaultFields
        idx = strmatch(lookupField{i}, colNames);
        if isempty(idx)
            fclose(fid);
            error(message('bioinfo:affysnpannotread:DefaultLookupFieldNotValid', lookupField{ i }, filename));
        end
        lookupCol(i) = idx(1);
    end
    
    if isempty(lookupCol)
        fclose(fid);
        error(message('bioinfo:affysnpannotread:NotExpectedSNPAnnot', filename));
    end
    
    output = struct('ProbeSetIDs', '',...
        'Chromosome', '',...
        'ChromPosition', '',...
        'Cytoband', '',...
        'Sequence', '',...
        'AlleleA', '',...
        'AlleleB', '',...
        'Accession', '',...
        'FragmentLength', '');
end



%== Format string
formatStr = repmat('%*q', 1, numCols);
%== Only read the lookup columns
formatStr(lookupCol*3-1)=[];
fPos = ftell(fid);
bufSize = 4095;
numTries = 0;
maxTries = 10;
while numTries < maxTries;
    try
        allData = textscan(fid,formatStr,'Delimiter',',','Bufsize',bufSize);
        break
    catch theErr % deal with the case where a very long line exists
        if strcmp(theErr.identifier,'MATLAB:textscan:BufferOverflow')
            bufSize = (bufSize*2)+1;
            numTries = numTries; %#ok
            fseek(fid,fPos,-1);
        else
            fclose(fid);
            rethrow(theErr);
        end
    end
end
fclose(fid);

%== Find snpID in probesetID
probesetID = allData{1};

if readOnlyProbeSetIds
    output = probesetID;
    return
end

[tf_snpid, loc] = ismember(snpID, probesetID);

% == Error if no snpID found in proberstID
if sum(tf_snpid)==0
   error(message('bioinfo:affysnpannotread:ProbeSetIDNotMatch'));
end

%== Put data into structure for output
outputFields = fieldnames(output);
numOutFields = numel(outputFields);

for i = 1:numOutFields
    tmp_data = repmat({' '}, length(snpID),1);
    data = allData{i};
    tmp_data(tf_snpid) = data(loc(tf_snpid));
    allData{i} = tmp_data;
end

% Extra only the fragment length 
if defaultFlag
    %== Convert chromosome to numbers with X to be 23, Y to be 24 and probe
    %   sets with no chromosome number to -1
    %== G751409: Convert MT to 25
    allData{2} = strrep(allData{2}, ' ', '---');
    tmpstr = strtrim(allData{2});
    tmpstr = strrep(tmpstr, 'X', '23');
    tmpstr = strrep(tmpstr, 'Y', '24');
    tmpstr = strrep(tmpstr, 'MT', '25');
    tmpstr = strrep(tmpstr, '---', '-1');
    allData{2} = int8(sscanf(sprintf('%s ',tmpstr{:}),'%d')); 
    
    %== Convert chromosome genome positions to numbers
    allData{3} = strrep(allData{3}, ' ', '---');
    tmpstr = strtrim(allData{3});
    tmpstr = strrep(strtrim(tmpstr), '---', '-1');
    allData{3} = sscanf(sprintf('%s ',tmpstr{:}),'%d'); 
    
    %==Get only the accession number from the 'Associated Gene' field
    allData{end-1} = regexp(allData{end-1}, '\w+_*\d+', 'match', 'once');
    
    %== Extract and convert only PCR fragment length
    tmpstr = strtrim(regexp(allData{end}, '\s\d+\s', 'match', 'once'));
       % There are ProbeSets with not fragment length (these retun a '')
    allData{end} = sscanf(sprintf('0%s ',tmpstr{:}),'%d');
end

%== Put data into structure for output
for i = 1:numOutFields
    output.(outputFields{i}) = allData{i};
end

end % affysnpannotread

function lookupfield = parse_inputs(varargin)
% Parse the varargin parameter/value inputs

% Check that we have the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:affysnpannotread:IncorrectNumberOfArguments', mfilename));
end

% The allowed inputs
okargs = {'lookupfield'};
lookupfield = [];

% Deal with the inputs
for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % columns
            if ischar(pval)
                lookupfield = pval;
            else
                error(message('bioinfo:affysnpannotread:LookupFieldNotChar'));
            end
    end
end
end % parse_input
