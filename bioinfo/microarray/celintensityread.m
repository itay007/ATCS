function [probeintenstruct, cdfStruct] = celintensityread(CELFiles, CDFFile, varargin)
%CELINTENSITYREAD reads probe intensities from Affymetrix CEL files.
%
%   PROBEINTENSTRUCT = CELINTENSITYREAD(CELFILES, CDFFILE) reads Affymetrix
%   CELFILES and the associated library CDFFILE, and creates a structure
%   containing probe intensities, probe indices, and the corresponding
%   probe set IDs. CELFILES can be a cell array of CEL file names. CDFFILE
%   is a string of the CDF file name. If CELFILES is set to wildcard '*',
%   then all the CEL files in the current directory will be read. If
%   CELFILES or CDFFILE is set to '', an Open File dialog will open for you
%   to select CEL files or a CDF file interactively. 
%   The returned structure PROBEINTENSTRUCT with these fields
%       CDFName
%       CELNames
%       NumChips
%       NumProbeSets
%       NumProbes
%       ProbeSetIDs
%       ProbeIndices
%       GroupNumbers
%       PMIntensities 
%       MMIntensities (optional)
%   The ProbeSetIDs and GroupNumbers are ordered as in the CDF file, and
%   ProbeIndices are ordered according to ProbesetIDs. The probe intensity
%   matrices, PMIntensities and MMIntensities, rows correspond to
%   ProbeIndices, and columns correspond to CEL files as read in.
%
%   CELINTENSITYREAD(..., 'CELPATH', CELPATH) allows you to specify the
%   directory where the CELFILES are stored.
%
%   CELINTENSITYREAD(..., 'CDFPATH', CDFPATH) allows you to specify the
%   directory where the CDFFILE is stored.
%
%   CELINTENSITYREAD(..., 'PMONLY', TF) returns a structure with perfect
%   match (PM) probe intensities only if TF is true. If FALSE is specified
%   both perfect match (PM) and mismatch (MM) probe intensities are
%   returned. Default is TRUE.
%
%   CELINTENSITYREAD(...,'VERBOSE',false) turns off verbose output. Default
%   is true.
%
%   Example:
%
%       % Read in all the CEL files in current directory and a CDF file in
%       % a specified directory
%       pmStruct = celintensityread('*','HG_U95Av2.CDF',...
%                            'cdfpath', 'D:\Affymetrix\LibFiles\HGGenome');
%
%       % Select CEL files and a CDF file using Open File dialogs
%       pmStruct = celintensityread('', '');
%
%   See also AFFYGCRMA, AFFYINVARSETNORM, AFFYPREPROCESSDEMO,
%   AFFYPROBESEQREAD, AFFYREAD, AFFYRMA, AGFEREAD, GCRMA, GCRMABACKADJ,
%   GPRREAD, PROBELIBRARYINFO, PROBESETLINK, PROBESETLOOKUP, PROBESETPLOT,
%   PROBESETVALUES, RMABACKADJ, RMASUMMARY, SPTREAD.
%
%   Affymetrix is a registered trademarks of Affymetrix, Inc.

%   [PROBEINTENSTRUCT, CDFSTRUCT] = CELINTENSITYREAD(...) returns the
%   CDF structure. 

% Copyright 2005-2012 The MathWorks, Inc.

probeintenstruct = [];
bioinfochecknargin(nargin,2,mfilename);

% Initialization
celpath = pwd; % current directory
libpath = '';
pmonlyFlag = true;
verbose = true;

if iscellstr(CELFiles)
    celfiles = CELFiles;
elseif ischar(CELFiles)
    celfiles = cellstr(CELFiles);
else
    error(message('bioinfo:celintensityread:InvalidCELFileNames'));
end

if ~ischar(CDFFile) || (iscellstr(CDFFile)&& length(CDFFile) > 1)
    error(message('bioinfo:celintensityread:CDFFILENOTASTRING'));
else
    cdffile = CDFFile;
end

% Deal with the various inputs
if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:celintensityread:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'celpath', 'cdfpath', 'pmonly','verbose'};
    for j=1:2:nargin-3
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:celintensityread:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:celintensityread:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % celpath
                    if ~ischar(pval)
                        error(message('bioinfo:celintensityread:CELPATHNOTASTRING', upper( char( okargs( k ) ) )));
                    end
                    celpath = pval;
                case 2 % cdfpath
                    if ~ischar(pval)
                        error(message('bioinfo:celintensityread:CDFPATHNOTASTRING', upper( char( okargs( k ) ) )));
                    end
                    libpath = pval;
                case 3 % pmonly flag
                    pmonlyFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 4 % verbose flag
                    verbose = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            end
        end
    end
end

% Open the celfiles UI get-file when celfiles = '' or ' '
if length(celfiles)==1 && isempty(strtrim(celfiles{:}))
    [celfiles, celpath] = uigetfile({'*.cel;*.CEL','CEL Files'}, 'Select CEL Files', 'MultiSelect', 'on');
    if iscell(celpath)
        celpath = celpath{:};
    end
    
    if isnumeric(celfiles) && celfiles == 0 % user canceled
        return;
    end
elseif length(celfiles)==1 && strcmpi('*', strtrim(celfiles{:})) %* wildcard
    c = dir([celpath,filesep,'*.cel']);
    % deal with case sensivity on non-Windows platforms
    if ~ispc
        c = [c;dir([celpath,filesep,'*.CEL'])];
    end
    if isempty(c)
        error(message('bioinfo:celintensityread:NotCELFilesExistInCELPath', celpath));
    end
    celfiles = {c.name};
end

%Open the cdffile UI get-file when cdfile = '' or ' ' 
if isempty(strtrim(cdffile))
    [cdffile, libpath] = uigetfile({'*.CDF;*cdf','CDF Files'}, 'Select CDF File');
    if cdffile == 0 % user canceled
        return;
    end
elseif isempty(libpath) % Figure out if the cdf file is in the current directory or on the path

    wcdffile = which(cdffile);
    if ~isempty(wcdffile)
        cdffile = wcdffile;
    end
    
    [thePath,theName,theExt] = fileparts(cdffile);
    if isempty(thePath)
        thePath = pwd;
    end
    
    if isempty(theExt)
        theExt = '.CDF';
    end
    cdffile = [theName,upper(theExt)];
    if ~exist(cdffile,'file') && exist([theName,lower(theExt)],'file');
        cdffile = [theName,lower(theExt)];
    end
    
    libpath = thePath;
    
end

% Figure out the cdf file path
fullcdf = [libpath,filesep,cdffile];

% Read cdf file
if verbose
    fprintf('Reading CDF file: %s\n', cdffile);
end

try
    cdfStruct = affyread(fullcdf);
catch theException
    % catch error throw by open_seqview and rethrow
    msgId = 'bioinfo:celintensityread:InvalidCDFFile';
    newException = MException(msgId,'%s',getString(message(msgId,fullcdf)));
    throw(addCause(newException,theException))
end
nProbeSets = cdfStruct.NumProbeSets;
nProbes = sum([cdfStruct.ProbeSets.NumPairs]);

probesetIDs = {cdfStruct.ProbeSets.Name}';
% Remove QC probesets names
probesetIDs(nProbeSets+1 : end)= [];

if iscell(celfiles)
    nArrays = numel(celfiles);
else
    nArrays = 1; % only one array
    celfiles = {celfiles}';
end

probeIndices = zeros(nProbes,1, 'uint8');
groupNumber = zeros(nProbes,1, 'uint8');
pmMatrix = zeros(nProbes, nArrays,'single');
if ~pmonlyFlag
    mmMatrix = zeros(nProbes, nArrays,'single');
end

% Create probesetindexes vector and probe indexes vector
probeCount=0;
for i = 1:nProbeSets
    numPairs = cdfStruct.ProbeSets(i).NumPairs;
    probeIndices(probeCount+1 : probeCount + numPairs) = (0:numPairs-1)';
    groupNumber(probeCount+1 : probeCount + numPairs) = cdfStruct.ProbeSets(i).ProbePairs(1:numPairs,1);
    probeCount = probeCount + numPairs;
end

% Create probe intensity matrix
for i=1:nArrays
    celname = celfiles{i};
    if verbose
        disp(sprintf('Reading file %d of %d: %s', i, nArrays, celname));
    end
    [~,~,theExt] = fileparts(celname);
    
    if isempty(theExt)
        celname = [celname, '.cel']; %#ok
    end
    
    fullcel = [celpath,filesep,celname];
    if ~exist(fullcel,'file') && exist([celpath,filesep,celfiles{i}, '.CEL'],'file')
        celname = [celfiles{i}, '.CEL'];
        fullcel = [celpath,filesep,celname];
    end
    try
        celStruct = affyread(fullcel);
    catch theException
        msgId = 'bioinfo:celintensityread:InvalidCELFile';
        newException = MException(msgId,'%s',getString(message(msgId,fullcel)));
        throw(addCause(newException,theException))
    end

    % Check the chiptypes are the same as the CDF file. Throw warning if
    % not match.
    if ~strcmpi(celStruct.ChipType, cdfStruct.ChipType)
        warning(message('bioinfo:celintensityread:ChipTypeDoestNotMatch', celStruct.ChipType, celname, cdfStruct.ChipType, cdffile));
    end
    
    % Check the numbers of Rows and Cols are the same in both files.
    if celStruct.Rows ~= cdfStruct.Rows || celStruct.Cols ~= cdfStruct.Cols
       error(message('bioinfo:celintensityread:RowsOrColsNotMatch', celStruct.Rows, celStruct.Cols, celname, cdfStruct.Rows, cdfStruct.Cols, cdffile));
    end 

    pmMatrix(:,i) = getProbeIntensity(celStruct, cdfStruct, nProbes, 1);

    if ~pmonlyFlag
        mmMatrix(:,i) = getProbeIntensity(celStruct, cdfStruct, nProbes, 2);
    end
end

if nargout >= 1
    probeintenstruct.CDFName = cdffile;
    probeintenstruct.CELNames = cellstr(strtok(celfiles, '.'))';
    probeintenstruct.NumChips = nArrays;
    probeintenstruct.NumProbeSets = nProbeSets;
    probeintenstruct.NumProbes =  nProbes;
    probeintenstruct.ProbeSetIDs = probesetIDs;
    %     probeintenstruct.GenBankIDs = genbankIDs;
    %     probeintenstruct.ProbeSetIndexes = probesetIndexes;
    probeintenstruct.ProbeIndices = probeIndices;
    probeintenstruct.GroupNumbers = groupNumber;
    probeintenstruct.PMIntensities = pmMatrix;

    if ~pmonlyFlag
        probeintenstruct.MMIntensities = mmMatrix;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%5
function intensity = getProbeIntensity(celStruct, cdfStruct, numProbes, type)
% GETPROBEINTENSITY is a help function to read out the probe intensity from a
% Affymetrix CEL structure.
intensity = zeros(numProbes, 1,'single');
numCols = cdfStruct.Cols;
pairCount = 0;

PMXCol = strcmp('PMPosX',cdfStruct.ProbeSetColumnNames);
PMYCol = strcmp('PMPosY',cdfStruct.ProbeSetColumnNames);
MMXCol = strcmp('MMPosX',cdfStruct.ProbeSetColumnNames);
MMYCol = strcmp('MMPosY',cdfStruct.ProbeSetColumnNames);

if type == 1 % For PM probe intensity
    XcolIdx = PMXCol;
    YcolIdx = PMYCol;
elseif type == 2
    XcolIdx = MMXCol; % for MM probe intensity
    YcolIdx = MMYCol;
end

for i = 1:cdfStruct.NumProbeSets
    numPairs = cdfStruct.ProbeSets(i).NumPairs;
    thePairs = cdfStruct.ProbeSets(i).ProbePairs;
    PX = thePairs(1:numPairs,XcolIdx);
    PY = thePairs(1:numPairs,YcolIdx);
    idx = pairCount+ (1:numPairs);
    intensity(idx, 1) = celStruct.Probes(PY*numCols + PX + 1, 3);
    intensity(idx( PX==0 | PY==0 ),1) = NaN; %% See g632121
    pairCount = pairCount + numPairs;
end


