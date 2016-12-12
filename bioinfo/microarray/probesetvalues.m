function [table,lookup] = probesetvalues(celStruct,cdfStruct,theID,varargin)
% PROBESETVALUES extracts probe set values from probe results
%
%   PSVALUES = PROBESETVALUES(CELSTRUCT,CDFSTRUCT,PS) creates a table of
%   values for probe set PS from the probe data in a CEL file structure
%   CELSTRUCT, where PS is a probe set index or probe set name from the CDF
%   library file structure CDFSTRUCT.
%
%   PSVALUES is a matrix with 20 columns and one row for each probe pair in
%   the probe set. The columns correspond to these fields:
%
%     'ProbeSetNumber'
%     'ProbePairNumber'
%     'UseProbePair'
%     'Background'
%     'PMPosX'
%     'PMPosY'
%     'PMIntensity'
%     'PMStdDev'
%     'PMPixels'
%     'PMOutlier'
%     'PMMasked'
%     'MMPosX'
%     'MMPosY'
%     'MMIntensity'
%     'MMStdDev'
%     'MMPixels'
%     'MMOutlier'
%     'MMMasked'
%     'GroupNumber'
%     'Direction'
%
%   PROBESETVALUES(...,'BACKGROUND',BGVAL) allows you to pass a vector of
%   pre-calculated background values or to disable the background
%   calculation. If BGVAL is a vector whose length is equal to the number
%   of probes in CELSTRUCT then the values in BGVAL will be used as the
%   background values. If BGVAL is false then the background values will
%   not be calculated. The default behavior is to calculate the background
%   value. This can be slow if you need to lookup the value of many probes.
%
%   If PROBESETVALUES is called with no input arguments it returns this list
%   of column names as a cell array.
%   The 'UseProbePair' column is for backwards compatibility and is not
%   currently used.%
%
%   Example:
%       celStruct = affyread('Ecoli-antisense-121502.cel');
%       cdfStruct = affyread(...
%                    'C:\Affymetrix\LibFiles\Ecoli_ASv2\Ecoli_ASv2.CDF');
%       % get the values for probe set 'argG_b3172_at'
%       [baData,zones,background] = zonebackadj(celStruct,'cdf',cdfStruct);
%       psvals = probesetvalues(celStruct,cdfStruct,'argG_b3172_at',...
%                      'background',background)
%
%   See also AFFYDEMO, AFFYREAD, CELINTENSITYREAD, PROBELIBRARYINFO,
%   PROBESETLINK, PROBESETLOOKUP, PROBESETPLOT, RMABACKADJ.

%   Affymetrix and NetAffx are registered trademarks of Affymetrix, Inc.

% Copyright 2003-2008 The MathWorks, Inc.


% Thank you to Andrew Knyazev and Abram Jujunashvili of  University of
% Colorado at Denver for improvements to the implementation of this
% function.


if nargin == 0
    table = {
        'ProbeSetNumber'
        'ProbePairNumber'
        'UseProbePair'
        'Background'
        'PMPosX'
        'PMPosY'
        'PMIntensity'
        'PMStdDev'
        'PMPixels'
        'PMOutlier'
        'PMMasked'
        'MMPosX'
        'MMPosY'
        'MMIntensity'
        'MMStdDev'
        'MMPixels'
        'MMOutlier'
        'MMMasked'
        'GroupNumber'
        'Direction'
        };
    lookup = [];
    return
end

bioinfochecknargin(nargin,3,mfilename)

% Now check that the CEL and CDF struct are cel and cdf structs

if ~isstruct(celStruct)
    error(message('bioinfo:probesetvalues:CelStructNotStruct'));
end

if  ~isstruct(cdfStruct)
    error(message('bioinfo:probesetvalues:CdfStructNotStruct'));
end

if ~isfield(celStruct,'Name') || ~isfield(celStruct,'ChipType') ||...
        ~isfield(celStruct,'Probes') || isempty(regexpi(celStruct.Name,'.cel$'))
    error(message('bioinfo:probesetvalues:BadCelStruct'));
end
if ~isfield(cdfStruct,'Name') || ~isfield(cdfStruct,'ChipType') ||...
        ~isfield(cdfStruct,'ProbeSets') || isempty(regexpi(cdfStruct.Name,'.CDF$'))
    error(message('bioinfo:probesetvalues:BadCdfStruct'));
end

% Check that the ChipType match
if strcmpi(celStruct.ChipType, cdfStruct.ChipType) == 0
    if  strncmpi(celStruct.ChipType, cdfStruct.ChipType, min(numel(celStruct.ChipType),numel(cdfStruct.ChipType)))
        warning(message('bioinfo:probesetvalues:ChipTypeMismatch', cdfStruct.ChipType, celStruct.ChipType));
    else
        error(message('bioinfo:probesetvalues:ChipTypeMismatch', cdfStruct.ChipType, celStruct.ChipType));
    end
end

backFlag = true;
backVals = [];

% deal with the various inputs
if nargin > 3
    if rem(nargin,2) == 0
        error(message('bioinfo:probesetvalues:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'background'};
    for j=1:2:nargin-3
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:probesetvalues:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:probesetvalues:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % noback- Undocumented, only used by PROBESETPLOT (DEFAULT=TRUE)
                    if numel(pval) == celStruct.NumProbes
                        backFlag = true;
                        backVals = pval;
                    elseif isscalar(pval)
                        backFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                    else
                        error(message('bioinfo:probesetvalues:BadBackgroundValue'));
                    end
            end
        end
    end
end
persistent  affyProbeSetNames  affyCDFStructName

if isempty(affyProbeSetNames) || ~isequal(affyCDFStructName, cdfStruct.Name)
    affyProbeSetNames = {cdfStruct.ProbeSets.Name};
    affyCDFStructName = cdfStruct.Name;
end

% get the ID
if ischar(theID)
    ID = find(strncmp(theID,affyProbeSetNames,numel(theID)));
    if isempty(ID)
        error(message('bioinfo:probesetvalues:UnknownProbeName', theID));
    elseif length(ID)>1
        warning(message('bioinfo:probesetvalues:AmbiguousProbeName', theID));
        ID = ID(1);
    end
else
    ID = theID;
end

numCols = cdfStruct.Cols;

% Check that the ID is valid
if ~isscalar(ID)
    error(message('bioinfo:probesetvalues:PSIDNotScalar'));
end
if ID > numel(cdfStruct.ProbeSets) || ID < 1 ||  floor(ID)~= ID
    error(message('bioinfo:probesetvalues:BadPSID', ID));
end
NumPairs = cdfStruct.ProbeSets(ID).NumPairs;
table = zeros(NumPairs,18);
lookup = zeros(NumPairs,2);
thePairs = cdfStruct.ProbeSets(ID).ProbePairs;

PMXCol = strcmp('PMPosX',cdfStruct.ProbeSetColumnNames);
PMYCol = strcmp('PMPosY',cdfStruct.ProbeSetColumnNames);
MMXCol = strcmp('MMPosX',cdfStruct.ProbeSetColumnNames);
MMYCol = strcmp('MMPosY',cdfStruct.ProbeSetColumnNames);
GroupCol = strcmp('GroupNumber',cdfStruct.ProbeSetColumnNames);
DirectionCol = strcmp('Direction',cdfStruct.ProbeSetColumnNames);

IntensityCol = strcmp('Intensity',celStruct.ProbeColumnNames);
StdDevCol = strcmp('StdDev',celStruct.ProbeColumnNames);
PixelsCol = strcmp('Pixels',celStruct.ProbeColumnNames);
OutlierCol = strcmp('Outlier',celStruct.ProbeColumnNames);
MaskedCol = strcmp('Masked',celStruct.ProbeColumnNames);


PMX = thePairs(:,PMXCol);
PMY = thePairs(:,PMYCol);
PMRow = PMY*numCols + PMX +1;
MMX = thePairs(:,MMXCol);
MMY = thePairs(:,MMYCol);
MMRow = MMY*numCols + MMX + 1;

table(:,1) = ID-1;
table(:,2) = 0:NumPairs-1;
% table(:,3) = 0; % UseProbePair -- currently unused
% table(:,4) = 0; % backgrounds -- fill in later if necessary
table(:,5) = PMX;
table(:,6) = PMY;
table(:,7) = celStruct.Probes(PMRow,IntensityCol);
table(:,8) = celStruct.Probes(PMRow,StdDevCol);
table(:,9) = celStruct.Probes(PMRow,PixelsCol);
table(:,10) = celStruct.Probes(PMRow,OutlierCol);
table(:,11) = celStruct.Probes(PMRow,MaskedCol);
table(:,12) = MMX;
table(:,13) = MMY;
table(:,14) = celStruct.Probes(MMRow,IntensityCol);
table(:,15) = celStruct.Probes(MMRow,StdDevCol);
table(:,16) = celStruct.Probes(MMRow,PixelsCol);
table(:,17) = celStruct.Probes(MMRow,OutlierCol);
table(:,18) = celStruct.Probes(MMRow,MaskedCol);
table(:,19) = thePairs(:,GroupCol);
table(:,20) = thePairs(:,DirectionCol);

% Create the lookup values
lookup(:,1) = PMRow;
lookup(:,2) = MMRow;

if backFlag
    if isempty(backVals)
        [adjVals,zones,backgrounds] = zonebackadj(celStruct,'cdf',cdfStruct,'bgindices',lookup(:,1));
        table(:,4) = backgrounds;
    else
        table(:,4) = backVals(lookup(:,1));
    end
end

