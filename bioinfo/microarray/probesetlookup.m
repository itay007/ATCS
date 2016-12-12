function outStruct = probesetlookup(affyStruct,theID)
% PROBESETLOOKUP looks up the gene name for a probe set
%
%   PROBESETLOOKUP(AFFYSTRUCT,ID) returns the probe set information for
%   probe set or gene ID from CHP or CDF structure affyStruct. ID can be a
%   string or cell array of strings.
%
%   This function requires that you have the GIN file that is associated
%   with the chip type you are using. The GIN file needs to be in the same
%   directory as the CHP or CDF file from which the AFFYSTRUCT was created.
%
%   Example:
%       cdfStruct = affyread(...
%                    'C:\Affymetrix\LibFiles\Ecoli_ASv2\Ecoli_ASv2.CDF');
%       % Look up based on the probe set ID
%       probesetlookup(cdfStruct,'argG_b3172_at')
%       % Or look up based on the gene ID
%       probesetlookup(cdfStruct,'3315278')
%
%   See also AFFYDEMO, AFFYREAD, CELINTENSITYREAD, PROBELIBRARYINFO,
%   PROBESETLINK, PROBESETPLOT, PROBESETVALUES, RMABACKADJ.

%   Affymetrix is a registered trademarks of Affymetrix, Inc.

% Copyright 2003-2008 The MathWorks, Inc.


bioinfochecknargin(nargin,2,mfilename)


if (iscellstr(theID)&& length(theID) > 1)
    numIDs = numel(theID);
    for idCount = numIDs:-1:1
        outStruct(idCount) = probesetlookup(affyStruct,theID{idCount}); %#ok<AGROW>
    end
return
end



% Now check that the AFFYSTRUCT is a struct and CDF or CEL structs
if  ~isstruct(affyStruct)
    error(message('bioinfo:probesetlookup:CdfStructNotStruct'));
end

if ~isfield(affyStruct,'Name') || ~isfield(affyStruct,'ChipType') || ~isfield(affyStruct,'ProbeSets')
    error(message('bioinfo:probesetlookup:BadAffyStruct'));
end

% Check that theID is a string
if ~ischar(theID)
    error(message('bioinfo:probesetlookup:BadIDInput'));
end


persistent ginStruct affyProbeSetNames affyStructName

libPath = affyStruct.LibPath;
chipType = affyStruct.ChipType;

if isempty(ginStruct) ||~isfield(ginStruct,'Name') ||...
        ~isequal(ginStruct.Name, chipType)
    try
        if exist(fullfile(libPath,[chipType '.GIN']),'file')
            ginStruct = affyread([chipType '.GIN'],libPath);
        else
            ginStruct = affyread([chipType '.gin'],libPath);
        end
    catch allExceptions
        error(message('bioinfo:probesetlookup:BadGINFile', chipType, libPath));
    end
end

if isempty(affyProbeSetNames)
    affyProbeSetNames = {affyStruct.ProbeSets.Name};
    affyStructName = affyStruct.Name;
else
    if ~isequal(affyStructName,affyStruct.Name) || ...
            ~isequal(numel(affyProbeSetNames),numel(affyStruct.ProbeSets))
        affyProbeSetNames = {affyStruct.ProbeSets.Name};
        affyStructName = affyStruct.Name;
    end
end


% initialize output values to empty
geneName = '';
probeSetName = '';
cIndex = '';
description = '';
source = '';
sourceUrl = '';

% see if this is a probe name
geneID = find(strncmp(theID,ginStruct.ProbeSetName,numel(theID)));
if isempty(geneID)
    geneID = find(strncmp(theID,ginStruct.ID,numel(theID)));
end

try
    geneName = ginStruct.ID{geneID};
    % if the ID is empty try the first chunk of the Description using / as
    % a delimiter.
    if isempty(geneName)
        [firstField,theRest] = strtok(ginStruct.Description{geneID},'/');
        if ~isempty(theRest)
            geneName = firstField;
        end
    end
    probeSetName = ginStruct.ProbeSetName{geneID};
    cIndex = find(strncmp(probeSetName,affyProbeSetNames,numel(probeSetName)));

    description = ginStruct.Description{geneID};
    if numel(ginStruct.SourceID) > 1
        sourceID = ginStruct.SourceID(geneID);
    else
        sourceID = 1;
    end
    source = ginStruct.SourceNames{sourceID};
    sourceUrl = sprintf(ginStruct.SourceURL{sourceID},geneName);
catch allExceptions %#ok<NASGU>
    if isempty(geneID)
        error(message('bioinfo:probesetlookup:UnknownProbeName', theID));
    elseif numel(geneID) > 1
        error(message('bioinfo:probesetlookup:AmbiguousProbeName', theID));
    end
end
outStruct.Identifier = geneName;
outStruct.ProbeSetName = probeSetName;
outStruct.CDFIndex = cIndex;
outStruct.GINIndex = geneID;
outStruct.Description = description;
outStruct.Source = source;
outStruct.SourceURL = sourceUrl;
