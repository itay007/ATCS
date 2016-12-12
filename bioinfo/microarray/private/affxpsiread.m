function out = affxpsiread(filename)
%AFFXPSIREAD reads Affymetrix PSI files.
%
%   psiStruct = AFFXPSIREAD(FILE) reads an Affymetrix PSI file FILE
%   and creates a structure psiStruct. PSI files contain information about
%   the names of the probe sets. This information is also stored in the CDF
%   file but the PSI file can be used for fast lookup.
%
%   See also AFFYREAD, AGFEREAD, CELINTENSITYREAD, GPRREAD,
%   PROBELIBRARYINFO, PROBESETLINK, PROBESETLOOKUP, PROBESETPLOT,
%   PROBESETVALUES, SPTREAD.
%
%   GeneChip and Affymetrix are registered trademarks of Affymetrix, Inc.

% Copyright 2007-2012 The MathWorks, Inc.


% get the parts of the file

[~,theName] = fileparts(filename);

out.ChipType = theName;

% open the file
fid = fopen(filename,'r');
if fid == -1
    error(message('bioinfo:affxpsiread:invalidPsiFile', filename));
end
% The first line should contain the number of probes
numProbeSetsCheck = [];
try
    header = fgetl(fid);
    if ~strncmpi(header,'#Probe Sets:',12)
        warning(message('bioinfo:affxpsiread:BadPsiHeader'));
        numProbeSetsCheck = strread(header,'%*s%d','delimiter',':');
    end
    theData = textscan(fid,'%*d%s%d');
    fclose(fid);
    numProbeSets = numel(theData{1});
    if ~isempty(numProbeSetsCheck) && numProbeSets ~= numProbeSetsCheck
        warning(message('bioinfo:affxpsiread:NumProbeSetsMismatch'));
    end
    
    out.NumProbeSets = numProbeSets;
    out.ProbeSets = theData{1};
    out.NumProbesUsed = theData{2};
    
catch ME
    error(message('bioinfo:affxpsiread:badPsiFile', filename))
end






