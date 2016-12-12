function table = probelibraryinfo(CELStruct,CDFStruct)
% PROBELIBRARYINFO extracts probe set library information for probe results
%
%   PROBEINFO = PROBELIBRARYINFO(CELSTRUCT,CDFSTRUCT) creates a table of
%   information linking the probe data in a CEL file structure, CELSTRUCT,
%   with probe set information from a CDF file structure, CDFSTRUCT.
%
%   PROBEINFO is a matrix with three columns and the same number of rows as
%   the probes field of the CELSTRUCT. The first column is the index of the
%   probe set to which the probe belongs. The second column contains the
%   probe pair number, and the third column indicates whether the probe is
%   a perfect match (1) or mismatch (-1) probe. Probes that do not
%   correspond to a probe set in the CDF library file have probe set ID 0.
%
%   Note that Affymetrix probe pair indexing is 0-based but MATLAB
%   indexing is 1-based. The output from PROBELIBRARYINFO is 1-based.
%
%   Example:
%       celStruct = affyread('Ecoli-antisense-121502.cel');
%       cdfStruct = affyread(...
%                    'C:\Affymetrix\LibFiles\Ecoli_ASv2\Ecoli_ASv2.CDF');
%       probeinfo = probelibraryinfo(celStruct,cdfStruct);
%       % Find out which probe set the 1104th probe belongs to.
%       cdfStruct.ProbeSets(probeinfo(1104,1)).Name
%
%   See also AFFYDEMO, AFFYREAD, CELINTENSITYREAD, PROBESETLINK,
%   PROBESETLOOKUP, PROBESETPLOT, PROBESETVALUES.

%   Affymetrix and NetAffx are registered trademarks of Affymetrix, Inc.

% Copyright 2003-2008 The MathWorks, Inc.


bioinfochecknargin(nargin,2,mfilename);

% Now check that the CEL and CDF struct are cel and cdf structs
if ~isstruct(CELStruct)
    error(message('bioinfo:probelibraryinfo:CELStructNotStruct'));
end

if  ~isstruct(CDFStruct)
    error(message('bioinfo:probelibraryinfo:CDFStructNotStruct'));
end

if ~isfield(CELStruct,'Name') || ~isfield(CELStruct,'ChipType') ||...
        ~isfield(CELStruct,'Probes') || isempty(regexpi(CELStruct.Name,'.cel$'))
    error(message('bioinfo:probelibraryinfo:BadCELStruct'));
end
if ~isfield(CDFStruct,'Name') || ~isfield(CDFStruct,'ChipType') ||...
        ~isfield(CDFStruct,'ProbeSets') || isempty(regexpi(CDFStruct.Name,'.CDF$'))
    error(message('bioinfo:probelibraryinfo:BadCDFStruct'));
end

% Check that the ChipType match
if strcmpi(CELStruct.ChipType, CDFStruct.ChipType) == 0
    error(message('bioinfo:probelibraryinfo:ChipTypeMismatch', CDFStruct.ChipType, CELStruct.ChipType));
end

numCols = double(max(CELStruct.Probes(:,1)))+1;

numProbeSets = numel(CDFStruct.ProbeSets);
table = zeros(size(CELStruct.Probes,1),3);

for count = 1:numProbeSets
    NumPairs = CDFStruct.ProbeSets(count).NumPairs;
    thePairs = CDFStruct.ProbeSets(count).ProbePairs;
    for inner = 1:NumPairs
        PMX = thePairs(inner,3);
        PMY = thePairs(inner,4);
        PMRow = PMY*numCols + PMX +1;
        table(PMRow,:) = [count,inner,1];
        MMX = thePairs(inner,5);
        MMY = thePairs(inner,6);
        MMRow = MMY*numCols + MMX + 1;
        table(MMRow,:) = [count,inner,-1];
    end
end
