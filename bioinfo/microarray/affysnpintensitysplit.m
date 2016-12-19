function outStruct = affysnpintensitysplit(probeStruct, varargin)
%AFFYSNPINTENSITYSPLIT splits the Affymetrix SNP probe intensity matrix
%  for allele A and allele B.
% 
%   PSOUT = AFFYSNPINTENSITYSPLIT(PSIN) splits the Affymetrix Mapping SNP
%   array probe intensity matrix into intensity matrices for allele A and
%   allele B from the probe intensity structure PSIN returned by the
%   CELINTENSITYREAD function. The returned structure PSOUT contains these
%   fields:
%       CDFName
%       CELNames
%       NumChips
%       NumProbeSets
%       NumProbes (changed to number of allele probes)
%       ProbeSetIDs
%       ProbeIndices (changed to allele probe indices)
%       PMAIntensities (PM intensities for allele A)
%       PMBIntensities (PM intensities for allele B)
%       MMAIntensities (If input structure contains MMIntensities.)
%       MMBIntensities (If input structure contains MMIntensities.)
%
%   Note: The output probes do not include the control probes.
% 
%   AFFYSNPINTENSITYSPLIT(PSIN,'CONTROLS', TF) returns the control probe
%   intensities if TF=TRUE. Control probes sometimes only contain
%   information for one allele. In this case the value for the
%   corresponding allele A or B that is not present is set to NaN. Default
%   TF=FALSE. 
%
%   Example:
%
%       % Read the CEL files in current directory and its CDF file. 
%       ps = celintensityread('*','Mapping250K_Nsp.CDF');
%       % Extract PM probe intensities for allele A and allele B. 
%       ps_snp = affysnpintensitysplit(ps);
% 
%   See also AFFYSNPANNOTREAD, AFFYREAD, CELINTENSITYREAD.
% 
%   Affymetrix is a registered trademark of Affymetrix, Inc.

%   AFFYSNPINTENSITYSPLIT(PSIN, 'GWARRAY', TF) indicates the array is an
%   Affymetrix GenomeWide SNP Arrays 5.0 or 6.0. By default, the function
%   checks if the CDFName field in the PSIN structure contains
%   'GenomeWide', and guesses the type of the SNP array.
% 
%   AFFYSNPINTENSITYSPLIT(PSIN, 'SNPONLY', TF) returns the SNP probe
%   intensities and an extra field called CNVData containing the probe
%   intensities of copy number variation, if TF=FALSE. This option only
%   works for Affymetrix GenomeWide SNP Arrays 5.0 or 6.0. Default TF=TRUE.
% 
%   AFFYSNPINTENSITYSPLIT(PSIN, 'CNVONLY', TF) returns only the probe
%   intensities for the detection of copy number variation (CNV) in the
%   CNVData field in the output structure if TF=TRUE. This option only
%   works for Affymetrix GenomeWide SNP Arrays 5.0 or 6.0. Default
%   TF=FALSE.
%
%   Note: For Affymetrix GenomeWide SNP Arrays 5.0 or 6.0, if the CONTROLS
%   option is true. the control probe intensities are returned in the
%   Controls field in the output structure.

%   Copyright 2008-2009 The MathWorks, Inc.


%== Check inputs
bioinfochecknargin(nargin, 1, mfilename)
%== Check the input structure
checkIntensityStruct(probeStruct)

%== Initializing
inPV = parse_inputs(varargin{:});
mmFlag = isfield(probeStruct,'MMIntensities');

%== Figure out the array type
if isempty(inPV.GWAFlag)
    inPV.GWAFlag = strncmpi(probeStruct.CDFName, 'GenomeWide', 10);
end

%==Find all the CN probesets which contain a single probe
if inPV.GWAFlag
    CNIdx = strncmp(probeStruct.ProbeSetIDs, 'CN_', 3);
    nCNProbes =  sum(CNIdx);
    
    if nCNProbes == 0
        warning(message('bioinfo:affysnpintensitysplit:NOCNVProbes'));
        inPV.GWAFlag = false;
    end
end

if ~inPV.GWAFlag
    inPV.CNVOnlyFlag = false;
end

%== Find probe set indices
probesetIndices = find(probeStruct.ProbeIndices == 0);
probePairs = diff([probesetIndices;probeStruct.NumProbes+1]);

%== Find probes for allele A and allele B.
row_A = probeStruct.GroupNumbers == 1 |...
        probeStruct.GroupNumbers == 3 |...
        probeStruct.GroupNumbers ==0;
row_B = probeStruct.GroupNumbers == 2 |...
        probeStruct.GroupNumbers == 4;
nProbe_A = sum(row_A);
nProbe_B = sum(row_B);

%== The number of probes for each allele and the number of probes
nProbes = max(nProbe_A, nProbe_B);
nProbeSets = probeStruct.NumProbeSets;

%==Create intensity holders
intClass = class(probeStruct.PMIntensities);
PMA = nan(nProbes, probeStruct.NumChips, intClass);
PMB = nan(nProbes, probeStruct.NumChips, intClass);

if mmFlag
    MMA = nan(nProbes, probeStruct.NumChips, intClass);
    MMB = nan(nProbes, probeStruct.NumChips, intClass);
end
%== Get the splitted probe indices
probeIndices = zeros(nProbes,1, 'uint8');
PM = probeStruct.PMIntensities;

if mmFlag
    MM = probeStruct.MMIntensities;
end

%== Loop through each probe set and split the probes for allele A and B
pstart = 1;
controlsProbeSet = false(nProbeSets,1);
controlsProbe = false(nProbes,1);

if inPV.GWAFlag
    CNVProbes = false(nProbes,1);
end

for i = 1 : nProbeSets
    ppstart = probesetIndices(i);  % start probe pair index in a probe set
    ppend = ppstart+probePairs(i)-1; % end prope pair index in a probe set
    pidx = ppstart:ppend;
    groups = probeStruct.GroupNumbers(pidx);
    if inPV.GWAFlag && CNIdx(i)
        npp= 1:(ppend-ppstart+1);
        pend = pstart + npp(groups == 1)-1;
        if isempty(pend)
            pend = pstart;
        end
        CNVProbes(pend) = true;
    elseif ~inPV.GWAFlag || ( inPV.GWAFlag && ~CNIdx(i))
        if all(groups==1) || all(groups == 0)
            controlsProbeSet(i) = true;
            pend = pstart+ppend-ppstart-1;
            controlsProbe(pstart:pend) = true;
        else
            row_A = groups == 1 | groups == 3;
            
            nA = sum(row_A);
            nB = probePairs(i) - nA;
            pend = pstart + max(nA,nB)-1;
            
            probeIndices(pstart:pend) = 0:(pend-pstart);
            PMA(pstart:pend, :) = PM(pidx(row_A), :);
            
            if nB > 0
                PMB(pstart:pend, :) = PM(pidx(~row_A), :);
            end
            
            if mmFlag
                MMA(pstart:pend, :) = MM(pidx(row_A), :);
                if nB > 0
                    MMB(pstart:pend, :) = MM(pidx(~row_A), :);
                end
            end
        end
    end
    pstart = pend + 1;
end

%== Update probeStruct for output
if inPV.GWAFlag
    outStruct.CDFName = probeStruct.CDFName;
    outStruct.CELNames = probeStruct.CELNames;
    outStruct.NumChips = probeStruct.NumChips;
    
    if ~inPV.CNVOnlyFlag
        outStruct.NumProbeSets = sum(~controlsProbeSet) - nCNProbes;
        outStruct.ProbeSetIDs = probeStruct.ProbeSetIDs(~controlsProbeSet & ~ CNIdx);
        outStruct.NumProbes = sum(~controlsProbe & ~CNVProbes );
        outStruct.ProbeIndices = probeIndices(~controlsProbe & ~CNVProbes);
        outStruct.PMAIntensities = PMA(~controlsProbe & ~CNVProbes, :);
        outStruct.PMBIntensities = PMB(~controlsProbe & ~CNVProbes, :);
    end
    
    if inPV.CNVOnlyFlag || ~inPV.SNPOnlyFlag
        outStruct.CNVData.NumProbeSets = nCNProbes;
        outStruct.CNVData.ProbeSetIDs = probeStruct.ProbeSetIDs(CNIdx);
        outStruct.CNVData.PMIntensities = PMA(CNVProbes);
    end
    
    if inPV.controlsFlag 
        outStruct.Controls.NumProbeSets = sum(controlsProbeSet);
        outStruct.Controls.ProbeSetIDs = probeStruct.ProbeSetIDs(controlsProbeSet);
        outStruct.Controls.NumProbes = sum(controlsProbe);
        outStruct.Controls.ProbeIndices = probeIndices(controlsProbe);
        outStruct.Controls.PMIntensities = PMA(controlsProbe);
    end
else
    probeStruct.NumProbeSets = sum(~controlsProbeSet);
    probeStruct.ProbeSetIDs(controlsProbeSet) = [];
    probeStruct.NumProbes = sum(~controlsProbe);
    probeStruct.ProbeIndices = probeIndices(~controlsProbe);
    probeStruct.PMAIntensities = PMA(~controlsProbe, :);
    probeStruct.PMBIntensities = PMB(~controlsProbe, :);
    if mmFlag
        probeStruct.MMAIntensities =  MMA(~controlsProbe, :);
        probeStruct.MMBIntensities = MMB(~controlsProbe, :);
    end
    
    probeStruct = rmfield(probeStruct, {'PMIntensities', 'GroupNumbers'});    
    if mmFlag
        probeStruct = rmfield(probeStruct, 'MMIntensities');
    end

    outStruct = probeStruct;
end
end % affysnpintensitysplit

%----------------Helper functions ----------
function checkIntensityStruct(cbStruct)
validFields = {'PMIntensities',...
               'GroupNumbers',...
               'ProbeIndices',...
               'NumChips',...
               'NumProbeSets'};

m = isfield(cbStruct, validFields); 

if ~all(m)
    missingFields = validFields(~m);
    if numel(missingFields)== numel(validFields)
        error(message('bioinfo:affysnpintensitysplit:MissingAllRequiredFields'));
    else
        mf = ''; % list with required fields
        for i=1:numel(missingFields)
            mf = [mf missingFields{i} ', ']; %#ok<AGROW>
        end
        mf = mf(1:end-2); % remove extra comma and space from the end
        error(message('bioinfo:affysnpintensitysplit:MissingRequiredField', mf));
    end
end      
end % checkIntensityStruct

function inPV = parse_inputs(varargin)
% Parse the varargin parameter/value inputs

% Check that we have the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:affysnpintensitysplit:IncorrectNumberOfArguments', mfilename));
end

% The allowed inputs
okargs = {'controls', 'gwarray',  'snponly', 'cnvonly'};
inPV.controlsFlag = false;
inPV.SNPOnlyFlag = true;
inPV.CNVOnlyFlag = false;
inPV.GWAFlag = [];

% Deal with the inputs
for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % controls
            inPV.controlsFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 2 % gwarray
            inPV.GWAFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 3  % SNP only - output SNP intensities only
            inPV.SNPOnlyFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 4  % CNV only - output CNV intensities only
            inPV.CNVOnlyFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end
end
end % parse_input


