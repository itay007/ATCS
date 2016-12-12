function outStruct = affysnpquartets(celtruct,cdfstruct,ps)
% AFFYSNPQUARTETS Creates a table of SNP quartets for a probe set.
%
%   SNPQSTRUCT = AFFYSNPQUARTETS(CELSTRUCT,CDFSTRUCT,PS) creates a struct
%   of SNP results from the probe data in a CEL file structure CELSTRUCT,
%   where PS is a probe set index or probe set name from the CDF library
%   file structure CDFSTRUCT.
%
%   SNPQSTRUCT is a structure with the following fields:
%
%       'ProbeSet'
%       'AlleleA'
%       'AlleleB'
%       'Quartet'
%
%   The 'Quartet' field is a structure array containing PM and MM intensity
%   values for sense and anti-sense probes for alleles A and B.
%
%   Sample SNP data files are available from 
%   http://www.affymetrix.com/support/technical/sample_data/hapmap_trio_data.affx
%
%   Example:
%
%   CELStruct=affyread('NA06985_Hind_B5_3005533.CEL');
%   CDFStruct=affyread('Mapping50K_Hind240.cdf');
%   SNPQStruct = affysnpquartets(CELStruct,CDFStruct,'SNP_A-1684395')
%   SNPQStruct.Quartet(1)
%
%   See Also AFFYREAD, PROBESETVALUES.

% Copyright 2007 The MathWorks, Inc.



bioinfochecknargin(nargin,3,mfilename);
psvals = probesetvalues(celtruct,cdfstruct,ps,'background',false);
numProbePairs = size(psvals,1);
usedMask = false(numProbePairs,1);
PMCol = 7;
MMCol = 14;
groupCol = 19;
dirCol = 20;
numGroups = numel(unique(psvals(:,groupCol)));

quartetCount = 1;
outStruct.ProbeSet = cdfstruct.ProbeSets(psvals(1)+1).Name;
% loop over all entries
for count = 1:numProbePairs
    % figure out if this value has already been used
    if usedMask(count)
        continue
    end
    % make room for the quartet
    theQuartet = zeros(numGroups,2);
    groupCount = 1;
    theQuartet(groupCount,1) = psvals(count,groupCol);
    theQuartet(groupCount,2) = count;
    usedMask(count) = true;
    % now find the corresponding pairs
    for inner = count+1:numProbePairs
        if ~usedMask(inner) && ~any(psvals(inner,groupCol)==theQuartet(:,1))
            groupCount = groupCount+1;
            theQuartet(groupCount,1) = psvals(inner,groupCol);
            theQuartet(groupCount,2) = inner;
            usedMask(inner) = true;
        end
    end
    % set up the outputs first time through
    if quartetCount == 1
        outStruct.AlleleA = cdfstruct.ProbeSets(psvals(1)+1).GroupNames{1};
        outStruct.AlleleB = char(setxor(cdfstruct.ProbeSets(psvals(1)+1).GroupNames,outStruct.AlleleA));
        outStruct.Quartet = struct('A_Sense_PM',[],'B_Sense_PM',[],...
            'A_Sense_MM',[],'B_Sense_MM',[],...
            'A_Antisense_PM',[],'B_Antisense_PM',[],...
            'A_Antisense_MM',[],'B_Antisense_MM',[]);
    end
    % Fill out a table of the quartet
    for outLoop = 1:sum(theQuartet(:,1)~=0)
        if isequal(cdfstruct.ProbeSets(psvals(1)+1).GroupNames{theQuartet(outLoop,1)},outStruct.AlleleA)
            if psvals(theQuartet(outLoop,2),dirCol) == 1
                outStruct.Quartet(quartetCount).A_Sense_PM = psvals(theQuartet(outLoop,2),PMCol);
                outStruct.Quartet(quartetCount).A_Sense_MM = psvals(theQuartet(outLoop,2),MMCol);
            else
                outStruct.Quartet(quartetCount).A_Antisense_PM = psvals(theQuartet(outLoop,2),PMCol);
                outStruct.Quartet(quartetCount).A_Antisense_MM = psvals(theQuartet(outLoop,2),MMCol);
            end
        else
            if psvals(theQuartet(outLoop,2),dirCol) == 1
                outStruct.Quartet(quartetCount).B_Sense_PM = psvals(theQuartet(outLoop,2),PMCol);
                outStruct.Quartet(quartetCount).B_Sense_MM = psvals(theQuartet(outLoop,2),MMCol);
            else

                outStruct.Quartet(quartetCount).B_Antisense_PM = psvals(theQuartet(outLoop,2),PMCol);
                outStruct.Quartet(quartetCount).B_Antisense_MM = psvals(theQuartet(outLoop,2),MMCol);
            end
        end

    end
    quartetCount = quartetCount +1;
end
