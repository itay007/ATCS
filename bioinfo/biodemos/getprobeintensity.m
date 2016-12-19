function intensity = getprobeintensity(celStruct, cdfStruct, numProbes, type)
% GETPROBEINTENSITY is a help function to read out the probe intensity from a
% Affymetrix CEL structure.

% Copyright 2006 The MathWorks, Inc.

intensity = zeros(numProbes, 1);
numCols = cdfStruct.Cols;
pairCount = 0;

colIdx = 0;
if type == 1 % For PM probe intensity 
    colIdx = 3;
elseif type == 2
    colIdx = 5; % for MM probe intensity
end
    
for i = 1:cdfStruct.NumProbeSets
    numPairs = cdfStruct.ProbeSets(i).NumPairs;
    thePairs = cdfStruct.ProbeSets(i).ProbePairs;
    PX = thePairs(:,colIdx);
    PY = thePairs(:,colIdx + 1);
    intensity(pairCount+1 : pairCount + numPairs, 1) = celStruct.Probes(PY*numCols + PX + 1, 3);
    pairCount = pairCount + numPairs;
end