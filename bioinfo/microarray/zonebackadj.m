function [ adjVal,zones,backgrounds] = zonebackadj(celfile,varargin)
%ZONEBACKADJ performs the background adjustment of Affymetrix microarray
% probe-level data using a zone based method.
%
%   Background intensities for cells are estimated by using the standard
%   deviation of the lowest 2% of intensities in K rectangular zones. These
%   estimates are smoothed by using a weighted average of the distance from
%   the cells to the centers of the zones. These values are then adjusted
%   using a noise correction to ensure that the adjusted values do not
%   become zero.
%
%   ZONEBACKDATA = ZONEBACKADJ(DATA) returns the background adjusted
%   intensities of a chip file. DATA can be the CEL structure or array of
%   CEL structures returned by AFFYREAD.
%
%   [ZONEBACKDATA, ZONESTRUCT] = ZONEBACKADJ(DATA) returns a structure
%   containing the centers of the zones used and the estimates of the
%   background values at the centers of each zone.
%
%   [ZONEBACKDATA, ZONESTRUCT, BACKGROUND] = ZONEBACKADJ(DATA) returns the
%   estimated background values.
%
%   ZONEBACKADJ(...,'NUMZONES',K) allows you to specify the number of zones
%   used. If K is a scalar then it must be a square number. To specify a
%   non-square grid use a 2 element array where K(1) gives the number of
%   rows in the grid and K(2) the number of columns. The default value is
%   16.
%
%   ZONEBACKADJ(...,'PERCENT',P) allows you to specify P where the lowest P
%   percent values from each zone are used to estimate the background for
%   the zone. The default is 2.
%
%   ZONEBACKADJ(...,'SMOOTHFACTOR',SFACTOR) allows you to specify the
%   smoothing factor used in the calculation of the weighted average of the
%   contributions of each zone to the background of a point. The default
%   value is 100.
%
%   ZONEBACKADJ(...,'NOISEFRAC',NF) allows you to specify the noise
%   fraction NF such that the adjusted value is given by:
%      max((Intensity - weightedBackground), NF*localNoiseEstimate)
%   The default is 0.5.
%
%   ZONEBACKADJ(...,'CDF',CDF) allows you to specify a CDF file or CDF
%   struct created with AFFYREAD. This is used to identify which cells are
%   control cells or masked cells. Typically control and masked cells are
%   not used in the background estimates.
%
%   ZONEBACKADJ(...,'MASK',MASK) allows you to specify a logical vector of
%   the cells to be used in the background calculation. Typically control
%   and masked cells are not used in the background estimates. If the input
%   is a CEL structure then the Masked column in the Probes field is used.
%   Otherwise the default is false(size(DATA,1)).
%
%   ZONEBACKADJ(...,'SHOWPLOT',SP) creates an image of the background
%   estimates and zone centers.
%
%   Example:
%
%       % Read in a CEL file
%       celStruct = affyread('Ecoli-antisense-121502.CEL')
%
%       % Background adjust data from one chip
%       backadjData = zonebackadj(celStruct,'cdf',...
%                      'C:\Affymetrix\LibFiles\Ecoli_ASv2\Ecoli_ASv2.CDF');
%
%
%   See also AFFYREAD, AFFYINVARSETNORM, CELINTENSITYREAD, GCRMA,
%   GCRMABACKADJ, PROBELIBRARYINFO, PROBESETLINK,  PROBESETLOOKUP,
%   PROBESETVALUES, RMABACKADJ, RMASUMMARY.


% References:
% [1] Statistical Algorithms Description Document
% www.affymetrix.com/support/technical/whitepapers/sadd_whitepaper.pdf

% Copyright 2007 The MathWorks, Inc.


%
%   ZONEBACKADJ(...,'BGINDICES',INDICES) allows you to specify a the indices
%   of values for which you want the backgrounds.

% Validate input data
bioinfochecknargin(nargin,1,mfilename);

if ~isstruct(celfile) || ~isfield(celfile,'ChipType')
    error(message('bioinfo:zonebackadj:InputNotCelStruct'))

end
if numel(celfile) > 1
    % if all the ChipTypes are the same then return as array
    returnAsMatrix =  numel(unique({celfile.ChipType})) ==1;

    for outerLoop = numel(celfile):-1:1
        [ adjVal{outerLoop},zones(outerLoop),backgrounds{outerLoop}] =...
            zonebackadj(celfile(outerLoop),varargin{:}); 
    end
    if returnAsMatrix 
        adjVal = cell2mat(adjVal);
        backgrounds = cell2mat(backgrounds);
    else
        warning(message('bioinfo:zonebackadj:MixedChipType'));
    end
    return
end

numZones = 16; %#ok
zoneRows = 4;
zoneCols = 4;

numProbes = size(celfile.Probes(:,1),1);
maskField = strcmpi(celfile.ProbeColumnNames,'Masked');
if any(maskField)
    celmask = celfile.Probes(:,maskField);
else
    celmask = false(numProbes,1);
end


pct = .02;
smooth = 100;
noiseFrac = 0.5;
useCDF = false;
userMask = false;
showPlot = false;
customBGIndices = false;
% deal with the various inputs
if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:zonebackadj:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'numzones', 'percent','smoothfactor','noisefrac','cdffile','mask', 'showplot','bgindices'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:zonebackadj:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:zonebackadj:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % numZones
                    if ~isnumeric(pval)
                        error(message('bioinfo:zonebackadj:ZoneNotNumeric'));
                    end
                    if numel(pval)>2 || any(pval<=0) || any(pval ~= floor(pval))
                        error(message('bioinfo:zonebackadj:BadNumZones'));
                    end
                    numZones = pval;
                    if isscalar(numZones)
                        if sqrt(numZones) ~= ceil(sqrt(numZones))
                            warning(message('bioinfo:zonebackadj:NumZonesNotSquare', ceil( sqrt( numZones ) ) .^ 2, numZones));
                        end
                        zoneRows = ceil(sqrt(numZones));
                        zoneCols = zoneRows;
                    else
                        zoneRows = numZones(1);
                        zoneCols = numZones(2);
                    end


                case 2 % percent
                    if ~isnumeric(pval)
                        error(message('bioinfo:zonebackadj:PercentNotNumeric'));
                    end
                    if ~isscalar(pval) || pval<=0 || pval > 100
                        error(message('bioinfo:zonebackadj:BadPercent'));
                    end
                    pct = pval/100;
                case 3 % smoothingfactor
                    if ~isnumeric(pval)
                        error(message('bioinfo:zonebackadj:SmoothingFactorNotNumeric'));
                    end
                    if ~isscalar(pval) || pval<=0
                        error(message('bioinfo:zonebackadj:SmoothingFactorNotScalar'));
                    end
                    smooth = pval;
                case 4 % noisefrac
                    if ~isnumeric(pval)
                        error(message('bioinfo:zonebackadj:NoiseFracNotNumeric'));
                    end
                    if ~isscalar(pval) || pval< 0
                        error(message('bioinfo:zonebackadj:NoiseFracNotScalar'));
                    end
                    noiseFrac = pval;
                case 5 % cdffile
                    cdfStruct = pval;
                    useCDF = true;
                    if ~isstruct(cdfStruct)
                        try
                            cdfStruct = affyread(cdfStruct);
                        catch theException
                            error(message('bioinfo:zonebackadj:CDFReadError', cdfStruct, theException.message));
                        end
                    end
                    cdfmask = ~isExpProbe(cdfStruct);
                    if size(cdfmask,1) ~= numProbes
                        error(message('bioinfo:zonebackadj:CDFMismatch'));
                    end
                case 6 % mask
                    celmask = pval(:);
                    userMask = true;
                    if size(celmask,1) ~= numProbes
                        error(message('bioinfo:zonebackadj:MaskSizeMismatch'));
                    end
                    try
                        if ~islogical(celmask)

                            celmask= logical(celmask);
                        end
                    catch allExceptions 
                        error(message('bioinfo:zonebackadj:MaskNotLogical'));
                    end
                case 7 % showplot
                    showPlot = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 8 % bgindices
                    bgIndices = pval;
                    customBGIndices = true;
            end
        end
    end
end

if useCDF
    if userMask
        warning(message('bioinfo:zonebackadj:CDFIgnored'))
    else
        celmask = cdfmask;
    end
end

% set up the zones
numRows = celfile.Rows;
numCols = celfile.Cols;
theProbes = celfile.Probes(~celmask,:);
rowBounds = 0:numRows/zoneRows:numRows;
rowCenters = filter([1/2 1/2],1,rowBounds);
colBounds =  0:numCols/zoneCols:numCols;
colCenters = filter([1/2 1/2],1,colBounds);
[~,rowBin] = histc(theProbes(:,1),rowBounds);
[~,colBin] = histc(theProbes(:,2),colBounds);

% estimate the background for the zones and estimate the noise
vals = zeros(zoneRows,zoneCols);
localNoise =  zeros(zoneRows,zoneCols);
for iRow = 1:zoneRows
    for iCol = 1:zoneCols
        mask = (rowBin==iRow) & (colBin == iCol);
        localProbes = sort(theProbes(mask,3));

        vals(iRow,iCol) = mean(localProbes(1:floor(numel(localProbes)*pct)));
        localNoise(iRow,iCol) = std(localProbes(1:floor(numel(localProbes)*pct)));
    end
end

% turn everything back into vectors
zoneY = repmat(rowCenters(2:end),zoneCols,1);
zoneY = zoneY(:);

zoneX = repmat(colCenters(2:end)',zoneRows,1);
zoneX = zoneX(:);
zones.Centers = [zoneX zoneY];
vals = vals(:);
localNoise =localNoise(:);
zones.Background = vals;
cellVals = celfile.Probes;


% set up the output data and pre-allocate the outputs
if ~customBGIndices
    bgIndices = 1:size(cellVals,1);
end

backgrounds = bgIndices(:);
noise = backgrounds;
Iprime = max(celfile.Probes(:,3),0.5);
adjVal = backgrounds;

% loop over the required values
for count = 1:numel(bgIndices)
    weights = 1./(((cellVals(bgIndices(count),1)-zoneX).^2 + (cellVals(bgIndices(count),2)-zoneY).^2)+smooth);
    backgrounds(count) = (1./sum(weights)).*weights'*vals;
    noise(count)= (1./sum(weights)).*weights'*localNoise;
    adjVal(count) = max(Iprime(bgIndices(count))-backgrounds(count),noiseFrac*noise(count));
end

if showPlot
    plotStruct.FullPathName = 'backgrounds.cel';
    plotStruct.ProbeColumnNames = {'PosX','PosY','Backgrounds'};
    plotStruct.Probes = [celfile.Probes(:,[1 2]),backgrounds];
    h = maimage(plotStruct,'Backgrounds','Title','Estimated Background');
    set(h,'Tag','zoneBackAdjPlot');
    hold on
    line(zoneX,zoneY,'linestyle','none','marker','x','color','k');
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function isExpressionProbe = isExpProbe(cdfStruct)

numRows = cdfStruct.Rows;
numCols = cdfStruct.Cols;
numProbes = numRows.*numCols;

expProbes = vertcat(cdfStruct.ProbeSets(1:cdfStruct.NumProbeSets).ProbePairs);


% use the name information in case we ever change the row order
pmXCol = strcmpi(cdfStruct.ProbeSetColumnNames,'PMPosX');
pmYCol = strcmpi(cdfStruct.ProbeSetColumnNames,'PMPosY');
mmXCol = strcmpi(cdfStruct.ProbeSetColumnNames,'MMPosX');
mmYCol = strcmpi(cdfStruct.ProbeSetColumnNames,'MMPosY');

% allocate the output mask
isExpressionProbe = false(numProbes,1);
% Do the perfect match then mismatch
pM = 1+ expProbes(:,pmXCol) + numRows.*expProbes(:,pmYCol);
isExpressionProbe(pM) = true;
mM = 1+ expProbes(:,mmXCol) + numRows.*expProbes(:,mmYCol);
isExpressionProbe(mM) = true;
