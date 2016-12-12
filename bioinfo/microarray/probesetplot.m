function probesetplot(celStruct,cdfStruct,theID,varargin)
% PROBESETPLOT plots intensity values for an Affymetrix probe set
%
%   PROBESETPLOT(CELSTRUCT,CDFSTRUCT,PS) plots the PM and MM intensity
%   values for probe set PS. CELSTRUCT is a structure created by the
%   function AFFYREAD from an Affymetrix CEL. CDFSTRUCT is the structure
%   created by the function AFFYREAD from the CDF file corresponding to the
%   CEL file. PS can be the index of the probe set or the probe set name.
%   Note that the probe set numbers for a CDF file typically use 0-based
%   indexing, whereas MATLAB uses 1-based indexing, so
%   CDFStruct.ProbeSets(1) has ProbeSetNumber 0.
%
%   PROBESETPLOT(...,'GENENAME',true) uses the gene name, rather than probe
%   set name, for the title. This option requires that the GIN file
%   associated with the chip type that you are using is in the same
%   directory as the CDF file from which CDFSTRUCT was created.
%
%   PROBESETPLOT(...,'FIELD',FIELDNAME) shows the data for the field
%   FIELDNAME. Valid field names are 'Intensity', 'StdDev', 'Background',
%   'Pixels', and 'Outlier'.
%
%   PROBESETPLOT(...,'SHOWSTATS',true) adds mean and standard deviation
%   lines to the plot.
%
%   Example:
%       celStruct = affyread('Ecoli-antisense-121502.cel');
%       cdfStruct = affyread(...
%                    'C:\Affymetrix\LibFiles\Ecoli_ASv2\Ecoli_ASv2.CDF');
%       % plot the values for probe set 'argG_b3172_at'
%       probesetplot(celStruct,cdfStruct,'argG_b3172_at','showstats',true);
%
%   See also AFFYDEMO, AFFYREAD, CELINTENSITYREAD, PROBESETLINK, PROBESETLOOKUP,
%   PROBESETVALUES.

%   Affymetrix is a registered trademark of Affymetrix, Inc.

%   PROBESETPLOT(...,'PLOTONLY',TF) shows only the intensity plot without
%   the heat map image if TR is set to TRUE (default). 
% 
%   PROBESETPLOT(...,'IMAGEONLY',TF) shows only the false color image
%   without the intensity plot if TR is set to TRUE. The default is FALSE.
%   The false color image of PM and MM probe pairs are displayed in
%   parallel. For SNP array, Allele A PM and MM probe pairs and Allele B PM
%   and M probe pairs are shown in parallel. The left blocks are for sense
%   probes and the right blocks are for antisense probes. If IMAGEONLY
%   option is set to TRUE, the PLOTONLY option will be ignored.

% Copyright 2003-2008 The MathWorks, Inc.


% Check the number of inputs
bioinfochecknargin(nargin,3,mfilename);
% Now check that the CEL and CDF struct are cel and cdf structs

if ~isstruct(celStruct)
    error(message('bioinfo:probesetplot:CelStructNotStruct'));
end

if  ~isstruct(cdfStruct)
    error(message('bioinfo:probesetplot:CdfStructNotStruct'));
end

if ~isfield(celStruct,'Name') || ~isfield(celStruct,'ChipType') ||...
        ~isfield(celStruct,'Probes') || isempty(regexpi(celStruct.Name,'.cel$'))
    error(message('bioinfo:probesetplot:BadCelStruct'));
end
if ~isfield(cdfStruct,'Name') || ~isfield(cdfStruct,'ChipType') ||...
        ~isfield(cdfStruct,'ProbeSets') || isempty(regexpi(cdfStruct.Name,'.CDF$'))
    error(message('bioinfo:probesetplot:BadCdfStruct'));
end

% Check that the ChipType match
if strcmpi(celStruct.ChipType, cdfStruct.ChipType) == 0
    if  strncmpi(celStruct.ChipType, cdfStruct.ChipType, min(numel(celStruct.ChipType),numel(cdfStruct.ChipType)))
        warning(message('bioinfo:probesetplot:ChipTypeMismatch', cdfStruct.ChipType, celStruct.ChipType));
    else
        error(message('bioinfo:probesetplot:ChipTypeMismatch', cdfStruct.ChipType, celStruct.ChipType));
    end
end


% Set defaults
fieldNames = { 'Background',...
    'Intensity',...
    'StdDev',...
    'Pixels',...
    'Outlier',...
    'Masked'};

fieldNamesMapping = {'ProbeSetNumber',...
    'ProbePairNumber',...
    'UseProbePair',...
    'Background',...
    'PMPosX',...
    'PMPosY',...
    'PMIntensity',...
    'PMStdDev',...
    'PMPixels',...
    'PMOutlier',...
    'PMMasked',...
    'MMPosX',...
    'MMPosY',...
    'MMIntensity',...
    'MMStdDev',...
    'MMPixels',...
    'MMOutlier',...
    'MMMasked'};

legendFieldNames =    {'Probe Set Number',...
    'Probe Pair Number',...
    'Use Probe Pair',...
    'Background',...
    'Perfect Match PosX',...
    'Perfect Match PosY',...
    'Perfect Match Intensity',...
    'Perfect Match StdDev',...
    'Perfect Match Pixels',...
    'Perfect Match Outlier',...
    'Perfect Match Masked',...
    'Mismatch PosX',...
    'Mismatch PosY',...
    'Mismatch Intensity',...
    'Mismatch StdDev',...
    'Mismatch Pixels',...
    'Mismatch Outlier',...
    'Mismatch Masked'};
fieldNum = 2;
PMfieldNum = 7;
MMfieldNum = 14;
statsFlag = false;
backgroundFlag = false;
genenameFlag = false;
plotFlag = true;
imageFlag = false;

% get the ID
if ischar(theID)
    ID = strmatch(theID,{cdfStruct.ProbeSets.Name});
    if isempty(ID)
        error(message('bioinfo:probesetplot:UnknownProbeName', theID));
    elseif length(ID)>1
        warning(message('bioinfo:probesetplot:AmbiguousProbeName', theID));
        ID = ID(1);
    end
else
    ID = theID;
end

% deal with the various inputs
if nargin > 3
    if rem(nargin,2) == 0
        error(message('bioinfo:probesetplot:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'fieldname','showstats','genename','plotonly', 'imageonly'};
    for j=1:2:nargin-3
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:probesetplot:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:probesetplot:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % field
                    if ischar(pval)
                        fieldNum = find(strncmpi(pval,fieldNames,length(pval)));
                        if isempty(fieldNum)
                            error(message('bioinfo:probesetplot:UnknownFieldName', pval));
                        elseif length(fieldNum)>1
                            error(message('bioinfo:probesetplot:AmbiguousFieldName', pval));
                        end
                        if fieldNum == 1  % background is a special case
                            backgroundFlag = true;
                            PMfieldNum = find(strncmpi(pval,fieldNamesMapping,length(pval)));
                        else
                            backgroundFlag = false;

                            PMfieldNum = find(strncmpi(['PM' pval],fieldNamesMapping,length(pval)));
                            MMfieldNum = find(strncmpi(['MM' pval],fieldNamesMapping,length(pval)));
                        end
                    else
                        error(message('bioinfo:probesetplot:BadFieldName'));
                    end
                case 2 % showStats
                    statsFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 3 % genename
                    genenameFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 4 % plotonly
                    plotFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 5 % imageonly
                    imageFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                    if imageFlag 
                        plotFlag = false;
                    end
            end
        end
    end
end

% Extract the data from the big struct
try
    theSet = probesetvalues(celStruct,cdfStruct,ID,'background',backgroundFlag);
catch theErr
    if strcmpi(theErr.identifier,'bioinfo:probesetvalues:BadPSID')
        error(message('bioinfo:probesetplot:BadPSID', ID));
    elseif strcmpi(theErr.identifier,'bioinfo:probesetvalues:PSIDNotScalar')
        error(message('bioinfo:probesetplot:PSIDNotScalar'));
    else
        rethrow(theErr)
    end
end
if ~genenameFlag
    theName = cdfStruct.ProbeSets(ID).Name;
else
    try
        lookupResult = probesetlookup(cdfStruct,cdfStruct.ProbeSets(ID).Name);
        theName = lookupResult.Identifier;
    catch ME %#ok
        warning(message('bioinfo:probesetplot:GeneNameLookupFailed', [ cdfStruct.ChipType, '.GIN' ], cdfStruct.LibPath));
        theName = cdfStruct.ProbeSets(ID).Name;
    end
end
% Everything is zero based, so add one -- hopefully this won't confuse
% anyone. Li & Wong do this
if theSet(1,2) == 0
    xData = theSet(:,2)+1;
end

%--------------------
%==Get group infor
grpCol = 19;
groups = theSet(:,grpCol);
uni_groups = unique(groups);
numGrp = numel(uni_groups);
snpFlag = false;
if numGrp == 4
    snpFlag = true;
end 

PM_I = theSet(:,PMfieldNum);
MM_I = theSet(:,MMfieldNum);
NumPairs = cdfStruct.ProbeSets(ID).NumPairs;
%--------------------------
if snpFlag
    grpidx_A = strcmpi(cdfStruct.ProbeSets(ID).GroupNames,...
        cdfStruct.ProbeSets(ID).GroupNames{1});
    grpA = uni_groups(grpidx_A);
    grpB = uni_groups(~grpidx_A);
end
% plot the PM data
hAxis = [];
hfig = figure;
if ~imageFlag || plotFlag
    if snpFlag
        hAxis = plotSNPIntensity(xData, PM_I, MM_I, groups, grpA, grpB, backgroundFlag);
        legendString = {'Perfect Match Intensity',...
            'Mismatch Intensity',...
            'Allele A','Allele B'};
        if backgroundFlag
            legendString(2) = [];
        end
    else
        % plot the PM data
        hPlot = plot(xData,PM_I,'marker','x','markersize',10);
        hold on
        legendString = {legendFieldNames{PMfieldNum}};
        if ~backgroundFlag % plot the mismatch data
            hPlot = plot(xData,MM_I,'r','marker','x','markersize',10);
            legendString = {legendString{:} legendFieldNames{MMfieldNum}}; %#ok<*CCAT>
        end
        hAxis = get(hPlot,'parent');
    end

    set(hAxis,'xtick',xData,'xlim',[0 max(xData)+1],...
        'xticklabel',num2str(xData));
    xlabel('Probe Pair')
    ylabel(fieldNames{fieldNum});
end

titleAx = hAxis;
if snpFlag
    if find(groups == grpA(1), 1) > find(groups == grpA(2), 1)
        grpA = flipud(grpA);
    end
    if find(groups == grpB(1), 1) > find(groups == grpB(2), 1)
        grpB = flipud(grpB);
    end
    MM_A = [MM_I(groups == grpA(1))', MM_I(groups == grpA(2))'];
    PM_A = [PM_I(groups == grpA(1))' PM_I(groups == grpA(2))'];
    PM_B = [PM_I(groups == grpB(1))' PM_I(groups == grpB(2))'];
    MM_B = [MM_I(groups == grpB(1))' MM_I(groups == grpB(2))'];
    I=[MM_A; PM_A; PM_B; MM_B];
else
    I = [PM_I'; MM_I'];
end
if plotFlag && imageFlag
    hfig2 = figure;
    ha = plotImage([], hfig2, I, NumPairs, snpFlag);
    title(ha, sprintf('%s %s',theName,fieldNames{fieldNum}),...
        'interpreter','none');
elseif ~plotFlag && imageFlag
    ha = plotImage([], hfig, I, NumPairs, snpFlag);
    if imageFlag
        titleAx = ha;
    end
elseif ~plotFlag && ~imageFlag
    ha = plotImage(hAxis, hfig, I, NumPairs, snpFlag);
    if imageFlag
        titleAx = ha;
    end
end
title(titleAx, sprintf('%s %s',theName,fieldNames{fieldNum}),...
    'interpreter','none');

% Add on the mean and std and a legend
if statsFlag && ~isempty(hAxis)
    xData = [0;xData;xData(end)+1];
    rep = ones(size(xData));
    mVal = nanmean(theSet(:,PMfieldNum));
    plot(hAxis, xData,mVal(rep),'k-.');
    stdVal = nanstd(theSet(:,PMfieldNum));
    plusOne = mVal + stdVal;
    minusOne = mVal - stdVal;
    plot(hAxis,[xData;NaN;xData],[plusOne(rep);NaN;minusOne(rep)],'m-.');
    % display mismatch stats
    if backgroundFlag
        legendString =  {legendString{:} ,'Mean','STD'};
    else
        mVal = nanmean(theSet(:,MMfieldNum));
        plot(hAxis, xData,mVal(rep),'c-.');
        stdVal = nanstd(theSet(:,MMfieldNum));
        plusOne = mVal + stdVal;
        minusOne = mVal - stdVal;
        plot(hAxis, [xData;NaN;xData],[plusOne(rep);NaN;minusOne(rep)],'g-.');
        legendString =  {legendString{:} ,'PM Mean','PM STD','MM Mean','MM STD'};
    end
end

hold off
if ~isempty(hAxis)
    legend(hAxis, legendString);
end
end

function ha = plotImage(hAxis, hfig, I, NumPairs, snpFlag)
% Plot the image of the probe pairs.
% hAxis-plot axis

%== Grid pixels for each probe
gdPixels = 20;

ha = axes('Parent', hfig);
imageFlag = false;
if isempty(hAxis)
    pos = getpixelposition(ha);
    midPY = pos(2)+ pos(4)/2;
    imageFlag = true;
else
    pos = getpixelposition(hAxis);
    offset = gdPixels*3;
    if snpFlag
        offset = gdPixels*5;
    end
    setpixelposition(hAxis, [pos(1) pos(2)+offset pos(3) pos(4)-offset])
    midPY = 10+ offset/2;
end
midPX = pos(1)+pos(3)/2;

if snpFlag
    x = 1:NumPairs/2;
    y = 1:4;
    height = gdPixels * 4;
    width = gdPixels * NumPairs/2;
else
    x = 1:NumPairs;
    y = 1:2;
    height = gdPixels * 2;
    width = gdPixels * NumPairs;
end
setpixelposition(ha, [midPX-width/2, midPY-height/2 width height])
imagesc(x, y, I, 'Parent', ha)
axis(ha, 'off')
xl = get(ha, 'Xlim');
xt = xl(1)-0.5;
t_pros = {'Parent', ha,'HorizontalAlignment', 'right'};
if snpFlag
    probeLabelStr = {'mmA', 'pmA', 'pmB', 'mmB'};
    mmA_t = text(t_pros{:}, 'String', probeLabelStr{1},'Position', [xt, 1]);%#ok
    pmA_t = text(t_pros{:}, 'String', probeLabelStr{2},'Position', [xt, 2]);%#ok
    pmB_t = text(t_pros{:}, 'String', probeLabelStr{3},'Position', [xt, 3]);%#ok
    mmB_t = text(t_pros{:}, 'String', probeLabelStr{4},'Position', [xt, 4]);%#ok
else
   probeLabelStr = {'pm', 'mm'};
   mm_t = text(t_pros{:}, 'String', probeLabelStr{1},'Position', [xt, 1]);%#ok
   pm_t = text(t_pros{:}, 'String', probeLabelStr{2},'Position', [xt, 2]);%#ok
end

%== Add colorbar
cbarbutton = uigettool(hfig, 'Annotation.InsertColorbar');
set(cbarbutton, 'ClickedCallback', {@updateColorbarCB, hfig, ha, imageFlag})
end

function hAxis = plotSNPIntensity(xData, PM_I, MM_I, groups, grpA, grpB,backgroundFlag)
AS_Idx = groups == grpA(1); % Allele A, sense
AAS_Idx = groups == grpA(2); % Allele A, antisense
BS_Idx = groups == grpB(1); % Allele B, sense
BAS_Idx = groups == grpB(2); % Allele B, antisense
grpStr = cellstr(num2str(groups));
textProps = {'VerticalAlignment', 'Bottom'};
hPlot = plot(xData,PM_I, 'marker','x','markersize',9);
hold on
if ~backgroundFlag % plot the mismatch data
    plot(xData,MM_I, 'r','marker','x','markersize',8);
    text(xData(AS_Idx|AAS_Idx), MM_I(AS_Idx|AAS_Idx), grpStr(AS_Idx|AAS_Idx),...
        'Color','g',textProps{:});
    text(xData(BS_Idx|BAS_Idx), MM_I(BS_Idx|BAS_Idx), grpStr(BS_Idx|BAS_Idx),...
        'Color','m',textProps{:});
end

text(xData(AS_Idx|AAS_Idx), PM_I(AS_Idx|AAS_Idx), grpStr(AS_Idx|AAS_Idx),...
           'Color','g', textProps{:});
text(xData(BS_Idx|BAS_Idx), PM_I(BS_Idx|BAS_Idx), grpStr(BS_Idx|BAS_Idx),...
           'Color','m', textProps{:});
%== Hidden lines for legend
hiddenProps = {'linestyle', 'none',...
               'marker', '.',...
               'visible', 'off'};
line([1 1], [1 1], 'color', 'g', hiddenProps{:})
line([1 1], [1 1], 'color', 'm', hiddenProps{:})
hAxis = get(hPlot,'parent');
end

function updateColorbarCB(hsrc, event, hfig, cax, imageflag)  %#ok<INUSL>
% Insert Colorbar callback function to modify the color bar context menu
% items.
cbmenu = findall(hfig,'Tag','figMenuInsertColorbar');
cbtogg = uigettool(hfig,'Annotation.InsertColorbar');
p_pos = get(cax, 'Position');

if ~isempty(cax)
    if (~isempty(cbmenu) && isequal(cbmenu,gcbo) && strcmpi(get(cbmenu,'checked'),'off')) ||...
            (~isempty(cbtogg) && isequal(cbtogg,gcbo) && strcmpi(get(cbtogg,'state'),'on'))
        cbar = colorbar('peer',cax, 'location', 'EastOutside');
        pos = get(cbar, 'Position');
        if imageflag
            m_pos = [p_pos(1)+p_pos(3)+pos(3) p_pos(2)-p_pos(4) pos(3)*0.5 pos(4)*3];
        else
            m_pos = [p_pos(1)+p_pos(3)+pos(3) p_pos(2) pos(3)*0.5 pos(4)];
        end
        set(cbar, 'Position', m_pos);
        set(cax, 'position', p_pos)
    else
        colorbar('peer',cax,'off');
    end
else
    if ~isempty(cbtogg)
        set(cbtogg,'State','off');
    end
end
end
