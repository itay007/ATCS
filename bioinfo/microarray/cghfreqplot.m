function varargout = cghfreqplot(cghdata, varargin)
%CGHFREQPLOT displays frequency of copy number alterations across multiple
% samples.  
% 
%   S = CGHFREQPLOT(DATA) displays the frequency of copy number gain or
%   loss across multiple samples for each clone on an array against their
%   genomic position along the chromosomes. The frequency is calculated as
%   the fraction of samples with the log2 ratios greater (gain) or less
%   (loss) than a specified threshold. DATA is a structure containing the
%   following fields:
%       Sample (optional)
%       Chromosome
%       GenomicPosition
%       Log2Ratio
%   or a matrix of data with the first column corresponding to chromosome
%   number, the second column corresponding to the genomic position, and
%   the third and higher columns corresponding to the log2 ratios of test
%   to reference intensities. If the 'Sample' field is missing, a default
%   sample name will be assigned to output, S. S is the frequency data
%   structure with the following fields: 
%       Group
%       Chromosome
%       GenomicPosition
%   The Group field is a structure array containing alteration frequency
%   information for each group with following fields:
%       Sample
%       GainFrequency
%       LossFrequency
%   The output structures can be used as input DATA.
% 
%   CGHFREQPLOT(..., 'THRESHOLD', THRESH) sets the gain/loss threshold. A
%   clone is considered to be a gain if its log2 ratio is above THRESH and
%   loss if its log2 ratio is below negative THRESH. If THRESH is a
%   positive scalar, it applies to all the samples. If THRESH is a
%   two-element vector, the first element is the threshold for gain, and
%   the second element is the threshold for loss. If it is a vector of the
%   same length as number of samples,  each element in the vector is
%   considered as a unique gain and loss threshold for each sample. The
%   default is 0.25.
% 
%   CGHFREQPLOT(..., 'GROUP', GRP) specifies the sample groups to calculate
%   the frequency from.  GRP can be a cell array of vectors of sample
%   column indices. Each element in the cell array is considered a group.
%   GRP can also be a vector of sample column indices. The samples
%   specified in the vector are considered a group. The default group is
%   all the samples in DATA.
% 
%   CGHFREQPLOT(..., 'SUBGRP', TF) analyzes samples by subgroups if TF is
%   set to TRUE (default) or together if TF is FALSE. Default is TRUE. 
%  
%   CGHFREQPLOT(..., 'SUBPLOT', TF) displays all the plots in one figure
%   window when more than one subgroup is analyzed if TF is TRUE (default).
%   Otherwise, plots for each subgroup will be displayed in separate
%   windows. 
% 
%   CGHFREQPLOT(..., 'CUTOFF', CUTOFF) plots only the clones with frequency
%   of gains or losses greater than or equal to CUTOFF. If CUTOFF is a
%   scalar, it applies to both gain and loss frequency. If CUTOFF is a
%   two-element vector, the first element is the cutoff for gains, the
%   second element for losses. The default value is 0.
% 
%	CGHFREQPLOT(..., 'CHROMOSOME', CHR) displays the frequency plots of
%	only the chromosome(s) specified by CHR. It can be a single chromosome
%	number or a vector of chromosome numbers. By default all the
%	chromosomes in DATA will be analyzed.
% 
%   CGHFREQPLOT(..., 'INCLUDEX', TF) includes the X chromosome in the
%   analysis if TF is TRUE (default).
% 
%   CGHFREQPLOT(..., 'INCLUDEY', TF) includes the Y chromosome in the
%   analysis if TF is TRUE.  The default is FALSE.
% 
%   CGHFREQPLOT(..., 'CHROMINFO', S) provides cytogenetic banding
%   information for the chromosomes. S can be a structure returned by the
%   CYTOBANDREAD function or the file name of an NCBI ideogram text file or
%   a UCSC Genome Browser cytoband text file. The default is Homo sapiens
%   cytogenetic banding information. 
% 
%   Note: You can download cytogenetic banding information files from the
%   NCBI or the UCSC Genome Browser ftp sites. For example, the current
%   ideogram for Homo sapiens can be downloaded from:
%   ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/mapview/ideogram.gz or
%   ftp://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBandIdeo.txt
%   .gz
% 
%   CGHFREQPLOT(..., 'SHOWCENTR', TF) shows the centromere positions in the
%   plot if TF is TRUE (default). The centromere genomic positions are
%   obtained from CHROMINFO.
% 
%   CGHFREQPLOT(...,'COLOR',C) allows you to specify a color scheme for the
%   vertical lines of the plot indicating the frequency of the gains and
%   losses. C can be the name of or function handle to a function that
%   returns a colormap, or an M-by-3 matrix containing RGB values. The
%   first row of the RGB color matrix is for gain frequency and the second
%   row is for loss frequency.  For example, C = [0 1 0;1 0 0] specifies
%   green for gain and red for loss. The default color scheme is a range of
%   colors from pure green (gain = 1) through yellow (0) to pure red (loss
%   = -1).  
% 
%   CGHFREQPLOT(..., 'YLIM', YLIM) sets the y vertical limits for the plot.
%   YLIM is a two-element vector specifying the minimum and maximum values
%   of the vertical axis. The default is [-1 1].
% 
%   CGHFREQPLOT(..., 'TITLES', TITLES) specifies titles for plots of the
%   groups analyzed. TITLES can be a single string or a cell array of
%   strings.
% 
%   Examples:
%
%       % Plot the frequency of gain and loss for the Coriell study.
%       load coriell_baccgh
%       S = cghfreqplot(coriell_data);
% 
%       % Plot the frequency plot of the Pancreatic cancer study in green/red.
%       load pancrea_oligocgh
%       cghfreqplot(pancrea_data, 'Color', [0 1 0; 1 0 0])
% 
%       % Plot the frequency plots of the two groups in pancreatic data
%       grp1 = strncmp('PA.C', pancrea_data.Sample,4);
%       grp1_ind = find(grp1);
%       grp2 = strncmp('PA.T', pancrea_data.Sample,4);
%       grp2_ind = find(grp2);
%       SP = cghfreqplot(pancrea_data, 'Group', {grp1_ind, grp2_ind},...
%                        'Title', {'CL', 'PT'}, 'Cutoff', 0.25);
%        
%       % Frequency plot of group 'CL' on chromosome 4. 
%       cghfreqplot(pancrea_data, 'Group', grp1_ind,...
%                        'Title', 'CL Group on Chr 4', 'Chromosome', 4);
%       % Add the chromosome 4 ideogram to this frequency plot
%       chromosomeplot('hs_cytoBand.txt', 4, 'addtoplot', gca, 'unit', 2)
% 
%   See also AFFYSNPCNVDEMO, BACACGHDEMO, CGHCBS, CHROMOSOMEPLOT,
%   CYTOBANDREAD.

%   Copyright 2007-2012 The MathWorks, Inc.


%== Input check
bioinfochecknargin(nargin,1,mfilename)

%== Parse input options
appdata = parse_inputs(varargin{:});

%== Check input DATA and validate input parameters
[S, appdata] = check_inputdata(cghdata, appdata);

%== Analyze and plot
if ~appdata.plotonlyFlag
    S = computeFrequency(S, appdata);
end
plotFrequency(S, appdata);

%== Output
if nargout >  0
    varargout{1} = S;
end
end

function S = computeFrequency(cghdata, appdata)
% Compute the gain/loss frequency for each clone to be displayed
ngrp = numel(appdata.group);

gstruct.Sample=[];
gstruct.GainFrequency =[];
gstruct.LossFrequency =[];
gstruct = repmat(gstruct, 1, ngrp);

% Get valid chromosome indices
chrIdx = ismember(cghdata.Chromosome, appdata.dispChrom);
nPoints = sum(chrIdx); 
for i = 1:ngrp
    sampleIDs = appdata.group{i};
    nSamples = appdata.nsamples;
    
    if appdata.subgrpFlag
        nSamples = numel(sampleIDs); 
    end
    
    threshold = appdata.threshold;
    threshold_l = [];
    if numel(appdata.threshold) == 2 && nSamples ~=2
        threshold = appdata.threshold(1);
        threshold_l = appdata.threshold(2);
    elseif numel(appdata.threshold) > 2 || (numel(appdata.threshold) == 2 && nSamples ==2)
        threshold = repmat(appdata.threshold(sampleIDs), nPoints, 1);
    end
    
    gainData = sum((cghdata.Log2Ratio(chrIdx, sampleIDs) > threshold), 2)/nSamples;
    if ~isempty(threshold_l)
        threshold = abs(threshold_l);
    end
    lossData = -sum((cghdata.Log2Ratio(chrIdx, sampleIDs)< - threshold), 2)/nSamples;
    
    gstruct(i).Sample = cghdata.Sample(sampleIDs);
    gstruct(i).GainFrequency = gainData;
    gstruct(i).LossFrequency = lossData;
end

S.Group = gstruct;
S.Chromosome = cghdata.Chromosome(chrIdx);
S.GenomicPosition = cghdata.GenomicPosition(chrIdx);
end

function plotFrequency(S, appdata)
% Plot the frequency
cutoff = appdata.cutoff;
delta = 0.1;
initYL = [appdata.ylim(1)*(1+delta), appdata.ylim(2)*(1+delta)];

%== Chromosome locations
chromEnd = cumsum(appdata.chrominfo.chromLen(appdata.dispChrom));
chromStart = [0; chromEnd(1:end-1)];
%== Centromere locations
centrLoc = appdata.chrominfo.centrLoc(appdata.dispChrom)+ chromStart;

%== Check the genomic position unit in data, and transfer them to kb unit.
xLoc = double(S.GenomicPosition);
if max(xLoc) >= max(appdata.chrominfo.chromLen)
    xLoc = double(S.GenomicPosition/1000);
elseif max(xLoc) <= max(appdata.chrominfo.chromLen)/10^6
    xLoc = double(S.GenomicPosition*1000);
end

%== X point to be drawn
for i = 1:length(appdata.dispChrom)
   idx = S.Chromosome == appdata.dispChrom(i);
   xLoc(idx) = xLoc(idx) + chromStart(i);
end

% Figure
hfig = [];
if appdata.subplotFlag
    hfig = figure;
end

%== Groups and subplots
grpNum = numel(S.Group);
titleCnt = 0;
if ~isempty(appdata.titles)
    titleCnt = numel(appdata.titles);
end
    
for g = 1:grpNum
    %== Handle figure and axes
    if appdata.subplotFlag
        hAxes = subplot(grpNum, 1, g);
    else
        hfig = figure;
        hAxes = gca;
        set(hAxes, 'Parent', hfig);
    end

    if isempty(appdata.cscheme)
        hfigColorMap = rygCMap(64);
    else
        hfigColorMap = twoColorCMap(128, appdata.cscheme);
    end

    set(hfig,...
        'Visible','off',...
        'Colormap',hfigColorMap);

    set(hAxes,...
        'ALimMode', 'manual',...
        'XLimMode', 'manual',...
        'YLimMode', 'manual',...
        'ZLimMode', 'manual',...
        'CLimMode', 'manual',...
        'XTick', [],...
        'ALim',[0 1],...
        'XLim', [0 chromEnd(end)],...
        'YLim', initYL,...
        'ZLim', [0 1],...
        'CLim',[-1 1],...
        'Tag', ['Group' num2str(g)]);
    xlabel(hAxes, 'Chromosome number');
    ylabel(hAxes, 'Fraction of Samples');
    
    %== Filer the data with cutoff and nonzero for display
    gainData = S.Group(g).GainFrequency(S.Group(g).GainFrequency >= cutoff(1)&...
                                        S.Group(g).GainFrequency ~=0 )';
    lossData = S.Group(g).LossFrequency(S.Group(g).LossFrequency <= -cutoff(2)&...
                                        S.Group(g).LossFrequency ~=0)';

    xGain = repmat(xLoc(S.Group(g).GainFrequency >= cutoff(1)&...
                        S.Group(g).GainFrequency ~=0)', 2, 1);
    yGain = [gainData; zeros(size(gainData))];
    
    xLoss = repmat(xLoc(S.Group(g).LossFrequency <= -cutoff(2)&...
                        S.Group(g).LossFrequency ~=0)', 2, 1);    
    yLoss = [zeros(size(lossData)); lossData];

    lgColor = [0 1 0; 1 0 0];
    hg = hggroup('Parent',hAxes);
    if isempty(appdata.cscheme)
        hGainLines = line(xGain, yGain,...
                          'Tag', 'Gain',...
                          'Parent',hg);
        hLossLines = line(xLoss, yLoss,...
                          'Tag', 'Loss',...
                          'Parent',hg);

        for i = 1:length(hGainLines)
            set(hGainLines(i), 'Color', greenYellowRed(gainData(i)));
        end

        for i = 1:length(hLossLines)
            set(hLossLines(i), 'Color', greenYellowRed(lossData(i)));
        end
    else
        hGainLines = line(xGain, yGain, 'Color', appdata.cscheme(1,:),...
                          'Tag', 'Gain',... 
                          'Parent',hg); %#ok
        hLossLines = line(xLoss, yLoss, 'Color', appdata.cscheme(2,:),...
                          'Tag', 'Loss',...
                          'Parent',hg); %#ok
        lgColor = appdata.cscheme;
        setLegendDispStyle(hg, 'off')
    end

    % Lines for legend
    lgProp = {'Visible', 'off', 'HitTest', 'off'};
    line([0 0], [0 0], 'Color', lgColor(1,:), 'DisplayName', 'Gain', lgProp{:});
    line([0 0], [0 0], 'Color', lgColor(2,:), 'DisplayName', 'Loss', lgProp{:});

    %== Draw Chromosome border lines
    hChrom = [];
    if length(appdata.dispChrom) > 1 % if only one chromosome, no border
        xChrom = repmat(chromEnd(1:end)', 2, 1);
        yChrom = repmat(initYL(:), 1, length(appdata.dispChrom));
        hChrom = line(xChrom, yChrom, 'color', [0 0 0], 'Tag', 'Chrom');
        setLegendGroup(hAxes, hChrom, 'on', 'Chromosomes');
    end
    %== Draw centromere
    hCentr = [];
    if appdata.dispCentrFlag
        xCentr = repmat(centrLoc', 2, 1);
        yCentr = repmat(initYL', 1, length(appdata.dispChrom));
        hCentr = line(xCentr, yCentr, 'color', [0.7 0.7 0.7],...
                     'LineStyle', ':', 'Tag', 'Centr');
        %set legend
        setLegendGroup(hAxes, hCentr, 'on', 'Centromere');
    end
    
    %== Chromsome labels
    % Label the autosomes with their chromosome numbers, and the sex
    % chromosome with X or Y.
    xLabel = chromStart + appdata.chrominfo.chromLen(appdata.dispChrom)/2;
    yLabel = zeros(1, length(xLabel)) + initYL(1);
    chromLabel = cellstr(num2str(appdata.dispChrom(:)));
    if any(appdata.dispChrom == appdata.chrominfo.XNum)
        chromLabel{appdata.dispChrom == appdata.chrominfo.XNum} = 'X';
    end
    if any(appdata.dispChrom == appdata.chrominfo.YNum)
        chromLabel{appdata.dispChrom == appdata.chrominfo.YNum} = 'Y';
    end

    hChrLabel = text(xLabel, yLabel, chromLabel,...
                    'Fontsize', 8,...
                    'HorizontalAlignment', 'center',...
                    'Clipping', 'on',...
                    'VerticalAlignment', 'Bottom',...
                    'FontWeight', 'demi');

    %== Draw a zero baseline
    hZero = line([0 chromEnd(end)], [0 0], 'Color', [0 0 0], 'Hittest', 'off');
    setLegendDispStyle(hZero, 'off')

    %== Handle data curor
    dcm_obj = datacursormode(hfig);
    set(dcm_obj,'UpdateFcn',{@dataCursorUpdatefcn, appdata, chromEnd, centrLoc,...
        gainData, xGain(1,:), lossData, xLoss(1,:)})
    
    %== Handle post zoom behavoir
    zoom_obj = zoom(hfig);
    set(zoom_obj, 'ActionPostCallback', {@zoomPostUpdatefcn, hChrLabel})
    
    %== Rotate 3d off
    rot_obj = rotate3d(hfig);
    setAllowAxesRotate(rot_obj,hAxes,false)
     %== Handle post pan behavoir
    pan_obj = pan(hfig);
    set(pan_obj, 'ActionPostCallback', {@panPostUpdatefcn, hChrLabel, hChrom, hCentr})
    
    %== Colorbar object (axes)
    cbarbutton = uigettool(gcf, 'Annotation.InsertColorbar');
    set(cbarbutton, 'ClickedCallback', {@updateColorbarCB, hfig, hAxes}) 
    %== Titles    
    if titleCnt > 0
        title(appdata.titles{g}, 'Interpreter', 'none');
        titleCnt = titleCnt - 1;
    end

    %== Set lines visible
    set(hAxes, 'Box', 'on')
    set(hfig, 'Visible', 'on')

end
end

function updateColorbarCB(hsrc, event, hfig, cax) %#ok
% Insert Colorbar callback function to modify the color bar context menu
% items.
hfig = gcbf;
cbmenu = findall(hfig,'Tag','figMenuInsertColorbar');
cbtogg = uigettool(hfig,'Annotation.InsertColorbar');
cax = get(hfig,'CurrentAxes');
if ~isempty(cax)
    if (~isempty(cbmenu) && isequal(cbmenu,gcbo) && strcmpi(get(cbmenu,'checked'),'off')) ||...
            (~isempty(cbtogg) && isequal(cbtogg,gcbo) && strcmpi(get(cbtogg,'state'),'on'))
        cbar = colorbar('peer',cax);
        updateColorbarContextmenu(cbar)
    else
        colorbar('peer',cax,'off');
    end
else
    if ~isempty(cbtogg)
        set(cbtogg,'State','off');
    end
end
end

function updateColorbarContextmenu(cbar)
% Set colorbar context menu items "Visible' off.
uic = get(cbar, 'UIContextmenu');
set(findall(uic,'Type','uimenu','Tag','scribe:colorbar:colormap'), 'Visible', 'off');
set(findall(uic,'Type','uimenu','Tag','scribe:colorbar:interactivecolormapshift'), 'Visible', 'off');
set(findall(uic,'Type','uimenu','Tag','scribe:colorbar:editcolormap'), 'Visible', 'off');
set(findall(uic,'Type','uimenu','Tag','scribe:colorbar:propedit'), 'Visible', 'off');
set(findall(uic,'Type','uimenu','Tag','scribe:colorbar:mcode'), 'Visible', 'off');
end

function hg = setLegendGroup(hAxes, hObjs, iconDispStyle, dispName)
% Add hObjs to a hggroup and set the legend for this group
% iconDispStyle - 'on' show legend, 'off' not show.
% Label - legend display string
hg = hggroup('Parent',hAxes);
set(hObjs, 'Parent', hg);
set(get(get(hg,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle',iconDispStyle); % Include this hggroup in the legend
set(hg, 'DisplayName', dispName);
end

function setLegendDispStyle(hObj, iconDispStyle)
% Set hObj legend displaystyle
% iconDispStyle - 'on' show legend, 'off' not show.
set(get(get(hObj,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle',iconDispStyle); % 
end

function datacursorLabel = dataCursorUpdatefcn(obj, event_obj, appdata,...
                           chromEnd, centrLoc,...
                           gainData, gainLoc, lossData, lossLoc) %#ok<INUSL>
% Display the data cursor.
datacursorLabel = [];
tg = get(event_obj, 'Target');
if strcmpi(get(tg, 'Type'), 'line')
    type = get(tg, 'Tag');
    if strcmpi(type, 'Chrom')
        pos = get(event_obj, 'Position');
        chromn = appdata.dispChrom(pos(1) == chromEnd);
        if ~isempty(chromn)
            chromlen = appdata.chrominfo.chromLen(chromn);
            datacursorLabel = {['Chromosome ' num2str(chromn)],...
                ['Length: ' num2str(chromlen) ' kb']};
        end
    elseif strcmpi(type, 'Centr')
        pos = get(event_obj, 'Position');
        chromn = appdata.dispChrom(pos(1) == centrLoc);
        if ~isempty(chromn)
            centrloc = appdata.chrominfo.centrLoc(chromn);
            datacursorLabel = {['Chromosome ' num2str(chromn)],...
                ['Centromere: ' num2str(centrloc) ' kb']};
        end
    elseif strcmpi(type, 'Gain')
        pos = get(event_obj, 'Position');
        k = find(pos(1) == gainLoc, 1);
        chromn = appdata.dispChrom(find(pos(1) <= chromEnd, 1, 'first'));
        if ~isempty(k)
            datacursorLabel = {['Gain: ' num2str(gainData(k))],...
                ['Location: ' num2str(pos(1)) ' kb on Chr ' num2str(chromn)]};
        end
     elseif strcmpi(type, 'Loss')
        pos = get(event_obj, 'Position');
        k = find(pos(1) == lossLoc, 1);
        chromn = appdata.dispChrom(find(pos(1) <= chromEnd, 1, 'first'));
        if ~isempty(k)
            datacursorLabel = {['Loss: ' num2str(lossData(k))],...
                ['Location: ' num2str(pos(1)) '(kb) on Chr ' num2str(chromn)]};
        end
    end
end
end

function zoomPostUpdatefcn(obj,event_obj, hLabel) %#ok
% Update the chromosome label positions post zoom actions
yLim = get(event_obj.Axes,'YLim');

% Update Label position
for i = 1:length(hLabel)
   pos = get(hLabel(i), 'Position');
   
   pos(2) = yLim(1);
   set(hLabel(i), 'Position', pos)
end
end

function panPostUpdatefcn(obj,event_obj, hLabel, hChrom, hCentr) %#ok
% After pan the chromosome and centromere lines should meet the limits
yLim = get(event_obj.Axes,'YLim');

% Update Label position
for i = 1:length(hLabel)
   pos = get(hLabel(i), 'Position');
   
   pos(2) = yLim(1);
   set(hLabel(i), 'Position', pos)
end
% Update chromsome line length
for i = 1:length(hChrom)
    set(hChrom(i), 'YData', yLim)
end

if ~isempty(hCentr)
    for i = 1:length(hCentr)
        set(hCentr(i), 'YData', yLim)
    end
end

end

function chrominfo = getChromosomeInfo(cytoStruct)
% Return a structure with the length of chromosomes in kb unit, centromere
% positions and autosome chromosome numbers, X and Y numbers

% Get chromsome length
[chromLabels, endIdx] = unique(cytoStruct.ChromLabels);

chromLength = cytoStruct.BandEndBPs(endIdx)/1000;
% Find the centromere positions for the chromosomes.
acenIdx = strcmpi(cytoStruct.GieStains, 'acen');
acenEnds = cytoStruct.BandEndBPs(acenIdx);

% Convert the cytoband data from bp to kilo bp.
centrLoc = acenEnds(1:2:end)/1000;

chrs = str2double(chromLabels);
chrs(strcmpi('X', chromLabels)) = max(chrs)+1;
chrs(strcmpi('Y', chromLabels)) = max(chrs)+1;

[chromNumbers, I] = unique(chrs);

XNum = 0;
YNum = 0;
if any(strncmp(chromLabels, 'Y', 1));
    YNum = chromNumbers(end);
end

if any(strncmp(chromLabels, 'X', 1));
    if YNum > 0
        XNum = YNum -1;
    else
        XNum = chromNumbers(end);
    end
end
chrominfo.autosomeNum = chromNumbers(chromNumbers ~= XNum & chromNumbers ~= YNum);
chrominfo.XNum = XNum;
chrominfo.YNum = YNum;
chrominfo.centrLoc = double(centrLoc);
chrominfo.chromLen = double(chromLength(I));
end

%=== Color scheme functions
function h = greenYellowRed(m)
%greenYellowRed generate green-yellow-red color RGB 1x3 vector
%   greenYellowRed(M) returns an 1-by-3 vector containing a green-yellow-red RGB
%   values based on value M. -1<M<1. red(-1) - yellow(0) - green(1) 

if m <= 0
    h=[1 m+1 0];
else
    h =[1-m 1 0];
end

end

function p = rygCMap(m, varargin)
%RYGCMAP creates a red-yellow-green colormap for colorbar.
% RYGCMAP(M) returns an M-by-3 matrix containing a red and green colormap.
% Low values are bright red, values in the center of the map are yellow,
% and high values are gree. If M is empty, the length of the map will be
% the same as the length of the colormap of the current figure.

if nargin < 1 || isempty(m) || m < 3
    m = size(get(gcf,'colormap'),1);
end

coloredLength = floor((m-1)/2);

% create an appropriately long linearly interpolated chunk
interpCol = (1/(coloredLength):1/(coloredLength):1);

if coloredLength == ((m-1)/2)
    fillerZeros = 1;
else
    fillerZeros = [1 1];
end

% red is linear for red half
red = [ones(size(interpCol)) fillerZeros 1-interpCol];
% green is opposite of red
green = [interpCol fillerZeros ones(size(interpCol))];
% blue is lowintensity for red half
blue = zeros(size(red));

p = [red',green',blue'];
end

function p = twoColorCMap(m, cmap)
%TWOCOLORCMAP creates a two-color colormap for colorbar.
% TWOCOLORCMAP(M) returns an M-by-3 matrix containing two color colormap.
% If M is empty, the length of the map will be the same as the length of
% the colormap of the current figure.

if isempty(m) || m < 3
    m = size(get(gcf,'colormap'),1);
end

hLength = floor((m-1)/2);
fillerZeros = [0 0 0];

p = [repmat(cmap(2,:), hLength, 1);fillerZeros;repmat(cmap(1,:), hLength, 1)];
end

%=== Input checking functions
function [cghdata, instruct]= check_inputdata(data, instruct)
% Check input data, and validate the input parameters

% Persistent Homo sapiens chromosome information
persistent hsChromInfo;

%=== Check data
if isstruct(data)
    datatype = checkInputStruct(data); 
    if datatype < 0
        error(message('bioinfo:cghfreqplot:InvalidDataFieldNames'))
    end
    % Input structure contain frequency, only need to plot them
    plotOnly = datatype == 1;
    cghdata = data;    
% matrix with one sample data	
elseif isnumeric(data) && isreal(data) && size(data, 2) >= 3
    cghdata = struct('Sample', [], 'Chromosome', data(:, 1),...
                     'GenomicPosition', data(:, 2),...
                     'Log2Ratio', data(:, 3:end)); 
    plotOnly = false;
else
    error(message('bioinfo:cghfreqplot:WrongInputDataFormat'))
end

if ~plotOnly
    [~, nsamples] = size(cghdata.Log2Ratio);

    % More than one sample, should have Sample field
    if isfield(cghdata,'Sample') && ~isempty(cghdata.Sample) && nsamples ~= numel(cghdata.Sample)
       error(message('bioinfo:cghfreqplot:SampleDataSizeNotMatch')) 
    end

    %=== Add sample names if there are not sample names
    if ~isfield(cghdata,'Sample') || isempty(cghdata.Sample)
        sid = num2str((1:nsamples)');
        cghdata.Sample = strcat({'Sample'},cellstr(sid))';
    end

%=== Validate the input parameters with DATA
% Threshold
ntmp = size(instruct.threshold(:), 1);
if ntmp > 2 && ntmp ~= nsamples
    error(message('bioinfo:cghfreqplot:ThresholdSizeNotMatchSampleSize'))
end

% Group
if isempty(instruct.group)
    instruct.group = {1:nsamples};
elseif ~iscell((instruct.group))
    group = instruct.group; 
    instruct.group = {getSampleMemberID(group, 1:nsamples)};
elseif iscell(instruct.group) % cell array
    for i = 1:numel(instruct.group)
       group = instruct.group{i};
       instruct.group{i} = getSampleMemberID(group, 1:nsamples);
    end
end
else
    ngrp = numel(cghdata.Group);
    nsamples = 0;
    for i = 1:ngrp
        nsamples = nsamples + numel(cghdata.Group(i).Sample);
    end
end % if ~plot

% Chrominfo Get default cytoband information
if isempty(instruct.chrominfo)
    if isempty(hsChromInfo)
        % get Homo sapien's cytoband information
        instruct.chrominfo = cytobandread('hs_cytoBand.txt');
        hsChromInfo = getChromosomeInfo(instruct.chrominfo);
    end
    chrominfo = hsChromInfo;
else
    chrominfo = getChromosomeInfo(instruct.chrominfo);
end

% Chromosome
chromIDs = unique(cghdata.Chromosome);
autosomeIDs = chromIDs(ismember(chromIDs, chrominfo.autosomeNum));
if instruct.dispChrom == 0
    instruct.dispChrom = autosomeIDs;
    if instruct.XFlag
        if any(chromIDs == chrominfo.XNum)
            instruct.dispChrom = [instruct.dispChrom;chrominfo.XNum];
        else
            warning(message('bioinfo:cghfreqplot:XNotInData'));
        end
    end
    
    if instruct.YFlag
        if any(chromIDs == chrominfo.YNum)
            instruct.dispChrom = [instruct.dispChrom; chrominfo.YNum];
        else
           warning(message('bioinfo:cghfreqplot:YNotInData')); 
        end
    end
else
    instruct.dispChrom = getChromosomeMemberID(instruct.dispChrom, chromIDs);
end

if ~isempty(instruct.titles)
    titles = instruct.titles;
    if ischar(titles)
        titles = {titles};
    elseif isnumeric(titles)
        titles = {titles(:)};
    end
    instruct.titles = titles;
end

instruct.chrominfo = chrominfo;
instruct.nsamples = nsamples;
instruct.plotonlyFlag = plotOnly;
end

function valid = checkInputStruct(dataStruct)
% Checks for input data structure field names

valid = isfield(dataStruct,'Chromosome')&& ...
    isfield(dataStruct,'GenomicPosition') && ...
    isfield(dataStruct,'Log2Ratio');
if ~valid
    valid = isfield(dataStruct,'Group')&& ...
        isfield(dataStruct,'Chromosome')&& ...
        isfield(dataStruct,'GenomicPosition');
    if valid
        valid = isfield(dataStruct.Group,'Sample') && ...
                isfield(dataStruct.Group,'GainFrequency') && ...
                isfield(dataStruct.Group,'LossFrequency');
    end
    if valid
        valid = 1;
    else
        valid = -1;
    end
    return;
end

valid = 0;
end % end of function


function inputStruct = parse_inputs(varargin)
% Parse input PV pairs.

% Check for the right number of inputs
if rem(nargin,2)== 1
     error(message('bioinfo:cghfreqplot:IncorrectNumberOfArguments', mfilename))
end

% Allowed inputs
okargs = {'threshold', 'group', 'subgrp', 'subplot',...
          'cutoff', 'chromosome','includex', 'includey',...
          'chrominfo', 'showcentr', 'color', 'ylim', 'titles'};

% Defaults
inputStruct.dispChrom = 0;        % ChromosomeIDs
inputStruct.threshold = 0.25;   % threshold of Log2Ratio
inputStruct.cutoff = [0 0];     % don't show freq below cutoff
inputStruct.group = [];          % sample groups
inputStruct.subgrpFlag = true;  % analyze by subgroup 
inputStruct.subplotFlag = true; % plot subgroup in subplots
inputStruct.XFlag = true;  % include X chromosome
inputStruct.YFlag = false; % include Y chromosome
inputStruct.dispCentrFlag = true; % show centromeres on chromosomes
inputStruct.chrominfo = []; % cytoband information file
inputStruct.cscheme = []; % Red-yellow-green 
inputStruct.ylim = [-1 1];
inputStruct.titles = [];

for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1 % threshold
            if isvector(pval) && isnumeric(pval) && isreal(pval)
                inputStruct.threshold = abs(pval(:)');
            else
                error(message('bioinfo:cghfreqplot:ThresholdOptionNotRealNumericVector'))
            end
        case 2 % group
            if iscell(pval)
                res1 = cellfun('isreal', pval);
                res2 = cellfun(@isnumeric,pval);
                res3 = cellfun(@isvector,pval);
                tf = all(res1(:)) && all(res2(:)) && all(res3(:));
                if tf 
                    inputStruct.group = pval;
                else
                    error(message('bioinfo:cghfreqplot:GroupCellOptionNotRealNumericVector'))
                end
            elseif isvector(pval) && isnumeric(pval) && isreal(pval)
                inputStruct.group = pval;
            else
                error(message('bioinfo:cghfreqplot:GroupOptionNotRealNumericVector'))
            end
        case 3 % subgroup
            inputStruct.subgrpFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 4 % subplot
            inputStruct.subplotFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 5 % cutoff
            if isvector(pval) && isnumeric(pval) && isreal(pval) && (size(pval(:), 1) <= 2)
                cutoff = abs(pval);
                if ~all(cutoff <= 1)
                    error(message('bioinfo:cghfreqplot:InvalidCutoffValue'))
                end

                if numel(cutoff) < 2
                    inputStruct.cutoff = ones(1,2) * cutoff;
                end
            else
                error(message('bioinfo:cghfreqplot:InvalidCutoffFormat'))
            end
        case 6  % chromosome
            if isvector(pval) && isnumeric(pval) && isreal(pval)
                inputStruct.dispChrom = fix(pval(:));
            else
                error(message('bioinfo:cghfreqplot:ChromosomeOptionNotRealNumericVector'))
            end            
        case 7 % includeX
            inputStruct.XFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 8 % includeY
            inputStruct.YFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 9 % chromosomeinfo
            if ~isstruct(pval)
                try
                    inputStruct.chrominfo = cytobandread(pval);
                catch theException
                    msgId = 'bioinfo:cghfreqplot:UnableReadInputFile';
                    newException = MException(msgId,getString(message(msgId)));
                    throw(addCause(newException,theException))
                end
            else
                if checkCytoBandStruct(pval)
                    inputStruct.chrominfo = pval;
                else
                    error(message('bioinfo:cghfreqplot:BadFieldNames'))
                end
            end
        case 10 % showacent
            inputStruct.dispCentrFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 11 % color
            validscheme = true;
            if ischar(pval) || isa(pval,'function_handle')
                try
                    colorscheme = eval(pval);
                catch allExceptions %#ok<NASGU>
                    validscheme = false;
                    colorscheme = [];
                end
            elseif isnumeric(pval)
                colorscheme = pval;
            else
                validscheme = false;
            end

            if ~validscheme || ~ismatrix(colorscheme) ||...
                  size(colorscheme,2) ~= 3
              error(message('bioinfo:cghfreqplot:InvalidColorFormat'))
            end
            if size(colorscheme, 1) > 1
                inputStruct.cscheme = colorscheme(1:2, :);
            else
                inputStruct.cscheme = repmat(colorscheme, 2, 1);
            end
        case 12 % ylim
            if isvector(pval) && isnumeric(pval) && isreal(pval) && (size(pval(:), 1)==2)
                inputStruct.ylim = pval(:)';
            else
                error(message('bioinfo:cghfreqplot:InvalidYlimFormat'))
            end
        case 13 % titles
            inputStruct.titles = pval;
    end
end
end

function valid = checkCytoBandStruct(cbStruct)
% Checks for cytoband structure field names

valid =  isfield(cbStruct,'ChromLabels')&& ...
         isfield(cbStruct,'BandStartBPs') && ...
         isfield(cbStruct,'BandEndBPs') && ...
         isfield(cbStruct,'BandLabels') && ...
         isfield(cbStruct,'GieStains');
end % checkCytoBandStruct function

function  memID = getSampleMemberID(memID, allID)
    % Check for member ids from a set of IDs, if any id in memID not in
    % allID, throw a warning with type in warning ID.
    memID = unique(memID(:))';
    idx = (ismember(memID, allID));
    if ~all(idx)
        warning(message('bioinfo:cghfreqplot:InvalidSampleId', num2str( memID( ~idx ) )));
        memID = memID(idx);
    end
end

function  memID = getChromosomeMemberID(memID, allID)
    % Check for member ids from a set of IDs, if any id in memID not in
    % allID, throw a warning with type in warning ID.
    memID = unique(memID(:))';
    idx = (ismember(memID, allID));
    if ~all(idx)
        warning(message('bioinfo:cghfreqplot:InvalidChromosomeId', num2str( memID( ~idx ) )));
        memID = memID(idx);
    end
end