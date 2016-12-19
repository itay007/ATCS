function [diffStruct, hFig] = mavolcanoplot(X, Y, PV, varargin)
%MAVOLCANOPLOT creates a significance vs expression ratio scatter plot.
%
%   MAVOLCANOPLOT(X,Y,PV) creates a significance vs expression ratio
%   scatter plot. X and Y are the log2 gene expression values. X or Y can
%   be a vector, a matrix with column-wise values or a DataMatrix object.
%   PV can be a vector of or a DataMatrix object containing p-values
%   resulting from statistical tests, e.g., the Student's t tests, of X and
%   Y. 
%   
%   Notes: 1) The statistical difference is measured by the negative log10
%   of the p-value for the corresponding statistical t-test of differences
%   between X and Y.  The biological effect is the log2 of the fold change
%   between the two groups. 2) The genes with both statistical significance
%   and biological significance will be located in the upper left and upper
%   right of the volcano plot.   
% 
%   MAVOLCANOPLOT(...,'LOGTRANS', true) will log2 transform X and Y from
%   linear scale. Otherwise X and Y are assumed to be log2 transformed.
%
%   MAVOLCANOPLOT(...,'PCUTOFF', P) adds a line to indicate the cutoff
%   p-value P. By default P=0.05. The cutoff line can be changed
%   interactively in the plot.
%
%   MAVOLCANOPLOT(...,'FOLDCHANGE', N) adds lines showing N fold-change. By
%   default N=2. The fold-change lines can be changed interactively in the
%   plot.
% 
%   MAVOLCANOPLOT(...,'LABELS',LABELS) allows you to specify a cell array
%   of gene labels or names for each row of the input data. If LABELS are
%   defined, then clicking on a data point on the plot will show the LABEL
%   corresponding to that data point. If LABELS are not defined at input,
%   the row indices will be used as LABELS.
% 
%   Note: If the input data contains DataMatrix objects, and 'LABELS'
%   option is not defined. The row names of the first input DataMatrix
%   object are used as LABELS. 
% 
%   MAVOLCANOPLOT(...,'PLOTONLY',TR) only displays the plot without the UI
%   components when TF is TRUE. The default is FALSE.
%
%   DIFFSTRUCT = MAVOLCANOPLOT(...) mavolcanoplot(...) returns a structure
%   containing information for differentially expressed genes. 
%
%   Examples:
%           load prostatecancerexpdata
%           pvals=mattest(dependentData, independentData);
%           diff=mavolcanoplot(dependentData,independentData,pvals, 'Labels', probesetIDs);
%           
%           % View volcano plot without UI components
%           mavolcanoplot(dependentData,independentData,pvals,...
%                         'Labels', probesetIDs, 'Plotonly', true);
%
%   See also CNSGENEEXPDEMO, MABOXPLOT, MAIRPLOT, MALOGLOG, MAPCAPLOT, MATTEST.

%   MAVOLCANOPLOT(...,'COLORS', COLORS) allows user to set the color scheme
%   for the volcano plot. Options are 'blue' (default), 'cyan', 'green',
%   'orange', 'pink', or 'redgrey'. If COLORS is a matrix with exactly two
%   rows, three columns, and elements in the interval 0 to 1, The first row
%   is the RGB color for the significant genes and the second row is the
%   RGB color insignificant ones.  The columns of the input matrix
%   represent intensity of red, blue and green, respectively.
% 
%   [DIFFSTRUCT, H] = MAVOLCANOPLOT(...) returns the figure handle.

% Copyright 2003-2012 The MathWorks, Inc.


% References: 
% [1] X. Cui, G.A.Churchill. Statistical tests for differential expression
%     in cDNA microarray experiments. Genome Biol. 2003;4(4):210. Epub 2003
%     Mar 17. Review. PMID: 12702200.

bioinfochecknargin(nargin,3,mfilename);

%== Check input parameters
paramStruct = parse_inputs(varargin{:});

dataMatrixFlag = false;

if isa(X, 'bioma.data.DataMatrix')
    if isempty(paramStruct.labels)
        paramStruct.labels = X.RowNames;
    end
    X = X.(':')(':');
    dataMatrixFlag = true;
end

if isa(Y, 'bioma.data.DataMatrix')
    if isempty(paramStruct.labels)
        paramStruct.labels = Y.RowNames;
    end
    Y = Y.(':')(':');
    dataMatrixFlag = true;
end

if isa(PV, 'bioma.data.DataMatrix')
    if isempty(paramStruct.labels)
        paramStruct.labels = PV.RowNames;
    end
    PV = PV.(':')(':');
    dataMatrixFlag = true;
end

%== Check input data and validate parameters
[X, Y, PV, paramStruct] = check_inputdata(X, Y, PV, paramStruct);

%== Do plot
figName = 'Volcano Plot';
units = 'Normalized';

hFig = figure('Units',units,...
              'Tag', 'volcanofig',...
              'Visible', 'on');

%== Set appdata
appdata = localGetAppData(hFig);
appdata.hfig = hFig;
appdata.plotonly = paramStruct.plotOnly;
appdata.glnumericFlag = paramStruct.glnumericFlag;
appdata = setColorScheme(appdata, paramStruct.colorScheme,...
                                  paramStruct.sigRGB, paramStruct.insRGB);
appdata.pvalues = PV(paramStruct.goodVals);
appdata.labels = paramStruct.labels(paramStruct.goodVals);
appdata.significance = -log10(PV(paramStruct.goodVals));
appdata.effect = X(paramStruct.goodVals) - Y(paramStruct.goodVals);

hwmsg = [];
if log2(abs(paramStruct.fcThd)) > max(abs(appdata.effect)) ||...
        log2(abs(paramStruct.fcThd)) < min(abs(appdata.effect))
    msg = sprintf(['%d cannot be displayed in the plot for it is out of range.'...
        '\nA default two-fold change is used.'], paramStruct.fcThd);
    hwmsg = warndlg(msg, 'Volcano Plot');
    paramStruct.fcThd = 2;
end

appdata.xlim = [min(appdata.effect)-0.2, max(appdata.effect)+0.2];
appdata.ylim = [0, max(appdata.significance)];
appdata.upMax = numel(appdata.effect);
appdata.dnMax = numel(appdata.effect);
appdata.DataMatrixFlag = dataMatrixFlag;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the plot
appdata.xThd = paramStruct.fcThd;
appdata.yThd = paramStruct.pvCutoff;
appdata.newx = log2(appdata.xThd);
appdata.newy = -log10(appdata.yThd);
appdata.ltype = 0;

%********* Layout UI *************
if ~appdata.plotonly
    fw = 0.475; % the figure width in relation to the screen width
    fh = 0.5; % the figure width in relation to the screen width
    if strncmpi(computer, 'GLNX', 4)
        fw = 0.75;
        fh = 0.65;
    end
   set(hFig, 'Name', figName,...
        'Numbertitle','off', ...
        'PaperUnits','normalized', ...
        'PaperPosition',[.55 .30 .5 .5], ...
        'Position',[.15 .25 fw fh]);
    % Reset figure menu items
    resetFigureTools(hFig);
    appdata = createUIs(hFig, appdata);
else
    dcm_obj = datacursormode(hFig);
    set(dcm_obj,'UpdateFcn',{@locaDataCursorUpdate, hFig})
end
    
%************* Create interactive scatter plot********
appdata = createPlot(hFig, appdata);
appdata = updateLists(appdata);

% Save state
localSetAppData(hFig,appdata);

set(hFig,'WindowButtonMotionFcn',{@movelines, 'motion',0},...
         'Visible','on');
localSetAppData(hFig,appdata);

if ishandle(hwmsg)
    figure(hwmsg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if output is requested, send out handle
if nargout > 0   
    idx = getDiffIndexes(appdata, 0);
    diffPV = appdata.pvalues(idx);

    if appdata.glnumericFlag
        diffL = conv2numric(appdata.labels(idx));
    else
        diffL = appdata.labels(idx);
    end
    
    diffFC = appdata.effect(idx);
    
    [diffPV, idx] = sort(diffPV);
    diffFC  = diffFC(idx);
    
    diffStruct.Name = 'Differentially Expressed';
    diffStruct.PVCutoff = appdata.yThd;
    diffStruct.FCThreshold = appdata.xThd;
    diffStruct.GeneLabels = diffL(idx);
    diffStruct.PValues = diffPV;
    diffStruct.FoldChanges = 2 .^ (abs(diffFC)); 
    
    idx = diffFC < 0;
    diffStruct.FoldChanges(idx) = - diffStruct.FoldChanges(idx);
    if appdata.DataMatrixFlag 
        diffStruct.PValues = bioma.data.DataMatrix(diffStruct.PValues,...
                                diffStruct.GeneLabels, {'p-values'});
        diffStruct.FoldChanges = bioma.data.DataMatrix(diffStruct.FoldChanges,...
                                diffStruct.GeneLabels, {'FoldChanges'});
    end
end
end

%----------------------------------------------
function txt = locaDataCursorUpdate(hobj,eobj, hfig) %#ok<*INUSL>
% Update the datacursor text.
tg = get(eobj, 'Target');
appdata = localGetAppData(hfig);

pos = get(eobj,'Position');
dataidx = get(eobj, 'DataIndex');

if strcmpi(get(tg, 'Type'), 'patch')    
    txt = {['Label: ', appdata.labels{dataidx}],...
        ['Fold-change: ',num2str(pos(1))],...
        ['p-value: ',num2str(10^(-pos(2)))]};
elseif strcmpi(get(tg, 'Type'), 'line')
    type = get(tg, 'Tag');
    if strcmpi(type, 'horzline')
         txt = {['PVCutoff: ', num2str(10^(-pos(2)))]};
    elseif strcmpi(type, 'vertline')
         txt = {['FoldChange: ', num2str(pos(1))]};
    else
        txt = [];
        return;
    end
else
    txt = [];
    return;
end
end
%-------------------------------------
function clickonplot(varargin)
% callback function highlights selected element and displays label.
[hObject,hfig] = gcbo;

% clean up any old labels and get appdata
appdata = clearSelections([], [], hfig);

if( hObject == appdata.hAxis)
    hObject = appdata.hPlot;
end

xvals = get(hObject,'XData');
yvals = get(hObject,'YData');
point = get(appdata.hAxis,'CurrentPoint');

% find the closest point
[~, index] = min((((xvals - point(1,1)).^2)./(point(1,1)^2)) +...
    ((yvals - point(1,2)).^2)./(point(1,2)^2));

%highlight the point and set the label
if ~isempty(index)
    hHighlight = line(xvals(index),yvals(index),...
                        'color',appdata.markerRGB,...
                        'marker','o',...
                        'linewidth', 1.5,...
                        'Tag','volHighlight');%#ok

    htext = text(point(1,1),point(1,2) , appdata.labels(index));
            
    % give it an off white background, black text and grey border
    set(htext,	'HorizontalAlignment','left',...
                'VerticalAlignment', 'bottom',... 
                'interpreter','none',...
                'Tag','volDataTip',...
                'Visible','on');
            
    %Check extends
    setlabelposition(htext, appdata.hAxis)
    
    set(hfig,'WindowButtonUpFcn',@clearClickLabel);
    localSetAppData(hfig,appdata);
end
end
%%%%%%%%%%%%%%%%%% Callback functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hFig = clearClickLabel(varargin) 
%clearClickLabel callback function to remove label from plot
hFig = gcbf;
% delete the old label if it exists
oldLabel = findobj(hFig,'Tag','volDataTip');
if ~isempty(oldLabel)
    delete(oldLabel);
end
% delete the old label if it exists
oldHighlight = findobj(hFig,'Tag','volHighlight');
if ~isempty(oldHighlight)
    delete(oldHighlight);
end
end

function appdata = clearSelections(hsrc, evet, hfig)
appdata = localGetAppData(hfig);

appdata.upSelectedText = deleteSelectedHandle(appdata.upSelectedText);
appdata.upSelectedMarkers = deleteSelectedHandle(appdata.upSelectedMarkers);
appdata.dnSelectedText = deleteSelectedHandle(appdata.dnSelectedText);
appdata.dnSelectedMarkers = deleteSelectedHandle(appdata.dnSelectedMarkers);
set(appdata.upllist, 'Value', []);
set(appdata.upplist, 'Value', []);
set(appdata.dnllist, 'Value', []);
set(appdata.dnplist, 'Value', []);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function movelines(varargin)
[~,hfig] = gcbo;

appdata = localGetAppData(hfig);
if ishandle(appdata.hAxis)
    axes  = appdata.hAxis;
else
    return;
end

cp = get(axes, 'CurrentPoint');
action = 'wait';

if nargin > 2
    action = varargin{3};
end

if nargin > 3
    appdata.ltype = varargin{4};
end
xrange = get(axes,'Xlim');
yrange = get(axes,'Ylim');
if strcmp(action,'down'),
    if ~isinaxes(cp, axes)
         set(hfig, 'Pointer', 'arrow');
        return;
    end
    clearSelections([], [], hfig);
    % Set cursor on the figure
    set(hfig, 'Pointer', 'crosshair');
    set(hfig,'WindowButtonMotionFcn',{@movelines, 'motion', appdata.ltype });
    set(hfig,'WindowButtonUpFcn',{@movelines, 'up'});
elseif strcmp(action, 'motion')
    if ~isinaxes(cp, axes)
        set(hfig, 'Pointer', 'arrow');
        return;
    end
    if appdata.ltype == 0
        cursorstate = get(hfig,'Pointer');
        cp = get(axes,'CurrentPoint');
        cx = cp(1,1);
        cy = cp(1,2);
        fuzzx = 0.01 * (xrange(2) - xrange(1));
        fuzzy = 0.01 * (yrange(2) - yrange(1));
        online = cy > yrange(1) && cy < yrange(2) && ...
            cx > xrange(1) && cx < xrange(2) &&...
            ((cy > appdata.newy - fuzzy && cy < appdata.newy + fuzzy) ||...
            (cx > appdata.newx - fuzzx && cx < appdata.newx + fuzzx) ||...
            (cx > -appdata.newx - fuzzx && cx < -appdata.newx + fuzzx));
        
        if online && strcmp(cursorstate,'arrow'),
            set(hfig,'Pointer','crosshair');
        elseif ~online && strcmp(cursorstate,'crosshair'),
            set(hfig,'Pointer','arrow');
        end
 
    elseif appdata.ltype == 2 || appdata.ltype == 3
        set(hfig, 'Pointer', 'crosshair');
        % Check the newpoints are within range
        appdata.newx = cp(1,1);
        if appdata.newx > xrange(2)
            appdata.newx = xrange(2);
        end
        if appdata.newx < xrange(1);
            appdata.newx = xrange(1);
        end

        if appdata.ltype == 2 && appdata.newx <= 0
            appdata.newx = 0.05;
        end
    
        if appdata.ltype == 3 && appdata.newx >= 0
            appdata.newx = 0.05;
        end
        
        appdata.newx = abs(appdata.newx);
        
        %     get the xfield handle
        set(appdata.xField, 'String', num2str(2^(appdata.newx), '%0.5g'));
        set(appdata.xField, 'Userdata', 2^(appdata.newx));
        set(appdata.vline1,'XData', repmat(appdata.newx, size(appdata.inityrange)),'YData', appdata.inityrange);
        set(appdata.vline2,'XData', repmat(-appdata.newx, size(appdata.inityrange)),'YData', appdata.inityrange);

        localSetAppData(hfig,appdata);
    elseif appdata.ltype == 1
        set(hfig, 'Pointer', 'crosshair');
        % Check the newpoints are within range
        appdata.newy = cp(1,2);
        if appdata.newy > yrange(2)
            appdata.newy = yrange(2);
        end
        if appdata.newy < yrange(1);
            appdata.newy = yrange(1);
        end
        %     get the yfield handle
        set(appdata.yField, 'String', num2str(10^(-appdata.newy), '%0.4g'));
        set(appdata.yField, 'Userdata', 10^(-appdata.newy));
        set(appdata.hline, 'Xdata', appdata.initxrange, 'Ydata', repmat(appdata.newy, size(appdata.initxrange)));
        localSetAppData(hfig,appdata);
    end
elseif  strcmp(action,'up'), 
%     update the colors of the scatter plot
    appdata = updateCData(appdata, 1);
    appdata = updateLists(appdata);

    localSetAppData(hfig,appdata);
    set(hfig,'WindowButtonMotionFcn',{@movelines, 'motion', 0});
    set(hfig,'WindowButtonUpFcn','');
    set(hfig,'Pointer','arrow');
elseif strcmp(action,'editx')
    clearSelections([], [], hfig);
    newx = str2double(get(appdata.xField,'String'));
    if isempty(newx) || isnan(newx) || newx <= 1
        newx = get(appdata.xField,'Userdata');
        set(appdata.xField,'String',num2str(newx, '%0.5g'));
        return;
    end
    appdata.newx = log2(newx);

    if appdata.newx > appdata.xlim(2)
        appdata.newx =  appdata.xlim(2);
        set(appdata.xField,'String',num2str(2^appdata.newx, '%0.5g'));
    end
    if appdata.newx <  appdata.xlim(1)
        appdata.newx =  appdata.xlim(1);
        set(appdata.xField,'String',num2str(2^appdata.newx, '%0.5g'));
    end
    set(appdata.xField,'String',2^(appdata.newx));
    set(appdata.xField,'Userdata',2^(appdata.newx));
    set(appdata.vline1,'XData', repmat(appdata.newx, size(appdata.inityrange)),'YData', appdata.inityrange);
    set(appdata.vline2,'XData', repmat(-appdata.newx, size(appdata.inityrange)),'YData', appdata.inityrange);
    appdata = updateCData(appdata, 1);
    appdata = updateLists(appdata);

    localSetAppData(hfig,appdata);
elseif strcmp(action,'edity'),
    clearSelections([], [], hfig);
    newy = str2double(get(appdata.yField,'String'));
    if isempty(newy) || isnan(newy) || newy < 0
        newy = get(appdata.yField,'Userdata');
        set(appdata.yField,'String',num2str(newy, '%0.4g'));
        return;
    end
    appdata.newy = -log10(newy);
     if appdata.newy > appdata.ylim(2)
        appdata.newy =  appdata.ylim(2);
        set(appdata.yField,'String',num2str(10^(-appdata.newy), '%0.4g'));
     end
    if appdata.newy <  appdata.ylim(1)
        appdata.newy =  appdata.ylim(1);
        set(appdata.yField,'String',num2str(10^(-appdata.newy), '%0.4g'));
    end
    set(appdata.yField,'Userdata',10^(-appdata.newy));
    set(appdata.hline, 'Xdata', appdata.initxrange, 'Ydata', repmat(appdata.newy, size(appdata.initxrange)));
    appdata = updateCData(appdata, 1);
    appdata = updateLists(appdata);

    localSetAppData(hfig,appdata);
elseif strcmp(action, 'update')
    movelines([], [], 'edity');
    movelines([], [], 'editx');
elseif strcmp(action, 'reset')
    clearSelections([], [], hfig);
    appdata.newx = log2(appdata.xThd);
    appdata.newy = -log10(appdata.yThd);
    set(appdata.hAxis, 'XLim',appdata.xlim,'YLim',appdata.ylim);
        
    % get the xfield handle
    set(appdata.xField, 'String', num2str(appdata.xThd, '%0.5g'));
    set(appdata.xField, 'Userdata', appdata.xThd);
    set(appdata.yField, 'String', num2str(appdata.yThd, '%0.5g'));
    set(appdata.yField, 'Userdata', appdata.yThd);

    set(appdata.vline1,'XData', repmat(appdata.newx, size(appdata.inityrange)),'YData', appdata.inityrange);
    set(appdata.vline2,'XData', repmat(-appdata.newx, size(appdata.inityrange)),'YData', appdata.inityrange);
    set(appdata.hline, 'Xdata', appdata.initxrange, 'Ydata', repmat(appdata.newy, size(appdata.initxrange)));
    
    appdata = updateCData(appdata, 1);
    appdata = updateLists(appdata);
    localSetAppData(hfig,appdata);
end
end
%%%%%%%%%%*******************************
function markselect(varargin) 
% The list boxes callback
[hobj,hfig] = gcbo;

appdata = localGetAppData(hfig);
tag = get(hobj,'Tag');

if length(tag)==2&&tag(2)=='l'
    if tag(1) == 'u';
      listt = appdata.upllist;
      listu = appdata.upplist;
      indices = appdata.upindices;
    else % down case
      listt = appdata.dnllist;
      listu = appdata.dnplist;
      indices = appdata.dnindices;
    end
    rowlabels = get(listt, 'String');
elseif length(tag)==2 && tag(2)=='r'
   if tag(1) == 'u';
      listt = appdata.upplist;
      listu = appdata.upllist;
      indices = appdata.upindices;
    else % down case
      listt = appdata.dnplist;
      listu = appdata.dnllist;
      indices = appdata.dnindices;
   end
    rowlabels = get(listu, 'String');
end

V = get(listt,'Value'); 
set(listu,'Value',V);

xdata = appdata.effect;
ydata = appdata.significance;

t = zeros(length(V),1);
p = zeros(length(V),1);

for i = 1:length(V)   
    ind = V(i);
    id = indices(ind);

    t(i) = text('Parent',appdata.hAxis,...
        'Position',[xdata(id) ydata(id)],...
        'String',rowlabels{ind},...
        'HorizontalAlignment','left',...
        'VerticalAlignment', 'bottom',...
        'interpreter','none',...
        'clipping', 'on');

    setlabelposition(t(i), appdata.hAxis)
    set(t(i), 'ButtonDownFcn',{@localDeleteMatching, tag});

    p(i) = line('Parent',appdata.hAxis,...
        'XData',xdata(id),...
        'YData',ydata(id),...
        'Color',appdata.markerRGB,...
        'Marker','o',...
        'linewidth', 1.5);
end
 
if tag(1) == 'u'
    appdata.upSelectedText = deleteSelectedHandle(appdata.upSelectedText);
    appdata.upSelectedMarkers = deleteSelectedHandle(appdata.upSelectedMarkers);
    appdata.upSelectedText = t;
    appdata.upSelectedMarkers = p;
elseif tag(1) == 'd'
    appdata.dnSelectedText = deleteSelectedHandle(appdata.dnSelectedText);
    appdata.dnSelectedMarkers = deleteSelectedHandle(appdata.dnSelectedMarkers);
    appdata.dnSelectedText = t;
    appdata.dnSelectedMarkers = p;
end

localSetAppData(hfig,appdata);
end
%*********************************
function localDeleteMatching(varargin)
% Delete the label text box when click on the box
[t,hfig] = gcbo;
appdata = localGetAppData(hfig);
if nargin > 2
    tag = varargin{3};
    if length(tag)==2&&tag(1)=='u'
        [appdata.upSelectedText, appdata.upSelectedMarkers, appdata.upllist, appdata.upplist] =...
            deleteSelectedList(t, appdata.upSelectedText, appdata.upSelectedMarkers, appdata.upllist, appdata.upplist);
    elseif length(tag)==2&&tag(1)=='d'
        [appdata.dnSelectedText, appdata.dnSelectedMarkers, appdata.dnllist, appdata.dnplist] =...
            deleteSelectedList(t, appdata.dnSelectedText, appdata.dnSelectedMarkers, appdata.dnllist, appdata.dnplist);
    end
end

localSetAppData(hfig,appdata);
end
%********* Export results *******************
function exportDiffResults(varargin) 
[~,hfig] = gcbo; 
    appdata = localGetAppData(hfig);

    labels = {'Up Regulated', 'Down Regulated', 'All Differential'};
    varnames = {'upstruct', 'downstruct', 'diffstruct'};
    
    pvc = 10^(-appdata.newy);
    fct = 2^(abs(appdata.newx));
    idx = getDiffIndexes(appdata, 0);
    diffPV = appdata.pvalues(idx);
    if appdata.glnumericFlag
        diffL = conv2numric(appdata.labels(idx));
    else
        diffL = appdata.labels(idx);
    end
    
    diffFC = appdata.effect(idx);
    [diffPV, idx] = sort(diffPV);
    
    diffstruct.Name = 'Differentially Expressed';
    diffstruct.PVCutoff = pvc;
    diffstruct.FCThreshold = fct;
    diffstruct.GeneLabels = diffL(idx);
    diffstruct.PValues = diffPV;
    diffstruct.FoldChanges = 2 .^ (abs(diffFC(idx))); 
    
    idx = diffFC < 0;
    
    diffstruct.FoldChanges(idx) = - diffstruct.FoldChanges(idx);
    
	upstruct.Name = 'Up regulated';
    upstruct.PVCutoff = pvc;
    upstruct.FCThreshold = fct;
    if appdata.glnumericFlag
        upstruct.GeneLabels = conv2numric(appdata.uplabels);
    else
        upstruct.GeneLabels = appdata.uplabels;
    end
    upstruct.PValues = appdata.uppv;
    upstruct.FoldChanges = 2.^appdata.upfc;
 
    downstruct.Name = 'Down regulated';
    downstruct.PVCutoff = pvc;
    downstruct.FCThreshold = fct;
    if appdata.glnumericFlag
        downstruct.GeneLabels = conv2numric(appdata.dnlabels);
    else
        downstruct.GeneLabels = appdata.dnlabels;
    end
    downstruct.PValues = appdata.dnpv;
    downstruct.FoldChanges = -2.^(abs(appdata.dnfc));
    
    if appdata.DataMatrixFlag 
        diffstruct.PValues = bioma.data.DataMatrix(diffstruct.PValues,...
                                diffstruct.GeneLabels, {'p-values'});
        diffstruct.FoldChanges = bioma.data.DataMatrix(diffstruct.FoldChanges,...
                                diffstruct.GeneLabels, {'FoldChanges'});
                            
        upstruct.PValues = bioma.data.DataMatrix(upstruct.PValues,...
                                upstruct.GeneLabels, {'p-values'});
        upstruct.FoldChanges = bioma.data.DataMatrix(upstruct.FoldChanges,...
                                upstruct.GeneLabels, {'FoldChanges'});
        
        downstruct.PValues = bioma.data.DataMatrix(downstruct.PValues,...
                                downstruct.GeneLabels, {'p-values'});
        downstruct.FoldChanges = bioma.data.DataMatrix(downstruct.FoldChanges,...
                                downstruct.GeneLabels, {'FoldChanges'});
    end
  
    items = {upstruct, downstruct, diffstruct};
    
    export2wsdlg(labels, varnames, items, 'Export to Workspace');
end
%*************************************
function closeplot(h,varargin) %#ok<INUSD>
% Callback function to close the plot
close(gcbf)
end
%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%%
function localSetAppData(hfig,appdata)
setappdata(hfig,'MAVolvanoPlot',appdata);
end
%************************
function [appdata] = localGetAppData(hfig)

if isappdata(hfig,'MAVolvanoPlot')
    appdata = getappdata(hfig,'MAVolvanoPlot');
else
    appdata = guihandles(hfig);
    appdata.glnumericFlag = false;
    
    appdata.pvalues = [];
    appdata.labels = [];
    appdata.significance = [];
    appdata.effect = [];
    

    appdata.xlim = [];
    appdata.ylim = [];
    appdata.upMax = 1;
    appdata.dnMax = 1;
    appdata.uplabels = [];
    appdata.dnlabels = [];
    appdata.uppv = [];
    appdata.dnpv = [];
    appdata.upfc = [];
    appdata.dnfc = [];
  
    appdata.upSelectedText = [];
    appdata.upSelectedMarkers = [];
    appdata.dnSelectedText = [];
    appdata.dnSelectedMarkers = [];
    appdata.upllist = [];
    appdata.upplist = [];
    appdata.dnllist =[];
    appdata.dnplist = [];
    appdata.upindices = [];
    appdata.dnindices = [];
    
    appdata.sigRGB = [0.7 0.7 0.9];
    appdata.insRGB = [0.4 0.4 0.7];
    appdata.markerRGB = [1.0 0.3 0.0];
    appdata.lineRGB = [1.0 0.3 0.0];
    appdata.colormap = [];
    
    appdata.DataMatrixFlag = false;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do the plot
    appdata.hAxis = [];
    appdata.hPlot = [];
    appdata.vline1 = [];
    appdata.vline2 = [];
    appdata.hline = [];

    appdata.xThd = [];
    appdata.yThd = [];
    appdata.newx = [];
    appdata.newy = [];
    appdata.ltype = 0;
    appdata.tempfig = [];
    %%%%%%%UI%%%%%%%%%%
    appdata.xField = [];
    appdata.yField = [];
    appdata.initxrange = [];
    appdata.inityrange = [];
end
end
%*****************************************************
function r = isinaxes(point, axes)
%ISINAXES determines whether or not a point is in an axes.
%   ISINAXES returns true if POINT is in AXES and false
%   if it is not. POINT is a CurrentPoint. This utility routine
%   works only checks X and Y coordinates.
xr = get(axes,'Xlim');
yr = get(axes,'Ylim');
cx = point(1,1);
cy = point(1,2);
if cx >= xr(1) && cx <= xr(2) && cy >= yr(1) && cy <= yr(2)
    r = true;
else
    r = false;
end
end
%************************************************************
function appdata = updateCData(appdata, action)
% Updates the CData property of the plot
diffid = getDiffIndexes(appdata,0);

appdata.colorCol = ones(numel(appdata.effect), 1);
appdata.colorCol(diffid) = 2;
if action == 1 % set the plot cdata
    set(appdata.hPlot, 'CData', appdata.colorCol);
end
end
%************************************************************
function diffidx = getDiffIndexes(appdata, action)
% Get the indices of differentially express values
if action == 0
diffidx = ((appdata.effect > appdata.newx | appdata.effect < -appdata.newx)...
    & appdata.significance > appdata.newy);
elseif action == 1
    diffidx = ((appdata.effect > appdata.newx )...
    & appdata.significance > appdata.newy);
elseif action == 2
    diffidx = ((appdata.effect < -appdata.newx)...
    & appdata.significance > appdata.newy);
else
    diffidx = [];
end
end
%************************************************************
function appdata = updateLists(appdata)
% Update the listbox strings when the cutoff values and threshold changed
%Up list
upidx = getDiffIndexes(appdata, 1);
uppv = appdata.pvalues(upidx);
upfc = appdata.effect(upidx);
indices = find(upidx);

[uppv, idx] = sort(uppv);
appdata.uppv = uppv;
appdata.upfc = upfc(idx);
appdata.upindices = indices(idx);
set(appdata.upplist,'String',num2str(uppv(:), '%0.8g'), 'Max', appdata.upMax, 'Value',[]);

uplabels = appdata.labels(upidx);
appdata.uplabels = uplabels(idx);
set(appdata.upllist,'String',appdata.uplabels(:), 'Max', appdata.upMax, 'Value',[]);

%Down list
dnidx = getDiffIndexes(appdata, 2);
dnpv = appdata.pvalues(dnidx);
dnfc = appdata.effect(dnidx);
indices = find(dnidx);

[dnpv, idx] = sort(dnpv);
appdata.dnpv = dnpv;
appdata.dnfc = dnfc(idx);
appdata.dnindices = indices(idx);
set(appdata.dnplist,'String',num2str(dnpv(:), '%0.8g'), 'Max', appdata.dnMax, 'Value',[]);

dnlabels = appdata.labels(dnidx);
appdata.dnlabels = dnlabels(idx);
set(appdata.dnllist,'String',appdata.dnlabels(:), 'Max', appdata.dnMax, 'Value',[]);
end
%*******************************************************
function obj = deleteSelectedHandle(obj)
% Delete the selected handle objects
if ~isempty(obj)
    for n = 1:numel(obj)
        if ishandle(obj(n))
            delete(obj(n));
        end
    end
    obj = [];
end
end
%*****************************************************
function [selText, selMarkers, llist, rlist] = deleteSelectedList(t, selText, selMarkers, llist, rlist)
% Delete the list selection of the text object t, and from selText and
% selMarkers in appdata.
[r,c] = find(selText == t,1);%#ok

delete(selText(r));
delete(selMarkers(r));
selText(r) = [];
selMarkers(r) = [];

v = get(llist,'Value');
v(r) = [];
set(llist,'Value',v);
set(rlist,'Value',v);
end
%*************** Create plot **************
function appdata = createPlot(hFig, appdata)
% Create the plot. Type = 1 for plot to figure

set(hFig, 'colormap', appdata.colormap);
% Define axis for plot
appdata.hAxis = axes;

set(appdata.hAxis, 'Parent', hFig,...
    'XLimMode', 'manual',...
    'YLimMode', 'manual',...
    'XLim',appdata.xlim,'YLim',appdata.ylim,...
    'Box', 'on');

if ~appdata.plotonly
    axisPos = [0.065 0.2 0.57 0.7];
    set(appdata.hAxis, 'Position',axisPos,...
                        'XAxisLocation', 'top')
end


appdata = updateCData(appdata, 0);

% Due to a performance bottleneck with scatter plot. Use Patch instead. See
% G306305.
% hPlot = scatter(effect, significance, 9, colorCol, 'filled');
appdata.hPlot=patch('xdata',appdata.effect,'ydata',appdata.significance,...
    'cdata',appdata.colorCol,...
    'linestyle','none','markeredgecolor','flat',...
    'facecolor','none','marker','.','MarkerSize',9,...
    'DisplayName', 'Score');
set(get(get(appdata.hPlot,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
set(appdata.hPlot, 'parent', appdata.hAxis);

xlabel('log2 ( ratio )')
ylabel('-log10(p-values)')

appdata.initxrange = get(appdata.hAxis,'Xlim');
appdata.inityrange = get(appdata.hAxis,'Ylim');
% Create cutofflines
xvertical1  = repmat(appdata.newx,  size(appdata.inityrange));
xvertical2  = repmat(-appdata.newx,  size(appdata.inityrange));

yhorizontal = repmat(appdata.newy,size(appdata.initxrange));
appdata.vline1 = line(xvertical1,appdata.inityrange,...
                                'LineStyle','-.',...
                                'Linewidth', 1.2,...
                                'Color',appdata.lineRGB,...
                                'Tag', 'vertline',...
                                'DisplayName', 'Fold change threshold');
appdata.vline2 = line(xvertical2,appdata.inityrange,...
                        'LineStyle','-.',...
                        'Linewidth', 1.2,...
                        'Tag', 'vertline',...
                        'Color',appdata.lineRGB);
set(get(get(appdata.vline2,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
appdata.hline = line(appdata.initxrange,yhorizontal,...
                        'LineStyle','-.',...
                        'Linewidth', 1.2,...
                        'Color',appdata.lineRGB,...
                        'Tag', 'horzline',...
                        'DisplayName', 'p-value cutoff');

set(appdata.vline1,'ButtonDownFcn',{@movelines, 'down', 2});
set(appdata.vline2,'ButtonDownFcn',{@movelines, 'down', 3});
set(appdata.hline,'ButtonDownFcn',{@movelines, 'down', 1});

% If labels exist, set up a buttondown function
if ~appdata.plotonly
    set(appdata.hAxis,'ButtonDownFcn',@clickonplot);
end
end
%***********Update print**************
function resetFigureTools(fig)
% helper function to set UI menus
oldSH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

% delete figure menus not used
%h1 = findall(fig,'Type','uimenu', 'Label','&Edit');
h1 = findall(fig,'Type','uimenu', 'Tag','figMenuEdit');
%h2 = findall(fig,'Type','uimenu', 'Label','&View');
h2 = findall(fig,'Type','uimenu', 'Tag','figMenuView');
%h3 = findall(fig,'Type','uimenu', 'Label','&Insert');
h3 = findall(fig,'Type','uimenu', 'Tag','figMenuInsert');
%h4 = findall(fig,'Type','uimenu', 'Label','&Desktop');
h4 = findall(fig,'Type','uimenu', 'Tag','figMenuDesktop');
delete([h1,h2,h3,h4])

% Repair "File" menu
%hw = findall(fig,'Type','uimenu', 'Label','&File');
hw = findall(fig,'Type','uimenu', 'Tag','figMenuFile');
hf = get(hw,'children');
%h1 = findall(hw,'Label','Expo&rt Setup...');
h1 = findall(hw,'Tag','figMenuFileExportSetup');
%h2 = findall(hw,'Label','Print Pre&view...');
h2 = findall(hw,'Tag','figMenuFilePrintPreview');
%h3 = findall(hw,'Label','&Print...');
h3 = findall(hw,'Tag','printMenu');
delete(setxor(hf,[h1,h2,h3]))

set(h2, 'Callback', {@printfcns, fig, 'preview'});
set(h3, 'Callback', {@printfcns, fig, 'print'});

uimenu(hw,'Label','Export to Workspace...','Position',1, 'Callback', @exportDiffResults);
uimenu(hw ,'Label','Print to Figure','Position',2, 'Callback', {@printttofigure, fig});
uimenu(hw,'Label','Exit','Separator','on','Position',6,'Callback',@closeplot)
set(h1,'Separator','on')

%  Repair "Tool" Menu
%hw = findall(fig,'Type','uimenu','Label','&Tools');
hw = findall(fig,'Type','uimenu','Tag','figMenuTools');
hf = get(hw,'children');
h1 = findall(hw,'Tag','figMenuZoomIn');
h2 = findall(hw,'Tag','figMenuZoomOut');
h3 = findall(hw,'Tag','figMenuPan');
h4 = findall(hw,'Tag','figMenuResetView'); 
h5 = findall(hw,'Tag','figMenuOptions');
set([h1,h4],'separator','off')
delete(setxor(hf,[h1,h2,h3,h4,h5]))
delete(findall(h5,'Tag','figMenuOptionsDatatip'))
delete(findall(h5,'Tag','figMenuOptionsDataBar'))

% Repair "Help" menu
%hw = findall(fig,'Type','uimenu','Label','&Help');
hw = findall(fig,'Type','uimenu','Tag','figMenuHelp');
delete(get(hw,'children'));
uimenu(hw,'Label','Bioinformatics Toolbox Help','Position',1,'Callback',...
       'helpview(fullfile(docroot,''toolbox'',''bioinfo'',''bioinfo.map''),''bioinfo_product_page'')')
uimenu(hw,'Label','Volcano Plot Help','Position',2,'Callback',...
       ['helpview(fullfile(docroot,''toolbox'',''bioinfo'', ''bioinfo.map'')'...
        ',''mavolcanoplot_refpage'')' ])
uimenu(hw,'Label','Examples','Position',3,'Separator','on',...
       'Callback','demo(''toolbox'',''bioinfo'')')   
tlbx = ver('bioinfo');
mailstr = ['web(''mailto:bioinfo-feedback@mathworks.com?subject=',...
           'Feedback%20for%20MAVolcanoPlot%20in%20Bioinformatics',...
           '%20Toolbox%20',tlbx(1).Version,''')'];
uimenu(hw,'Label','Send Feedback','Position',4,'Separator','on',...
       'Callback',mailstr);
   
% Turn off toolbar
set(fig, 'Toolbar','none');
set(0,'ShowHiddenHandles',oldSH)
end

function printfcns(hSrc, event, hfig, action)  
%PRINTFCN Create and display a print preview of MSViewer spectra.

appdata = localGetAppData(hfig);

if ~isempty(appdata.tempfig) && ishandle(appdata.tempfig)
    delete(appdata.tempfig)
    appdata.tempfig = [];
end
appdata.tempfig = figure('Visible','off','colormap', appdata.colormap);
h_axes = copyobj(appdata.hAxis, appdata.tempfig);
set(h_axes,'XAxisLocation', 'bottom');

switch action
    case 'pagesetup'
        pagesetupdlg(appdata.tempfig);
    case 'preview'
        printpreview(appdata.tempfig);    
    case 'print'
        printdlg(appdata.tempfig); %#ok<MCPRD>
        close(appdata.tempfig);
end
localSetAppData(hfig, appdata)
end
%**********************************
function printttofigure(hSrc, event, hFigure)
%PRINTTOFIGURE Print volcanoplot to a figure.  

% xlabels warn about rgb color not supported in painter's mode, but it did not
% have waning id.
ws = warning;
warning('off'); %#ok

appdata = localGetAppData(hFigure);
old_axes = appdata.hAxis;
new_fig = figure('Visible', 'off');

% Get the new axes position in the center of the figure.
new_axes = axes('Parent',new_fig, 'Units',get(old_axes,'Units'));
newpos = get(new_axes,'Position');
delete(new_axes)

h_axes = copyobj(old_axes, new_fig);
xlim = get(old_axes, 'Xlim');
ylim = get(old_axes,'Ylim');
set(h_axes,'Position', newpos,...
    'HandleVisibility','on',...
    'XAxisLocation', 'bottom',...
    'XLim', xlim,'YLim', ylim,...
    'ButtonDownFcn', '');

warning(ws);
set(new_fig, 'Visible', 'on', 'colormap', appdata.colormap);
end

%*************** layout UI ****************
function appdata = createUIs(hFig, appdata)
units = get(hFig, 'Units');
fgColor = get(hFig, 'Color');

% Create cutoff panel
hCutoffP = uipanel('parent', hFig, ...
            'units', units, ...
            'BackgroundColor',fgColor,...
            'title', 'Cutoff Values',...
            'titlePosition', 'centertop',...
            'BorderType', 'beveledout',...
            'position', [0.065 0.04 0.57 0.12]);
        
xFieldPos = [0.58 0.23 0.12 0.45];
xtp = [0.405 0.17 0.17 0.45];

yFieldPos = [0.225 0.23 0.15 0.45];
ytp = [0.008 0.16 0.21 0.45];
ytp2 = [0.375 0.16 0.02 0.45];

%  Pushbutton Reset : reset the values to default
uBPos = [0.72 0.2 0.13 0.55];
rBPos = [0.86 0.2 0.13 0.55];

% Create the ui compoenets in order
yText1  = uicontrol('Parent',hCutoffP,...
    'Style','text','Units',units,...
    'Position',ytp,...
    'BackgroundColor',fgColor,...
    'HorizontalAlignment','center',...
    'String','-log10(p-value ',...
    'Fontsize', 10);%#ok

appdata.yField = uicontrol('Parent',hCutoffP,...
    'Style','edit','Units',units,...
    'Position',yFieldPos,...
    'BackgroundColor','white',...
    'HorizontalAlignment', 'left',...
    'String',num2str(appdata.yThd),...
    'UserData', appdata.yThd,...
    'Tag', 'pvField',...
    'CallBack',{@movelines, 'edity'});

yText2  = uicontrol('Parent',hCutoffP,...
    'Style','text','Units',units,...
    'Position',ytp2,...
    'BackgroundColor',fgColor,...
    'HorizontalAlignment','left',...
    'String',')');%#ok

xText  = uicontrol('Parent',hCutoffP,...
    'Style','text','Units',units,...
    'Position',xtp,...
    'BackgroundColor',fgColor,...
    'HorizontalAlignment','center',...
    'String','Fold change',...
    'Fontsize', 10);%#ok

appdata.xField = uicontrol('Parent',hCutoffP,...
    'Style','edit','BackgroundColor','white',...
    'Units',units,...
    'Position',xFieldPos,...
    'HorizontalAlignment', 'left',...
    'String',num2str(appdata.xThd),...
    'UserData', appdata.xThd,...
    'Tag', 'fcField',...
    'CallBack',{@movelines, 'editx'});

hUpdateBtn = uicontrol('Parent',hCutoffP, ...
    'Units',units, ...
    'Position',uBPos, ...
    'String','Update', ...
    'TooltipString', 'Update cutoff values', ...
    'Callback','', ...
    'Tag','PushUpdate',...
    'Callback', {@movelines, 'update'});%#ok

hResBtn = uicontrol('Parent',hCutoffP, ...
    'Units',units, ...
    'Position',rBPos, ...
    'String','Reset', ...
    'TooltipString', 'Reset to default cutoff values', ...
    'Callback','', ...
    'Tag','PushReset',...
    'Callback', {@movelines, 'reset'});%#ok

% Create label panel
hBaseP = uipanel('parent', hFig, ...
    'units', units, ...
    'BackgroundColor',fgColor,...
    'BorderType', 'beveledout',...
    'position', [0.65 0.04 0.34 0.92]);

% Pushbutton to export: export to matlab workspace
eBPos = [0.6 0.018 0.25 0.05];
cBPos = [0.15 0.018 0.25 0.05];

hClrBtn = uicontrol('Parent',hBaseP, ...
    'Units',units,...
    'Position',cBPos,...
    'String', 'Clear',...
    'TooltipString', 'Clear selected gene labels', ...
    'Tag','PushReset',...
    'Callback', {@clearSelections, hFig});%#ok

hExpBtn = uicontrol('Parent',hBaseP, ...
    'Units',units,...
    'Position',eBPos,...
    'String', 'Export...',...
    'TooltipString', 'Export to MATLAB workspace', ...
    'Tag','PushReset',...
    'Callback', @exportDiffResults);%#ok

% Create Listboxes, frames and tags for the frames
frametags = {'Down Regulated', 'Up Regulated'};

[appdata.upllist, appdata.upplist] = createframes(hBaseP, fgColor, units, frametags{2}, 2);
[appdata.dnllist, appdata.dnplist] = createframes(hBaseP, fgColor, units, frametags{1}, 1);
end
%******************************************************
function hlist = createlistbox(hFig, units, pos, tag)
 hlist = uicontrol('Parent',hFig, ...
         'Units',units, ...
         'BackgroundColor','w', ...
         'Position', pos, ...
         'HorizontalAlignment', 'left',...
         'String',[], ...
         'Style','listbox', ...
         'Tag',tag, ...
         'Max', 3,...
         'Visible','on',...
         'Callback', @markselect);
end    
 %******************************************************
function [hlist1, hlist2] = createframes(hParent, fgColor, units, frametag, fn)
 
 pP =[0.025 0.16 0.95 0.38];
 lP = [0.02 0.02, 0.48, 0.83];
 labelP = [lP(1)+lP(3)/4, 0.86, 0.3, 0.1];
 
 hpanel = uipanel('parent', hParent, ...
            'units', units, ...
            'position', [pP(1), (fn-1)*0.43 + pP(2), pP(3), pP(4)],...
            'title', frametag,...
            'titlePosition', 'centertop',...
            'BackgroundColor', fgColor,...
            'Tag',  ['frame' num2str(fn)]);

hllistname = uicontrol('Parent',hpanel, ... 
    'BackgroundColor',fgColor, ...
    'Units',units, ...
    'Position', labelP, ...
    'String', 'Genes', ...
    'HorizontalAlignment','left', ...
    'Style','text');%#ok

hrlistname = uicontrol('Parent',hpanel, ...
    'BackgroundColor',fgColor, ...
    'Units',units, ...
    'Position', [labelP(1)+0.5 labelP(2:4)], ...
    'String', 'p-values', ...
    'HorizontalAlignment','left',...
    'Style','text');%#ok


pos1 = lP;
pos2 = [lP(1)+0.48, lP(2:4)];
if fn == 1 % down regulated genes
    tag1 = 'dl';
    tag2 = 'dr';
elseif fn == 2 % up regulated genes
    tag1 = 'ul';
    tag2 = 'ur';
end

hlist1 = createlistbox(hpanel, units, pos1, tag1);
hlist2 = createlistbox(hpanel, units, pos2, tag2);
end
%****************************************
function numv = conv2numric(cellv)
    numv = char(cellv(:));
    numv = str2num(numv); %#ok
end   
%******************************************
function setlabelposition(htext, haxes)
% If the text extent outside the axis, move it inside

% Get text extent,a nd axes range
ext = get(htext, 'Extent');
xrange = get(haxes,'Xlim');
yrange = get(haxes,'Ylim');

if ext(1)+ext(3) > xrange(2)
    set(htext, 'HorizontalAlignment','right')
elseif ext(2) > yrange(2)
    set(htext, 'VerticalAlignment','top')
else
    return;
end
end
%-------------------------------------------------------
function appdata = setColorScheme(appdata, colorScheme, sigRGB, insRGB)
switch(colorScheme)
    case 0 % user input rgb 
        appdata.sigRGB = sigRGB;
        appdata.insRGB = insRGB;
    case 1 % blue
        appdata.sigRGB = [0.7 0.7 0.9];
        appdata.insRGB = [0.4 0.4 0.7];
    case 2 %cyan
        appdata.sigRGB = [0.8 1.0 1.0];
        appdata.insRGB = [0.6 0.6 0.9];
    case 3 %green
        appdata.sigRGB = [0.2 0.8 0.2];
        appdata.insRGB = [0.8 1.0 0.8];
    case 4 %orange
        appdata.sigRGB = [1.0 0.85 0.6];
        appdata.insRGB = [1.0 0.55 0.3];
    case 5 %pink
        appdata.sigRGB = [1.0 0.8 1.0];
        appdata.insRGB = [1.0 0.4 0.0];
    case 6 % orangegrey
        appdata.sigRGB = [0.7 0.7 0.9];
        appdata.insRGB = [1.0 0.4 0.6];
end
appdata.colormap = [appdata.sigRGB; appdata.insRGB];
end
%------------------------------------------
function inputStruct = parse_inputs(varargin)
% Parse input PV pairs.

% Check for the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:mavolcanoplot:IncorrectNumberOfArguments', mfilename))
end

% Allowed inputs
okargs = {'logtrans', 'pcutoff', 'foldchange', 'labels', 'colors','plotonly'};

% Defaults
inputStruct.logFlag = false;
inputStruct.pvCutoff = 0.05;    % p = 0.05
inputStruct.fcThd = 2;          % Two fold line
inputStruct.labels = [];
inputStruct.colorScheme = 1;
inputStruct.sigRGB = [];
inputStruct.insRGB = [];
inputStruct.plotOnly = false;

for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1
            inputStruct.logFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 2  % p-value cutoff
            if ~isnumeric(pval) || pval < 0
                error(message('bioinfo:mavolcanoplot:PValueCutoffMustBeNumeric'));
            end
            inputStruct.pvCutoff = pval;
        case 3  % foldchangethresh
            if ~isnumeric(pval) || pval < 0
                error(message('bioinfo:mavolcanoplot:FoldChangeThreshMustBeNumeric'));
            end
            inputStruct.fcThd = pval;
        case 4 % labels
            inputStruct.labels = pval;
        case 5 % colors
            if isnumeric(pval)
                if isequal(size(pval), [2 3]) & all(pval<=1)& all(pval >=0)%#ok
                    inputStruct.sigRGB = pval(1,:);
                    inputStruct.insRGB = pval(2,:);
                    inputStruct.colorScheme = 0;
                else
                    warning(message('bioinfo:mavolcanoplot:BadRGBColorValues'));
                end
            elseif ischar(pval)
                okcolors = {'blue', 'cyan', 'green', 'orange', 'pink','redgrey'};
                nc = strcmpi(pval, okcolors);
                if isempty(nc)
                    warning(message('bioinfo:mavolcanoplot:BadColors'));
                elseif length(nc) > 1
                    warning(message('bioinfo:mavolcanoplot:AmbiguousColors', pval));
                else
                    inputStruct.colorScheme = nc;
                end
            else
                warning(message('bioinfo:mavolcanoplot:BadColorInput'));
            end
        case 6 % plotonly
            inputStruct.plotOnly = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end
end
end

function [X, Y, PV, paramStruct] = check_inputdata(X, Y, PV, paramStruct)
% Check input data type is numerical and contain the same number of genes
% or rows

% Validate X, Y, and PV and GL
if ~isnumeric(X) || ~isreal(X) || ~isnumeric(Y) || ~isreal(Y)
   error(message('bioinfo:mavolcanoplot:ExpressionValuesNotNumericAndReal')) 
end

if ~isnumeric(PV) || ~isreal(PV) || ~isvector(PV)
   error(message('bioinfo:mavolcanoplot:PValueNotNumericAndReal')) 
end

if size(X,1) ~= size(Y,1)
   error(message('bioinfo:mavolcanoplot:NotEqualNumberOfExpressionValues'))
end

if size(X,1) ~= size(PV,1)
   error(message('bioinfo:mavolcanoplot:InconsistentNumberOfExpressionValues'))
end

if ~isempty(paramStruct.labels) &&...
        (~isvector(paramStruct.labels) || numel(paramStruct.labels) < size(X,1))
    error(message('bioinfo:mavolcanoplot:NotEqualNumberOfGeneLabels'))
end

if ~isempty(paramStruct.labels)
    if ~ischar(paramStruct.labels) && isnumeric(paramStruct.labels)
        paramStruct.glnumericFlag = true;
        paramStruct.labels = textscan(sprintf('%-1d ',paramStruct.labels(:)),'%s');
        paramStruct.labels = paramStruct.labels{:};
    else
        paramStruct.glnumericFlag = false;
        paramStruct.labels = paramStruct.labels(:);
    end
else
    labels=1:numel(PV);
    paramStruct.labels = textscan(sprintf('%-1d ',labels(:)),'%s');
    paramStruct.labels = paramStruct.labels{:};   
    paramStruct.glnumericFlag = false;
end


% Handle the matrix input. Use its mean values per row
if size(X, 2) > 1
    X = mean(X,2);
end

if size(Y, 2) > 1
    Y = mean(Y,2);
end

% convert inputs to column vectors
X = X(:);
Y = Y(:);
PV = PV(:);

% discard any zero elements
allZeros = ((X == 0) | (Y == 0) | (PV == 0));
allNegative = ((X < 0) | (Y < 0) | (PV < 0));
allNaN = (isnan(X) | isnan(Y));

if any(allZeros)
    warning(message('bioinfo:mavolcanoplot:ZeroValues'));
end

if any(allNegative)
    warning(message('bioinfo:mavolcanoplot:NegativeValues'));
end

if paramStruct.pvCutoff < 0 || paramStruct.pvCutoff > max(PV) || paramStruct.pvCutoff < min(PV)
    warning(message('bioinfo:mavolcanoplot:InvalidPVCutoff', paramStruct.pvCutoff));
    paramStruct.pvCutoff = 0.05;
end

paramStruct.goodVals = ~(allZeros | allNegative | allNaN);

if paramStruct.logFlag
    X = log2(X);
    Y = log2(Y);
end

end
