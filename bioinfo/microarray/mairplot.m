function [outIntensity,outRatio,hout] = mairplot(Xdata,Ydata,varargin)
%MAIRPLOT creates log Intensity vs Ratio scatter plot of microarray
% data.
%
%   MAIRPLOT(X,Y) creates a log10 of product mean intensities vs log2 of
%   intensity ratio scatter plot of X vs Y. X and Y can also be DataMatrix
%   objects.
%
%   MAIRPLOT(..., 'TYPE', TYPE) sets the plot type. TYPE can be 'IR'
%   (default) or 'MA'. 'MA' plot creates the log2 mean intensity vs log2
%   ratio.
% 
%   MAIRPLOT(...,'LOGTRANS', true) will log transform X and Y from linear
%   scale (default). Otherwise X and Y are assumed to be log transformed
%   accordingly.
% 
%   MAIRPLOT(...,'FACTORLINES',N) adds lines showing a factor of N change.
%
%   MAIRPLOT(...,'TITLE',TITLE) allows you to specify a title for the plot.
%
%   MAIRPLOT(...,'LABELS',LABELS) allows you to specify a cell array of
%   labels for the data. If LABELS are defined, then clicking on a point on
%   the plot will show the LABEL corresponding to that point.
%
%   MAIRPLOT(...,'NORMALIZE',true) displays the lowess normalized ratio
%   values.
%
%   MAIRPLOT(...,'LOWESSOPTIONS',LOWESSOPTS) allows you to pass a cell
%   array of parameter/value pairs options to the lowess normalization. See
%   the help for MALOWESS for more details.
% 
%   MAIRPLOT(...,'SHOWPLOT',TF) does not display the plot when TF is FALSE.
%   The default is TRUE.
% 
%   MAIRPLOT(...,'PLOTONLY',TF) only displays the plot without the UI
%   components when TF is TRUE. The default is FALSE. If SHOWPLOT option is
%   FALSE, this option will be ignored.
%
%   [INTENSITY, RATIO] = MAIRPLOT(...) returns the intensity products and
%   ratio values. If the NORMALIZE option is used then the output RATIO
%   values are normalized values.  If X and Y are a DataMatrix objects,
%   INTENSITY and RATIO are DataMatrix objects with the same properties as
%   X or Y.
%
%   [INTENSITY, RATIO, H] = MAIRPLOT(...) returns the handle of the plot.
%
%   Examples:
%       maStruct = gprread('mouse_a1wt.gpr');
%       cy5data = magetfield(maStruct,'F635 Median');
%       cy3data = magetfield(maStruct,'F532 Median');
%       mairplot(cy5data,cy3data,'title','R vs G IR plot')
% 
%       % Add labels
%       mairplot(cy5data,cy3data,'title','R vs G  IR plot',...
%                'Labels',maStruct.Names)
% 
%       % Normalize the plot using lowess normalization
%       mairplot(cy5data,cy3data, 'type', 'MA',... 
%                'title','Normalized R vs G MA plot',...
%                'Normalize',true, 'Labels',maStruct.Names)
% 
%       % Return intensity products and ratios without showing the plot
%       [intensities, ratios] = mairplot(cy5data,cy3data,'Showplot', false);
% 
%       % View the normalized MA plot without UI components
%       [M, A] = mairplot(cy5data,cy3data, 'Normalize', true,...
%                                          'Plotonly', true,...
%                                          'Type', 'MA');  
% 
%   See also MABOXPLOT, MAGETFIELD, MAIMAGE, MAINVARSETNORM, MALOGLOG,
%   MALOWESS, MANORM, MATTEST, MAVOLCANOPLOT, MOUSEDEMO.

%   MAIRPLOT(...,'COLORS', COLORS) allows user to set the color scheme for
%   the ir plot. Options are 'blue', 'cyan', 'green', 'orange', 'pink', or
%   'redblue'(default). If COLORS is a matrix with exactly two rows, three
%   columns, and elements in the interval 0 to 1, The first row is the RGB
%   color for the significant genes and the second row is the RGB color
%   insignificant ones.  The columns of the input matrix represent
%   intensity of red, blue and green, respectively.

% Copyright 2003-2011 The MathWorks, Inc.


% References: 
% [1] Quackenbush J. "Microarray Data Normalization and Transformation" Nature
%     Genetics Suppl. 2002, Vol.32, pp496-501.
% [2] Dudoit S, Yang YH, Callow MJ, Speed TP. "Statistical Methods for
%     Identifying Differentially Expressed Genes in Replicated cDNA Microarray
%     Experiments" Statistica Sinica 2002, 12, pp 111-139.  

bioinfochecknargin(nargin,2,mfilename)

% Check input parameters
paramStruct = parse_inputs(varargin{:});

% No plot and no output
if nargout == 0 && ~paramStruct.showplotFlag
    return;
end
%== DataMatrix input
paramStruct.DataMatrixFlag = false;

if isa(Xdata, 'bioma.data.DataMatrix')
    if isempty(paramStruct.labels)
        paramStruct.labels = Xdata.RowNames;
    end
    Xdata = Xdata.(':')(':');
    paramStruct.DataMatrixFlag = true;
end

if isa(Ydata, 'bioma.data.DataMatrix')
    if isempty(paramStruct.labels)
        paramStruct.labels = Ydata.RowNames;
    end
    Ydata = Ydata.(':')(':');
    paramStruct.DataMatrixFlag = true;
end

% Check input data and validate parameters
[Xdata, Ydata, paramStruct] = check_inputdata(Xdata, Ydata, paramStruct);
%-------------------------Start here---------------------
units = 'Normalized';
hFig = figure('Units',units, 'Tag', 'mairplot',...
                             'Visible', 'off');
                         
%==Init appdata
appdata = localGetAppData(hFig);
appdata.hfig = hFig;
appdata.origcolormap = get(hFig, 'colormap');
appdata.plotonly = paramStruct.plotOnlyFlag;
appdata.labelnumericFlag = paramStruct.labelnumericFlag;
appdata.logflag = paramStruct.logFlag;
appdata = setColorScheme(appdata, paramStruct.colorScheme, paramStruct.rgbMap);
appdata.irflag = paramStruct.irFlag;
appdata.normflag = paramStruct.normFlag;
appdata.titlestr = paramStruct.titleString;

appdata.labels = paramStruct.labels(paramStruct.goodVals);
appdata.normoptions = getnormoptions(paramStruct.normOpts);
appdata.X = Xdata(paramStruct.goodVals);
appdata.Y = Ydata(paramStruct.goodVals);

%== Computation
appdata.logint = calculateXPlotData(appdata);

if paramStruct.normFlag
    appdata = donormalization(appdata);
else
    appdata.logratio = calculateYPlotData(appdata);
    appdata.smoothr=appdata.logratio;
end

%== Limits
appdata = updateAxesLim(appdata, 0);
appdata.upMax = numel(appdata.logint);
appdata.dnMax = numel(appdata.logint);

appdata.yThd = paramStruct.twoXScale;
appdata.newy = log2(appdata.yThd);
appdata.ltype = 0;

%== Do Plot and UI
if paramStruct.showplotFlag
    if ~paramStruct.plotOnlyFlag
        fw = 0.4; % The figure width in relation to the screen width
        fh = 0.5;
        if strncmpi(computer, 'GLNX', 4)
            fw = 0.7;
            fh = 0.6;
        end
        set(hFig, 'PaperUnits','normalized', ...
            'PaperPosition',[.55 .30 .5 .5], ...
            'Position',[.15 .25 fw fh],...
            'NumberTitle','on',...
            'Name', 'MAIRPlot',...
            'DeleteFcn', @deletePlot);

        % Reset figure menu items
        resetFigureTools(hFig);

        % create UI
        appdata = createUIs(hFig, appdata);
    else
        dcm_obj = datacursormode(hFig);
        set(dcm_obj,'UpdateFcn',{@locaDataCursorUpdate, hFig})
    end

    % create plot
    appdata = updateLists(appdata);

    appdata = createPlot(hFig, appdata);

    % Save state
    localSetAppData(hFig,appdata);

    set(hFig,'WindowButtonMotionFcn', {@movelines, 'motion', appdata.ltype },...
        'Visible','on');
    drawnow;
end

% if output is requested, send out handle
if nargout > 0
    hout = hFig;
    
    if paramStruct.DataMatrixFlag 
        if appdata.irflag
            colNames = {'I', 'R'};
        else
            colNames = {'M', 'A'};
        end
        outIntensity= bioma.data.DataMatrix(appdata.logint,...
                                appdata.labels, colNames(1));
        outRatio = bioma.data.DataMatrix(appdata.logratio,...
                                appdata.labels, colNames(2));
    else
        outIntensity = appdata.logint;
        outRatio = appdata.logratio;
    end
end
end
%--------------- Helper and callback functions -------------------%
function txt = locaDataCursorUpdate(hobj,eobj, hfig)  %#ok<*INUSL>
% Update the datacursor text.
tg = get(eobj, 'Target');
appdata = localGetAppData(hfig);

if appdata.irflag
   xl = 'I: ';
   yl = 'R: ';
else
   xl = 'A: ';
   yl = 'M: ';
end

pos = get(eobj,'Position');
dataidx = get(eobj, 'DataIndex');

if strcmpi(get(tg, 'Type'), 'patch')    
    txt = {['Label: ', appdata.labels{dataidx}],...
        [xl,num2str(pos(1))],...
        [yl,num2str(pos(2))]};
elseif strcmpi(get(tg, 'Type'), 'line')
    type = get(tg, 'Tag');
    if strcmpi(type, 'horzline')
         txt = {[yl, num2str(pos(2))]};
    else
        txt = [];
        return;
    end
else
    txt = [];
    return;
end
end
% ---------------------
function localSetAppData(hfig,appdata)
setappdata(hfig,'MAIRPlot',appdata);
end
%---------------------------------------------------
function [appdata] = localGetAppData(hfig)
if isappdata(hfig,'MAIRPlot')
    appdata = getappdata(hfig,'MAIRPlot');
else
    appdata = guihandles(hfig);
    appdata.labelnumericFlag = false;
    
    appdata.irflag = true;
    appdata.labels = [];    
    appdata.logratio = [];
    appdata.logint = [];
    appdata.smoothr = [];
    
    appdata.xlim = [];
    appdata.ylim = [];
    appdata.upMax = 1;
    appdata.dnMax = 1;
    appdata.uplabels = [];
    appdata.dnlabels = [];
  
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
    
    appdata.colormap = [];
    appdata.colorscheme = 1;
    appdata.markerRGB = [0.5 1.0 0.3];
    appdata.lineRGB = [1.0 0.3 0.0];
    
    % Do the plot
    appdata.hAxis = [];
    appdata.hPlot = [];
    appdata.hline1 = [];
    appdata.hline2 = [];

    appdata.yThd = [];
    appdata.newy = [];
    appdata.ltype = 0;
    appdata.tempfig = [];
    %%%%%%%UI%%%%%%%%%%
    appdata.yField = [];
    appdata.initxrange = [];
    appdata.plotonly = true;
end
end
%----------------Reset Figure menu and toolbar ---------------------%
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
uimenu(hw,'Label','Exit','Separator','on','Position',6,'Callback', @deletePlot)
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
uimenu(hw ,'Label','Set LOWESS Normalization Options...','Position',6,...
       'Separator', 'on',...
       'Callback', {@getlowessnormoptions, fig});

% Repair "Help" menu
%hw = findall(fig,'Type','uimenu','Label','&Help');
hw = findall(fig,'Type','uimenu','Tag','figMenuHelp');
delete(get(hw,'children'));
uimenu(hw,'Label','Bioinformatics Toolbox Help','Position',1,'Callback',...
       'helpview(fullfile(docroot,''toolbox'',''bioinfo'',''bioinfo.map''),''bioinfo_product_page'')')
uimenu(hw,'Label','IR Plot Help','Position',2,'Callback',...
       ['helpview(fullfile(docroot,''toolbox'',''bioinfo'', ''bioinfo.map'')'...
        ',''mairplot_refpage'')' ])
uimenu(hw,'Label','Examples','Position',3,'Separator','on',...
       'Callback','demo(''toolbox'',''bioinfo'')')   
tlbx = ver('bioinfo');
mailstr = ['web(''mailto:bioinfo-feedback@mathworks.com?subject=',...
           'Feedback%20for%20MAIRPlot%20in%20Bioinformatics',...
           '%20Toolbox%20',tlbx(1).Version,''')'];
uimenu(hw,'Label','Send Feedback','Position',4,'Separator','on',...
       'Callback',mailstr);
   
% Turn off toolbar
set(fig, 'Toolbar','none');
set(0,'ShowHiddenHandles',oldSH)
end
%-------------------------------------------------------------%
function appdata = setColorScheme(appdata, colorScheme, rgbMap)
switch(colorScheme)
    case 0 % user input rgb
        appdata.colormap = rgbMap;
        appdata.colorscheme = rgbMap;
    case 1 % blue
        appdata.colormap = [0.7 0.7 0.9; 0.4 0.4 0.7];
        appdata.colorscheme = 'blue';
    case 2 %cyan
        appdata.colormap = [0.8 1.0 1.0; 0.6 0.6 0.9];
        appdata.colorscheme = 'cyan';
    case 3 %green
        appdata.colormap = [0.8 1.0 0.8; 0.2 0.8 0.2];
        appdata.colorscheme = 'green';
    case 4 %orange
        appdata.colormap = [1.0 0.85 0.6; 1.0 0.55 0.3];
        appdata.colorscheme = 'orange';
    case 5 %pink
        appdata.colormap = [1.0 0.8 1.0; 1.0 0.4 0.0];
        appdata.colorscheme = 'pink';
    case 6 % redblue
        appdata.colormap = [0.55 0.65 1; 1.0 0.6 0.4];
        appdata.colorscheme = 'redblue';
end
end
%************************************************************
    function appdata = updateCData(appdata, action)
        % Updates the CData property of the plot
        appdata.cdata = ones(numel(appdata.logint), 1);
        if appdata.plotonly 
            chkBoxVal = 1;
        else
            chkBoxVal = get(appdata.checkbox, 'value');
        end

        if chkBoxVal == 0
            appdata.cdata(1) = 2;
        else
            diffid = getDiffIndexes(appdata,0);
            appdata.cdata(diffid) = 2;
        end

        if ~isequal(action, 0) % set the plot cdata
            set(appdata.hPlot, 'CData', appdata.cdata);
        end
    end
%------------------------------------------------------------%
function appdata = updateAxesLim(appdata, action)
    if action == 0
        appdata.xlim = [min(appdata.logint)-0.2, max(appdata.logint)+0.2];
        appdata.ylim = [min(appdata.logratio)-0.2, max(appdata.logratio+0.2)];
    elseif action == 1
        appdata.xlim = [min(appdata.logint)-0.2, max(appdata.logint)+0.2];
    elseif  action == 2
        appdata.ylim = [min(appdata.logratio)-0.2, max(appdata.logratio+0.2)];
    end
end
%---------------------- Create plot -----------------------%
function appdata = createPlot(hFig, appdata)
% Create the plot. Type = 1 for plot to figure
set(hFig, 'colormap', appdata.colormap);

% Define axis for plot
appdata.hAxis = axes('Parent', hFig,...
    'XLimMode', 'manual',...
    'YLimMode', 'manual',...
    'XLim',appdata.xlim,'YLim',appdata.ylim,...
    'Box', 'on');

if ~appdata.plotonly
    axisPos = [0.065 0.27 0.65 0.65];
    set(appdata.hAxis, 'Position',axisPos);% ,...
end 

appdata = updateCData(appdata, 0);

% Due to a performance bottleneck with scatter plot. Use Patch instead.
appdata.hPlot=patch('xdata',appdata.logint,'ydata',appdata.logratio,...
                    'cdata',appdata.cdata,...
                    'linestyle','none','markeredgecolor','flat',...
                    'facecolor','none','marker','.',...
                    'DisplayName', 'Ratio');
set(get(get(appdata.hPlot,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
set(appdata.hPlot, 'parent', appdata.hAxis);

if appdata.irflag
    xlabel('log10(Intensity)');
    ylabel('log2(Ratio)');
else
    xlabel('A');
    ylabel('M');
end 

title(appdata.titlestr);

Xrange = get(appdata.hAxis,'Xlim');

% Create cutofflines
appdata.initxrange = [Xrange(1),Xrange(2)];
yhorizontal1 = repmat(appdata.newy, size(appdata.initxrange));
yhorizontal2 = repmat(-appdata.newy, size(appdata.initxrange));
yequal = zeros(size(appdata.initxrange));
line(appdata.initxrange,yequal,'LineStyle','-.',...
                               'Linewidth', 1.2,...
                               'Color','k',...
                               'Tag', 'horzline',...
                               'DisplayName', 'Zero ratio');

appdata.hline1 = line(appdata.initxrange,yhorizontal1,'LineStyle','-.',...
                                                      'Linewidth', 1.2,...
                                                      'Color',appdata.lineRGB,...
                                                      'Tag', 'horzline',...
                                                      'DisplayName', 'Ratio thredshold');
appdata.hline2 = line(appdata.initxrange,yhorizontal2,'LineStyle','-.',...
                                                      'Linewidth', 1.2,...
                                                      'Color',appdata.lineRGB,...
                                                      'Tag', 'horzline');
set(get(get(appdata.hline2,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
appdata.smoothline = line(appdata.logint,appdata.smoothr,'LineStyle','none',...
                          'marker', 'x',...
                          'Color','r',...
                          'visible', 'off',...
                          'Tag', 'horzline');
set(get(get(appdata.smoothline,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
set(appdata.hline1,'ButtonDownFcn',{@movelines, 'down', 1});
set(appdata.hline2,'ButtonDownFcn',{@movelines, 'down', 2});

% If labels exist, set up a buttondown function
if ~appdata.plotonly
    set(appdata.hAxis, 'ButtonDownFcn',@clickonplot);
    set(appdata.hPlot, 'ButtonDownFcn',@clickonplot);
end
end
%--------------------------------------------------------------------%
    function updateFactorLines(hSrc, hEvt) %#ok<*INUSD>
        hfig = gcbf;
        appdata = localGetAppData(hfig);
        chkVal = get(appdata.checkbox, 'Value');
        if chkVal == 1
            set(appdata.hline1, 'visible', 'on');
            set(appdata.hline2, 'visible', 'on');          
        elseif chkVal ==0
            set(appdata.hline1, 'visible', 'off');
            set(appdata.hline2, 'visible', 'off');
        end
        appdata = updateCData(appdata, 1);
        drawnow;
        localSetAppData(hfig, appdata);
    end
   %--------------------------------------------------------------------%
    function updateSmoothLine(hSrc, hEvt)
        hfig = gcbf;
        appdata = localGetAppData(hfig);
        chkVal = get(appdata.smoothcbox, 'Value');
        if chkVal == 1
            [xdata, idx] = sort(appdata.logint);
            ydata = appdata.smoothr(idx);
            
            if appdata.normflag 
                ydata = zeros(length(appdata.logint), 1);
            end
            
            set(appdata.smoothline, 'xdata', xdata', 'ydata', ydata);
            set(appdata.smoothline, 'visible', 'on');       
        elseif chkVal ==0
            set(appdata.smoothline, 'visible', 'off');
        end
        localSetAppData(hfig, appdata);
    end
%--------------------------------------------------------------------%
function appdata = updateLists(appdata)
% Update the listbox strings when the threshold changed
%Up list
upidx = getDiffIndexes(appdata, 1);
appdata.upindices = find(upidx);
appdata.uplabels = appdata.labels(upidx);
set(appdata.upllist,'String',appdata.uplabels(:), 'Max', appdata.upMax, 'Value',[]);

%Down list
dnidx = getDiffIndexes(appdata, 2);
appdata.dnindices = find(dnidx);
appdata.dnlabels = appdata.labels(dnidx);
set(appdata.dnllist,'String',appdata.dnlabels(:), 'Max', appdata.dnMax, 'Value',[]);
end
%---------------------------------------------------------%
function diffidx = getDiffIndexes(appdata, action)
% Get the indices of differentially express values
if action == 0
    diffidx = (appdata.logratio > appdata.newy | appdata.logratio < -appdata.newy);
elseif action == 1
    diffidx = appdata.logratio > appdata.newy;
elseif action == 2
    diffidx = appdata.logratio < -appdata.newy;
else
    diffidx = [];
end
end
%----------------------------------------------------------%
function markselect(varargin) 
% The list boxes callback
[hobj,hfig] = gcbo;
appdata = localGetAppData(hfig);
tag = get(hobj,'Tag');

if tag(1) == 'u';
    listt = appdata.upllist;
    indices = appdata.upindices;
else % down case
    listt = appdata.dnllist;
    indices = appdata.dnindices;
end
rowlabels = get(listt, 'String');

V = get(listt,'Value'); 

xdata = appdata.logint;
ydata = appdata.logratio;

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
        'Interpreter','none',...
        'clipping', 'on','visible', 'on');

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
%----------------------------------------------------------------%
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
%------------------------------------------------------------%
    function changecursor(hsrc, hevt)
        hfig = gcbf;
        set(hfig, 'pointer', 'watch')
        return;
       
    end
%------------------------------------------------------------%
function movelines(varargin)
[hobj,hfig] = gcbo; 

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
    clearSelections([], [], appdata);
    
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
        fuzzy = 0.01 * (yrange(2) - yrange(1));
        online = cy > yrange(1) && cy < yrange(2) && ...
            cx > xrange(1) && cx < xrange(2) &&...
            ((cy > appdata.newy - fuzzy && cy < appdata.newy + fuzzy) ||...
             (cy > -appdata.newy - fuzzy && cy < -appdata.newy + fuzzy));
        
        if online && strcmp(cursorstate,'arrow'),
            set(hfig,'Pointer','crosshair');
        elseif ~online && strcmp(cursorstate,'crosshair'),
            set(hfig,'Pointer','arrow');
        end
 
    elseif appdata.ltype == 1 || appdata.ltype == 2
        set(hfig, 'Pointer', 'crosshair');
        % Check the newpoints are within range
        appdata.newy = cp(1,2);
        if appdata.newy > yrange(2)
            appdata.newy = yrange(2);
        end
        if appdata.newy < yrange(1);
            appdata.newy = yrange(1);
        end

        if appdata.ltype == 1 && appdata.newy <= 0
            appdata.newy = 1;
        end
    
        if appdata.ltype == 2 && appdata.newy >= 0
            appdata.newy = 1;
        end
        
        appdata.newy = abs(appdata.newy);
        
        %     get the xfield handle
        set(appdata.yField, 'String', num2str(2^(appdata.newy), '%0.5g'));
        set(appdata.yField, 'Userdata', 2^(appdata.newy));
        
        set(appdata.hline1,'XData',appdata.initxrange, 'YData',repmat(appdata.newy, size(appdata.initxrange)));
        set(appdata.hline2,'XData',appdata.initxrange, 'YData',repmat(-appdata.newy, size(appdata.initxrange)));
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
elseif strcmp(action,'edity'),
    clearSelections([], [], appdata);
    newy = str2double(get(appdata.yField,'String'));
    if isempty(newy) || isnan(newy) || newy <=1
        newy = get(appdata.yField,'Userdata');
        set(appdata.yField,'String',num2str(newy, '%0.4g'));
        return;
    end
    appdata.newy = log2(newy);
    if appdata.newy > appdata.ylim(2)
        appdata.newy =  appdata.ylim(2);
        set(appdata.yField,'String',num2str(2^appdata.newy, '%0.4g'));
    end
    if appdata.newy <  appdata.ylim(1)
        appdata.newy =  appdata.ylim(1);
        set(appdata.yField,'String',num2str(2^appdata.newy, '%0.4g'));
    end
    set(appdata.yField,'Userdata',2^(appdata.newy));
    set(appdata.hline1,'XData',appdata.initxrange, 'YData',repmat(appdata.newy, size(appdata.initxrange)));
    set(appdata.hline2,'XData',appdata.initxrange, 'YData',repmat(-appdata.newy, size(appdata.initxrange)));

    appdata = updateCData(appdata, 1);
    appdata = updateLists(appdata);
    localSetAppData(hfig,appdata);
    drawnow;
elseif strcmp(action, 'update')
    movelines([], [], 'edity');
elseif strcmp(action, 'reset')
    appdata = clearSelections([], [], appdata);

    appdata.newy = log2(appdata.yThd);
    set(appdata.hAxis, 'XLim',appdata.xlim,'YLim',appdata.ylim);
        
    % get the yfield handle
    set(appdata.yField, 'String', num2str(appdata.yThd, '%0.5g'));
    set(appdata.yField, 'Userdata', appdata.yThd);

    set(appdata.hline1,'XData',appdata.initxrange, 'YData',repmat(appdata.newy, size(appdata.initxrange)));
    set(appdata.hline2,'XData',appdata.initxrange, 'YData',repmat(-appdata.newy, size(appdata.initxrange)));
    
    appdata = updateCData(appdata, 1);
    appdata = updateLists(appdata);
    localSetAppData(hfig,appdata);
    drawnow;
end
end
%---------------------------------------------------------------%
function appdata = clearSelections(varargin) 
    if nargin > 2
        appdata = varargin{3};
    else
        hfig = gcbf;
        appdata = localGetAppData(hfig);
    end
    
appdata.upSelectedText = deleteSelectedHandle(appdata.upSelectedText);
appdata.upSelectedMarkers = deleteSelectedHandle(appdata.upSelectedMarkers);
appdata.dnSelectedText = deleteSelectedHandle(appdata.dnSelectedText);
appdata.dnSelectedMarkers = deleteSelectedHandle(appdata.dnSelectedMarkers);
set(appdata.upllist, 'Value', []);
set(appdata.dnllist, 'Value', []);
end
%-----------------Input normalization options--------------------%
function getlowessnormoptions(hSrc, event, hfig) 
   hdialog = findobj('Tag', 'Options for LOWESS');
   if ~isempty(hdialog)
       set(hdialog, 'visible', 'on');
       return;
   end
   
   appdata = localGetAppData(hfig);
 
   options = {num2str(appdata.normoptions{2}), num2str(appdata.normoptions{4}),...
              appdata.normoptions{6}};
   
   options = lowessoptiondlg(options);
   
   if isempty(options)
       return;
   end
   
   if options{1} ~= 1 && options{1} ~= 2
       options{1} = 1;
   end
   options{2} = str2double(options{2});
   options{3} = bioinfoprivate.opttf(options{3});
      
   err_str = '';
   
   if isempty(options{1}) || isnan(options{1}) || ~isreal(options{1})
       err_str = sprintf('Order must be numeric and real.\n');
   end 
   
   if isempty(options{3})
       err_str = strcat(err_str, sprintf('Use robust must be a logical value.\n'));
   end 
   
   if isempty(options{2}) || isnan(options{2}) || ~isreal(options{2})
       err_str = strcat(err_str, sprintf('Span must be numeric and real'));
   end 
   
   if isempty(err_str)
       appdata.normoptions{2} = options{1};
       appdata.normoptions{4} = options{2};
       appdata.normoptions{6} = options{3};
       appdata.normflag = false;
       set(appdata.hnormbtn, 'enable', 'on'); 
       localSetAppData(hfig, appdata)
       return
   end
   
   % wait for error dialog and go back to the input dialog
   uiwait(errordlg(err_str, name));
   getlowessnormoptions([], [], hfig)
end
%---------------------------------------------------------------%
function options = getnormoptions(optionpv)
    options = {'order', 1, 'span', 0.05, 'robust', false};
    
    if isempty(optionpv)
        return;
    end

    okpara = {'order',   'span', 'robust'};

    for i = 1:2:length(optionpv)-1
        pname = optionpv{i};
        pval = optionpv{i+1};
        k = strmatch(lower(pname), okpara);
        if ~isempty(k) && length(k)==1
            switch(k)
                case 1  % order
                    options{2} = pval;
                case 2 % span
                    options{4} = pval;
                case 3 % robust
                    options{6} = pval;
                
            end
        end
    end
end
    %------------------------------------------------------------------------%
     function normalize(hSrc, event, hfig)     
         ButtonName = questdlg('Show normalized plot in another window?', ...
             'Normalization', 'Yes'); 
         pause(.2);
         drawnow;
         
         switch ButtonName,
             case 'Yes',
                 cursorstate = get(hfig,'Pointer');
                 set(hfig, 'Pointer', 'watch')
                 set(hfig,'WindowButtonMotionFcn', {@changecursor});            
                 drawnow;
                 
                 appdata = localGetAppData(hfig);
                 type = 'ir';
                 
                 if ~appdata.irflag
                     type = 'ma';
                 end
                 
                 [int, ratio] = mairplot(appdata.X, appdata.Y, 'type', type,...
                     'logtrans', appdata.logflag,...
                     'factorline', appdata.yThd,...
                     'normalize', true,...
                     'LOWESSOPTIONS', appdata.normoptions,...
                     'title', appdata.titlestr,...
                     'labels', appdata.labels,...
                     'Colors', appdata.colorscheme);
                  
                  set(appdata.smoothcbox, 'Enable', 'on')
                  appdata.smoothr = appdata.logratio - ratio;
                  localSetAppData(hfig, appdata);
                  set(hfig, 'Pointer', cursorstate)
                  set(hfig,'WindowButtonMotionFcn', []);
                  drawnow;
                 return;
             case 'No',
                 cursorstate = get(hfig,'Pointer');
                 set(hfig, 'Pointer', 'watch')
                 set(hfig,'WindowButtonMotionFcn', {@changecursor});                             
                 drawnow;
                 
                 appdata = localGetAppData(hfig);
                 appdata = donormalization(appdata);
                 set(appdata.hPlot, 'YData', appdata.logratio);

                 appdata = updateAxesLim(appdata, 2);
                 appdata = updateLists(appdata);
                 appdata = updateCData(appdata, 1);
  
                 set(appdata.hAxis, 'XLim',appdata.xlim,'YLim',appdata.ylim)
                 
                 appdata.normflag = true;
                 set(appdata.hnormbtn, 'Enable', 'off');
                 set(appdata.smoothcbox, 'Enable', 'off')
                 
                 localSetAppData(hfig, appdata)
                 set(hfig, 'Pointer', cursorstate)
                 set(hfig,'WindowButtonMotionFcn', []);
                 drawnow;
                 return;                
             case 'Cancel',
                 return;
         end % switch
     end          
 %-------------------------------------------------------------------
         function appdata = donormalization(appdata)
             appdata = clearSelections([], [], appdata);
             ratio = calculateYPlotData(appdata);
             appdata.smoothr = malowess(appdata.logint, ratio, appdata.normoptions{:});
             appdata.logratio = ratio - appdata.smoothr;
         end
%--------------- Export results ---------------------------------%
function exportDiffResults(varargin) 
[obj,hfig] = gcbo; 
    appdata = localGetAppData(hfig);
   
    labels = {'Up Regulated', 'Down Regulated', 'All Differential'};
    varnames = {'upstruct', 'downstruct', 'diffstruct'};

    idx = getDiffIndexes(appdata, 0);
    diffstruct.Name = 'Differentially Expressed';
    diffstruct.FCThreshold = 2^(abs(appdata.newy));
    diffstruct.Indices = find(idx);
    if appdata.labelnumericFlag
        diffstruct.Labels = conv2numeric(appdata.labels(idx));
    else
        diffstruct.Labels = appdata.labels(idx);
    end
   

	upstruct.Name = 'Up regulated';
    upstruct.FCThreshold = 2^(abs(appdata.newy));
    upstruct.Indices = appdata.upindices;
    if appdata.labelnumericFlag
        upstruct.Labels = conv2numeric(appdata.uplabels);
    else
        upstruct.Labels = appdata.uplabels;
    end
 
    downstruct.Name = 'Down regulated';
    downstruct.FCThreshold = 2^(abs(appdata.newy));
    downstruct.Indices = appdata.dnindices;
    
    if appdata.labelnumericFlag
        downstruct.Labels = conv2numeric(appdata.dnlabels);
    else
        downstruct.Labels = appdata.dnlabels;
    end

    items = {upstruct, downstruct, diffstruct};
    
    export2wsdlg(labels, varnames, items, 'Export to Workspace');
end
%----------------------------------------------------------------%
function printttofigure(hSrc, event, hFigure)
%PRINTTOFIGURE Print irplot to a figure.  

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
    'XLim', xlim,'YLim', ylim,...
    'ButtonDownFcn', []);
% == Remove buttondownfcn from the patch
hplot=findobj(get(h_axes, 'Children'), 'Type', 'patch');
if ~isempty(hplot)
    set(hplot, 'ButtonDownFcn', []);
end
warning(ws);
set(new_fig, 'colormap', appdata.colormap, 'Visible', 'on');
end
%------------------------------------------------------------%
function printfcns(hSrc, event, hfig, action)
%PRINTFCN Create and display a print preview

appdata = localGetAppData(hfig);

if ~isempty(appdata.tempfig) && ishandle(appdata.tempfig)
    delete(appdata.tempfig)
    appdata.tempfig = [];
end
appdata.tempfig = figure('Visible','off', 'colormap', appdata.colormap);
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
%------------- layout UI -------------------------------------------%
function appdata = createUIs(hFig, appdata)
units = get(hFig, 'Units');
fgColor = get(hFig, 'Color');

% Create factor line panel
hThrdP = uipanel('parent', hFig, ...
            'units', units, ...
            'BackgroundColor',fgColor,...
            'title', 'Threshold',...
            'titlePosition', 'centertop',...
            'BorderType', 'beveledout',...
            'position', [0.065 0.04 0.65 0.13]);

yFieldPos = [0.225 0.15 0.15 0.45];
ytp = [0.008 0.10 0.21 0.47];

%  Pushbutton Reset : reset the values to default
uBPos = [0.5 0.15 0.20 0.47];
rBPos = [0.74 0.15 0.20 0.47];

% Create the ui compoenets in order
yText  = uicontrol('Parent',hThrdP,...
    'Style','text','Units',units,...
    'Position',ytp,...
    'BackgroundColor',fgColor,...
    'HorizontalAlignment','center',...
    'String','Fold change ',...
    'Fontsize', 10);%#ok

appdata.yField = uicontrol('Parent',hThrdP,...
    'Style','edit','Units',units,...
    'Position',yFieldPos,...
    'BackgroundColor','white',...
    'HorizontalAlignment', 'left',...
    'String',num2str(appdata.yThd),...
    'UserData', appdata.yThd,...
    'Tag', 'pvField',...
    'CallBack',{@movelines, 'edity'});

cBoxPos = [0.01, 0.75, 0.35, 0.3];
appdata.checkbox = uicontrol('Parent',hThrdP,...
    'Style','Checkbox','Units',units,...
    'Position',cBoxPos,...
    'BackgroundColor',fgColor,...
    'HorizontalAlignment', 'left',...
    'String', 'Show factor lines',...
    'Tag', 'checkBox',...
    'value', 1,...
    'CallBack', {@updateFactorLines});

hUpdateBtn = uicontrol('Parent',hThrdP, ...
    'Units',units, ...
    'Position',uBPos, ...
    'String','Update', ...
    'TooltipString', 'Update cutoff values', ...
    'Tag','PushUpdate',...
    'Callback', {@movelines, 'update'});%#ok

hResBtn = uicontrol('Parent',hThrdP, ...
    'Units',units, ...
    'Position',rBPos, ...
    'String','Reset', ...
    'TooltipString', 'Reset to default cutoff values', ...
    'Tag','PushReset',...
    'Callback', {@movelines, 'reset'});%#ok

% Create label panel
hBaseP = uipanel('parent', hFig, ...
    'units', units, ...
    'BackgroundColor',fgColor,...
    'BorderType', 'beveledout',...
    'position', [0.73 0.04 0.26 0.92]);

% Pushbutton to export: export to matlab workspace
eBPos = [0.58 0.02 0.35 0.05];
cBPos = [0.13 0.02 0.35 0.05];
nBPos = [0.25, 0.925, 0.45, 0.05];
sCBPos = [0.18, 0.87, 0.75, 0.05];

appdata.hnormbtn = uicontrol('Parent',hBaseP, ...
    'Units',units,...
    'Position',nBPos,...
    'String', 'Normalize',...
    'TooltipString', 'LOWESS Normalization', ...
    'Tag','PushReset',...
    'Callback', {@normalize, hFig});

if appdata.normflag
    set(appdata.hnormbtn, 'Enable', 'off');
end

appdata.smoothcbox = uicontrol('Parent',hBaseP,...
    'Style','Checkbox','Units',units,...
    'Position',sCBPos,...
    'BackgroundColor',fgColor,...
    'HorizontalAlignment', 'left',...
    'String', 'Show smooth curve',...
    'Tag', 'checkBox',...
    'Enable', 'off',...
    'value', 0,...
    'CallBack', {@updateSmoothLine});

hClrBtn = uicontrol('Parent',hBaseP, ...
    'Units',units,...
    'Position',cBPos,...
    'String', 'Clear',...
    'TooltipString', 'Clear selected gene labels', ...
    'Tag','PushReset',...
    'Callback', @clearSelections);%#ok

hExpBtn = uicontrol('Parent',hBaseP, ...
    'Units',units,...
    'Position',eBPos,...
    'String', 'Export...',...
    'TooltipString', 'Export to MATLAB workspace', ...
    'Tag','PushReset',...
    'Callback', @exportDiffResults);%#ok

% Create Listboxes, frames and tags for the frames
frametags = {'Down Regulated Genes', 'Up Regulated Genes'};

appdata.upllist = createframes(hBaseP, fgColor, units, frametags{2}, 2);
appdata.dnllist = createframes(hBaseP, fgColor, units, frametags{1}, 1);
end
 %******************************************************
function hlist = createframes(hParent, fgColor, units, frametag, fn)
 
 pP =[0.025 0.10 0.95 0.36];
 lP = [0.02 0.02, 0.95, 0.9];
 
 hpanel = uipanel('parent', hParent, ...
            'units', units, ...
            'position', [pP(1), (fn-1)*0.38 + pP(2), pP(3), pP(4)],...
            'title', frametag,...
            'titlePosition', 'centertop',...
            'BackgroundColor', fgColor,...
            'Tag',  ['frame' num2str(fn)]);

if fn == 1 % down regulated genes
    tag = 'dl';
elseif fn == 2 % up regulated genes
    tag = 'ul';
end

hlist = uicontrol('Parent',hpanel, ...
         'Units',units, ...
         'BackgroundColor','w', ...
         'Position', lP, ...
         'HorizontalAlignment', 'left',...
         'String',[], ...
         'Style','listbox', ...
         'Tag', tag, ...
         'Max', 3,...
         'Visible','on',...
         'Callback', @markselect);
end
%--------------------------------------------------------------%
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
%--------------------------------------------------------------%
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
%----------------------------------------------------%
function clickonplot(varargin)
% callback function highlights selected element and displays label.
[hObject,hfig] = gcbo;
% clean up any old labels
appdata = localGetAppData(hfig);
appdata = clearSelections([], [], appdata);

if( hObject == appdata.hAxis)
    hObject = appdata.hPlot;
end

xvals = get(hObject,'XData');
yvals = get(hObject,'YData');
point = get(appdata.hAxis,'CurrentPoint');

% find the closest point
[v, index] = min((((xvals - point(1,1)).^2)./(point(1,1)^2)) +...
    ((yvals - point(1,2)).^2)./(point(1,2)^2));

%highlight the point and set the title
if ~isempty(index)
    hHighlight = line(xvals(index),yvals(index),...
                        'color',appdata.markerRGB,...
                        'marker','o',...
                        'linewidth', 1.5,...
                        'Tag','logHighlight');%#ok
                    
    htext = text(point(1,1),point(1,2) , appdata.labels(index));

    % give it an off white background, black text and grey border
    set(htext,	'HorizontalAlignment','left',...
                'VerticalAlignment', 'bottom',... 
                'interpreter','none',...
                'Tag','logDataTip',...
                'Visible','on');
            
    %Check extends
    setlabelposition(htext, appdata.hAxis)
    
    set(hfig,'WindowButtonUpFcn',@clearClickLabel);
    localSetAppData(hfig,appdata);
end

end
%---------------------------------------------------------%
function hFig = clearClickLabel(varargin)
%clearClickLabel callback function to remove label from image

% get the handles to the figure, image and axis
hFig = gcbf;

% delete the old label if it exists
oldLabel = findobj(hFig,'Tag','logDataTip');
if ~isempty(oldLabel)
    delete(oldLabel);
end
% delete the old label if it exists
oldHighlight = findobj(hFig,'Tag','logHighlight');
if ~isempty(oldHighlight)
    delete(oldHighlight);
end
end
%----------------------------------------------------------%
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
%---------------------------------------------------------%
function logratio = calculateYPlotData(appdata)
% Calculate M, A pair for MA plot. 
% Y - the log2 ratio of X1 and X2

ratio = (appdata.X ./appdata.Y);

if appdata.logflag
    logratio = log2(ratio);
else
    logratio = ratio;
end
end
%------------------------------------------------------------%

function logint = calculateXPlotData(appdata)
% Calculate M, A pair for MA plot. 
% Y - the log2 ratio of X1 and X2
% X - the mean log2 (X1 * X2) for MA, log10 (X1*X2) for IR 

intensity = (appdata.X.*appdata.Y);

if appdata.logflag
    if appdata.irflag
        logint = log10(intensity);
    else
        logint = log2(intensity)/2;
    end
else
    logint = intensity;
end
end
%----------------------------------------------%
    function deletePlot(obj,event)
        hfig = gcbf;
        if ~isempty(hfig) && ishandle(hfig)
            localSetAppData(hfig, []);
            delete(hfig);
        end
    end

%****************************************
function numv = conv2numeric(cellv)
    numv = char(cellv(:));
    numv = str2num(numv); %#ok
end

function inputStruct = parse_inputs(varargin)
% Parse input PV pairs.

% Check for the right number of inputs
if rem(nargin,2)== 1
   error(message('bioinfo:mairplot:IncorrectNumberOfArguments', mfilename))
end

% Allowed inputs
okargs = {'type', 'logtrans', 'factorlines','title','labels',...
          'normalize', 'lowessoptions', 'colors', 'showplot', 'plotonly'};

% Defaults
inputStruct.hgargs = {};
inputStruct.twoXScale = 2;
inputStruct.titleString = '';
inputStruct.normFlag = false;
inputStruct.normOpts = {};
inputStruct.labels = '';
inputStruct.showplotFlag = true;
inputStruct.plotOnlyFlag = false;
inputStruct.irFlag = true;
inputStruct.logFlag = true;
inputStruct.colorScheme = 6; % redblue
inputStruct.rgbMap = [];

for j=1:2:nargin
    pname = varargin{j};
    pval = varargin{j+1};
    k = find(strncmpi(pname, okargs, numel(pname)));

    if isempty(k)
        % here we assume that these are handle graphics options
        inputStruct.hgargs{end+1} = pname; 
        inputStruct.hgargs{end+1} = pval; 
    elseif length(k)>1
        error(message('bioinfo:mairplot:AmbiguousParameterName', pname));
    else
        switch(k)
            case 1 % type
                if ischar(pval)
                    okmethods = {'ir','ma'};
                    nm = strmatch(lower(pval), okmethods);
                    if isempty(nm)
                        error(message('bioinfo:mairplot:TypeNotValid'));
                    elseif length(nm)>1
                        error(message('bioinfo:mairplot:AmbiguousType', pval));
                    else
                        inputStruct.irFlag = nm == 1;
                    end
                else
                    error(message('bioinfo:mairplot:TypeNotValid'));
                end
            case 2 % log transform
                inputStruct.logFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            case 3  % factorlines
                if ~isnumeric(pval) || pval < 0
                    error(message('bioinfo:mairplot:FactorMustBePositiveNumber'));
                end
                inputStruct.twoXScale = pval;
            case 4 % title
                inputStruct.titleString = pval;
            case 5 % labels
                inputStruct.labels = pval;
            case 6 % normalization
                inputStruct.normFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            case 7 % lowessoptions
                inputStruct.normOpts = pval;
                if ~iscell(inputStruct.normOpts)
                   error(message('bioinfo:mairplot:NormOptsNotCell'));
                end
            case 8 % colors
                if isnumeric(pval)
                    if isequal(size(pval), [2 3]) & all(pval<=1)& all(pval >=0)%#ok
                        inputStruct.rgbMap = pval;
                        inputStruct.colorScheme = 0;
                    else
                        warning(message('bioinfo:mairplot:BadRGBColorValues'));
                    end
                elseif ischar(pval)
                    okcolors = {'blue', 'cyan', 'green', 'orange', 'pink','redblue'};
                    nc = strmatch(lower(pval), okcolors);
                    if isempty(nc)
                        warning(message('bioinfo:mairplot:BadColors'));
                    elseif length(nc) > 1
                        warning(message('bioinfo:mairplot:AmbiguousColors', pval));
                    else
                        inputStruct.colorScheme = nc;
                    end
                else
                    warning(message('bioinfo:mairplot:BadColorInput'));
                end
            case 9 % showplot
                inputStruct.showplotFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            case 10 % plotonly
                inputStruct.plotOnlyFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        end
    end
end
end

function [Xdata, Ydata, paramStruct] = check_inputdata(Xdata, Ydata, paramStruct)
% Check input data type is numerical and contain the same number of genes
% or rows

% Validate Xdata, Ydata
if ~isnumeric(Xdata) || ~isreal(Xdata) ||~isvector(Xdata)||...
        ~isnumeric(Ydata) || ~isreal(Ydata) || ~isvector(Ydata)
   error(message('bioinfo:mairplot:ExpressionValuesNotNumericAndRealVector'))
end

if size(Xdata,1) ~= size(Ydata,1)
   error(message('bioinfo:mairplot:NotEqualNumberOfExpressionValues'))
end

% Handle labels input
if ~isempty(paramStruct.labels) &&... 
        (~isvector(paramStruct.labels) || numel(paramStruct.labels) < size(Xdata,1))
    error(message('bioinfo:mairplot:NotEqualNumberOfLabels'))
end

if ~isempty(paramStruct.labels)
    if ~ischar(paramStruct.labels) && isnumeric(paramStruct.labels)
        paramStruct.labelnumericFlag = true;
        paramStruct.labels = strread(sprintf('%-1d ',paramStruct.labels(:)),'%s');
    else
        paramStruct.labelnumericFlag = false;
        paramStruct.labels = paramStruct.labels(:);
    end
else
    labels=1:numel(Xdata);
    paramStruct.labels = strread(sprintf('%-1d ',labels(:)),'%s');
    paramStruct.labelnumericFlag = false;
end

% Handle the matrix input. Use its mean values per row
if size(Xdata, 2) > 1
    Xdata = mean(Xdata,2);
end

if size(Ydata, 2) > 1
    Ydata = mean(Ydata,2);
end

% convert inputs to column vectors
Xdata = Xdata(:);
Ydata = Ydata(:);

% discard any zero elements
allZeros = ((Xdata == 0) | (Ydata == 0));
allNegative = ((Xdata < 0) | (Ydata < 0));

if any(allZeros)
    warning(message('bioinfo:mairplot:ZeroValues'));
end

if any(allNegative)
    warning(message('bioinfo:mairplot:NegativeValues'));
end

paramStruct.goodVals = ~(allZeros|allNegative);

if isempty(paramStruct.titleString)
    if paramStruct.irFlag
        paramStruct.titleString = 'IR Plot';
    else
        paramStruct.titleString = 'MA plot';
    end
end
end
