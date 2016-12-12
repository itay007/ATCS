function hout = maloglog(Xdata,Ydata,varargin)
%MALOGLOG creates a loglog plot of microarray data.
%
%   MALOGLOG(X,Y) creates a loglog scatter plot of X vs Y. X and Y can be
%   MATLAB numeric matrices or DataMatrix objects.
%
%   MALOGLOG(...,'FACTORLINES',N) adds lines showing a factor of N change.
%
%   MALOGLOG(...,'TITLE',TITLE) allows you to specify a title for the plot.
%
%   MALOGLOG(...,'LABELS',LABELS) allows you to specify a cell array of
%   labels for the data. If LABELS are defined, then clicking on a point on
%   the plot will show the LABEL corresponding to that point.
%
%   MALOGLOG(...,Handle Graphics name/value) allows you to pass optional
%   Handle Graphics property name/property value pairs to the function. 
%
%   H = MALOGLOG(...) returns the handle to the plot.
%
%   Examples:
%
%       maStruct = gprread('mouse_a1wt.gpr');
%       Red = magetfield(maStruct,'F635 Median');
%       Green = magetfield(maStruct,'F532 Median');
%       maloglog(Red,Green,'title','Red vs Green');
%       % Add factorlines and labels
%       figure
%       maloglog(Red,Green,'title','Red vs Green',...
%               'FactorLines',2,'LABELS',maStruct.Names);
%       % Now create a normalized plot 
%       figure
%       maloglog(manorm(Red),manorm(Green),'title',...
%               'Normalized Red vs Green','FactorLines',2,...
%               'LABELS',maStruct.Names);
%
%   See also LOGLOG, MABOXPLOT, MAGETFIELD, MAIMAGE, MAINVARSETNORM,
%   MAIRPLOT, MALOWESS, MANORM, MATTEST, MAVOLCANOPLOT, MOUSEDEMO. 

% Copyright 2003-2008 The MathWorks, Inc.


%== Check input
import bioinfoprivate.*;
bioinfochecknargin(nargin,2,mfilename);
hgargs = {};
equalLine = false;
twoXLine = false;
twoXScale = 2;
titleString = '';
labels = '';

if nargin > 2
    
    if rem(nargin,2)== 1
        error(message('bioinfo:maloglog:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'factorlines','title','labels'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            % here we assume that these are handle graphics options
            hgargs{end+1} = pname;  %#ok
            hgargs{end+1} = pval; %#ok
        elseif length(k)>1
              error(message('bioinfo:maloglog:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % factorlines
                    equalLine = true;
                    twoXLine = true;
                    if pval == 0
                        equalLine = false;
                        twoXLine = false;
                    elseif pval == 1
                        twoXLine = false;
                    end
                    twoXScale = pval;
                case 2 % title
                    titleString = pval;
                case 3
                    labels = pval;
            end
        end
    end
end 

%== Check for DataMatrix
if isa(Xdata, 'bioma.data.DataMatrix')
    if isempty(labels)
        labels = Xdata.RowNames;
    end
    Xdata = Xdata.(':')(':');
end

if isa(Ydata, 'bioma.data.DataMatrix')
    if isempty(labels)
        labels = Ydata.RowNames;
    end
    Ydata = Ydata.(':')(':');
end

% convery inputs to column vectors
Xdata = Xdata(:);
Ydata = Ydata(:);

% discard any zero elements

allZeros = ((Xdata == 0) | (Ydata == 0));
allNegative = ((Xdata < 0) | (Ydata < 0));

if any(allZeros)
    warning(message('bioinfo:maloglog:ZeroValues'));
end
warnFlag = false;
if any(allNegative)
    warning(message('bioinfo:maloglog:NegativeValues'));
    warnFlag = true;
end

% make nicer bounds than the default loglog plot
upperBound = 1.5*max(max(Xdata(:)),max(Ydata(:)));
goodVals = ~(allZeros|allNegative);
lowerBound = 10^(floor(log10(min(min(Xdata(goodVals)),min(Ydata(goodVals))))));

% plot the figure
if warnFlag
    warnState = warning('off','MATLAB:Axes:NegativeDataInLogAxis');
end
hPlot = loglog(Xdata,Ydata,'+',hgargs{:});

% resize the axes
hAxis = get(hPlot,'parent');
set(hAxis,'Xlim',[lowerBound,upperBound],'Ylim',[lowerBound,upperBound]);

% plot line y = x;
if equalLine
    line( [lowerBound,upperBound],[lowerBound,upperBound],'color','k','linestyle','-.','linewidth',3)
end

% plot 2x line
if twoXLine
    line( twoXScale*[lowerBound,upperBound],[lowerBound,upperBound],'color','k','linestyle','-.','linewidth',2)
    line( [lowerBound,upperBound],twoXScale*[lowerBound,upperBound],'color','k','linestyle','-.','linewidth',2)
end

% Add a title
if ~isempty(titleString)
    title(titleString);
end

% If labels exist, set up a buttondown function
if ~isempty(labels)
    set(hAxis,'Userdata',hPlot,'ButtonDownFcn',@clickonplot);
    set(hPlot,'Userdata',labels,'ButtonDownFcn',@clickonplot);
end

if warnFlag
    warning(warnState);
end
% if output is requested, send out handle
if nargout > 0
    hout = hPlot;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clickonplot(varargin) 
% callback function highlights selected element and displays label.
[hObject,hFig] = gcbo;
% clean up any old labels
%hFig = localCleanUpLabel;
udata = get(hObject,'Userdata');
if ishandle(udata) 
    hObject  = udata;
end

xvals = get(hObject,'XData');
yvals = get(hObject,'YData');
hAxis = get(hObject,'parent');
point = get(hAxis,'CurrentPoint');
udata = get(hObject,'Userdata');

% find the closest point
[v, index] = min((((xvals - point(1,1)).^2)./(point(1,1)^2)) +...
    ((yvals - point(1,2)).^2)./(point(1,2)^2)); 

%highlight the point and set the title
if ~isempty(index)
    hHighlight = line(xvals(index),yvals(index),...
        'color','red','marker','d','Tag','loglogHighlight'); %#ok
    % get the value and label of the current point
    cpAct = get(hAxis,'CurrentPoint');
    
    % create a new text object -- start with it invisible
    htext = text(cpAct(1,1),cpAct(1,2) ,udata(index),'visible','off','interpreter','none');
    
    % give it an off white background, black text and grey border
    set(htext,	'BackgroundColor' , [1 1 0.933333],'Color' , [0 0 0],'EdgeColor' ,...
        [0.8 0.8 0.8], 'Tag','LogLogDataTip');
    % show the text
    set(htext, 'Visible','on')
    set(hFig,'WindowButtonUpFcn',@localCleanUpLabel);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hFig = localCleanUpLabel(varargin)
%localCleanUpLabel callback function to remove label from image

% get the handles to the figure, image and axis
hFig = gcbf;

% delete the old label if it exists
oldLabel = findobj(hFig,'Tag','LogLogDataTip');
if ~isempty(oldLabel)
    delete(oldLabel);
end
% delete the old label if it exists
oldHighlight = findobj(hFig,'Tag','loglogHighlight');
if ~isempty(oldHighlight)
    delete(oldHighlight);
end

