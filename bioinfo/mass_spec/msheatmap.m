function hFig = msheatmap(MZ,Y,varargin)
%MSHEATMAP creates a heat map image of a set of spectra
%
%   MSHEATMAP(MZ,Y) shows a heat map, or pseudocolor image, of the
%   spectrograms in Y. Y is a matrix with the ion intensities of several
%   spectra arranged columnwise, all sharing the same MZ scale. MZ is the
%   common mass-charge column vector with same number of elements as
%   rows in Y.
%
%   MSHEATMAP(MZ,T,Y) sets a vector with the retention time of a LCMS (or
%   GCMS) dataset. T is a numeric vector with the same number of elements
%   as rows in Y. 
%
%   MSHEATMAP(...,'RESOLUTION',RES) sets the horizontal resolution of the
%   image. MSHEATMAP analyzes the input MZ vector, when the spectra
%   contains large mass-charge ratios (m/z>2500) then RES defaults to 0.5,
%   otherwise RES defaults to 0.05. 
%
%   MSHEATMAP(...,'MIDPOINT',Q) sets the midpoint of the colormap at a
%   specified quantile of the ion intensities values. MSHEATMAP uses a 
%   custom colormap with cold colors for regions without peaks and hot
%   colors for peaks. The midpoint is white. When a retention time vector
%   (T) is provided (such as in LCMS and GCMS datasets), Q defaults to 0.99
%   w, i.e. only one percent of pixels are expected to represent molecular
%   peaks. Q defaults to 0.95 otherwise. Use the figure's interactive
%   colormap shift tool to further adjust the colormap in the pseudocolor
%   image.
%
%   MSHEATMAP(...,'MARKERS',M) sets a list of mass-charge markers.
%   Positions are marked along the top axis. Default is M = []. 
%
%   MSHEATMAP(...,'RANGE',R) uses a 1-by-2 vector (R) with the mass-charge
%   range for the desired heat map.
%
%   MSHEATMAP(...,'GROUP',G) sets the class label for every spectrogram
%   used to group the rows of the heat map. G can be a numeric vector or a
%   cell array of strings with the same number of elements as spectrograms
%   in Y. By default the spectrograms are not grouped.
%
%   MSHEATMAP(...,'SPECIDX',I) sets an arbitrary vector for the spectrogram
%   index. I can be a numeric vector or a cell array of strings with the
%   sample ID of every spectrogram. I must have the same number of elements
%   as spectrograms in Y. Default is 1:numberOfSpectrograms. 
%  
%   Examples:
%
%       % Load a SELDI-TOF dataset example.
%       load sample_lo_res
%       P = [3991.4 4598 7964 9160];
%         
%       % Show a heat map of a set of spectrograms limiting the
%       % viewing range and marking the location of the peaks.
%       msheatmap(MZ_lo_res,Y_lo_res,'markers',P,'range',[3000 10000])
%
%       % Show a heat map of a set of spectrograms grouping them
%       % into two sets.
%       groups = [1 1 2 2 1 1 2 2];
%       msheatmap(MZ_lo_res,Y_lo_res,'markers',P,'group',groups)
%
%       % Loading a LCMS dataset example.
%       load lcmsdata
%       [MZ,Y] = msppresample(ms_peaks,5000);
%       msheatmap(MZ,ret_time,log(Y))
%  
%   See also DIFFPROTDEMO, LCMSDEMO, MSALIGN, MSBACKADJ, MSLOWESS, MSNORM,
%   MSPREPRODEMO, MSRESAMPLE, MSSGOLAY, MSVIEWER. 

%   Copyright 2003-2008 The MathWorks, Inc.


% validate MZ, Y, and T (when necessary)

bioinfochecknargin(nargin,2,mfilename);

if ~isnumeric(MZ) || ~isreal(MZ) || ~isvector(MZ)
   error(message('bioinfo:msheatmap:MZNotNumericAndReal')) 
end

if numel(varargin) && isnumeric(varargin{1})
    T = Y;
    Y = varargin{1};
    varargin(1)=[];
    isLCMS = true;
else
    isLCMS = false;
end
    
if isLCMS && (~isreal(T) || ~isvector(T))
   error(message('bioinfo:msheatmap:TNotNumericAndReal')) 
end

if ~isnumeric(Y) || ~isreal(Y)
   error(message('bioinfo:msheatmap:IntensityNotNumericAndReal')) 
end

if numel(MZ) ~= size(Y,1)
   error(message('bioinfo:msheatmap:NotEqualNumberOfSamples'))
end

if isLCMS && (numel(T) ~= size(Y,2))
   error(message('bioinfo:msheatmap:NotEqualNumberOfTimePoints'))
end

[numSamples, numSpectrograms] = size(Y);  

% set defaults
colormapSize = 196;
B = [];
L = [0 inf];
classIdsProvided = false;
specidx = 1:numSpectrograms;
reMapYTickLabels = true;

if isLCMS 
    Q = 0.99;
else
    Q = 0.95;
end
if max(MZ)>2500
    res = 0.5;
else
    res = 0.05;
end
    
% get optional input arguments

if numel(varargin)
    if rem(numel(varargin),2)
        error(message('bioinfo:msheatmap:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'markers','group','range','resolution','specidx','midpoint','limits'};
    for j=1:2:numel(varargin)
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:msheatmap:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:msheatmap:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % markers
                    B = pval;
                    if ~isnumeric(B) || ~isreal(B) || ~isvector(B)
                        error(message('bioinfo:msheatmap:MarkersNotNumericAndReal'))
                    end
                    B=B(:);
                    if max(B)>max(MZ) || min(B)<min(MZ)
                        error(message('bioinfo:msheatmap:MarkersOutOfRange'))
                    end
                case 2 % class id
                    if isLCMS
                        warning(message('bioinfo:msheatmap:GroupingIgnored'))
                    else
                        ID = pval(:);
                        classIdsProvided = true;
                        if iscell(ID)
                            if ~iscellstr(ID)
                                error(message('bioinfo:msheatmap:InvalidCellForID'))
                            end
                            if ismember('',ID)
                                error(message('bioinfo:msheatmap:InvalidCellWithEmptyID'))
                            end
                            [ID,validGroups] = grp2idx(ID);
                        elseif isnumeric(ID) || islogical(ID)
                            if islogical(ID)
                                ID = double(ID);
                            end
                            if ~isreal(ID) || any(isnan(ID)) || any(rem(ID,1))
                                error(message('bioinfo:msheatmap:InvalidNumericForID'))
                            end
                            [ID,validGroups] = grp2idx(ID);
                        else
                            error(message('bioinfo:msheatmap:InvalidTypeForID'))
                        end
                        %check ID has the correct size
                        if numel(ID)~=numSpectrograms
                            error(message('bioinfo:msheatmap:InvalidSizeForID'))
                        end
                    end
                case 3 % limits (now is range)
                    if numel(pval)~=2 || diff(pval)<=0
                        error(message('bioinfo:msheatmap:badRange'))
                    end
                    L = pval;
                    L(1) = max(L(1),0);
                case 4 % resolution
                    if ~isnumeric(pval)|| ~isscalar(pval) || pval<=0
                        error(message('bioinfo:msheatmap:InvalidResolution'))
                    end
                    res = pval;
                case 5 % spec idx
                    if isLCMS
                        warning(message('bioinfo:msheatmap:IndexingIgnored'))
                    else
                        pval = pval(:);
                        if iscell(pval)
                            if ~iscellstr(pval)
                                error(message('bioinfo:msheatmap:InvalidCellForSPECIDX'))
                            end
                            if ismember('',pval)
                                error(message('bioinfo:msheatmap:InvalidCellWithEmptySPECIDX'))
                            end
                        elseif isnumeric(pval)
                            if ~isreal(pval) || any(isnan(pval))
                                error(message('bioinfo:msheatmap:InvalidNumericForSPECIDX'))
                            end
                        else
                            error(message('bioinfo:msheatmap:InvalidTypeForSPECIDX'))
                        end

                        if numel(pval)~=numSpectrograms
                            error(message('bioinfo:msheatmap:InvalidIndices'))
                        end
                    end
                    specidx = pval;
                case 6 % colormap midpoint
                    if ~isnumeric(pval)|| ~isscalar(pval) || pval<0 || pval>1
                        error(message('bioinfo:msheatmap:InvalidMidPoint'))
                    end
                    Q = pval;
                case 7 % limits (now is range, this case is for grandfathering)
                    if numel(pval)~=2 || diff(pval)<=0
                        error(message('bioinfo:msheatmap:badRange'))
                    end
                    L = pval;
                    L(1) = max(L(1),0);
            end
        end
    end
end

% new mz vector (for the image axis), make sure there is no negative points
newMZ = (round(max(L(1),max(0,min(MZ)))):res:round(min(L(2),max(MZ))))';

% remove any negative value
Y = max(0,Y);

% interpolate Y to the new mz values
if (numel(newMZ) == numel(MZ)+1) && max(MZ-newMZ(1:numel(MZ)))< sqrt(eps('single'))
    % Do nothing here. This is a special case due to behavior of
    % mspresample. We get N points from mspresample and N+1 points from
    % msheatmap when trying to resample to the same points.
else
    Y = single(interp1q(MZ,Y,newMZ));
end

% re-scaling to the size of the colormap
Y = ((Y - min(Y(:))) * (colormapSize-1)) ./ max(Y(:));

% calculate colormapmidpoint before reducing the precision of Y 
cmmp = quickQuantile(Y(:),Q)./(colormapSize-1);

if cmmp(1)<0.05
    warning(message('bioinfo:msheatmap:BadColormap'))
end

% reducing precision to the size of the colormap (saving memory space)
Y = uint8(Y);

% reorder Y if class ids were provided
if classIdsProvided 
    [ID,hID]=sort(ID);
    Y = Y(:,hID);
end

% create image
hFig = figure;
if isLCMS && max(abs(diff(diff(T))))<0.01
    hIm = image(double(newMZ), double(T), Y');
    reMapYTickLabels = false;
else
    hIm = image(double(newMZ), double(1:numSpectrograms), Y');
end
hAxes = get(hIm,'parent');
set(hAxes,'Tag','MSHeatMap');
hold on

xlabel('Mass/Charge (M/Z)')
if classIdsProvided 
    ylabel('Spectrogram Groups')
    q = find([diff(ID);1]);
    set(hAxes,'Ytick',q-[q(1);diff(q)]/2+0.5)
    set(hAxes,'YtickLabel',validGroups)
    for i = 1:length(q)-1
        hlin = plot([0 newMZ(end)],q([i i])+0.5,'r','Linewidth',2);
        setappdata(hlin,'legend_hgbehavior',0)
        set(hggetbehavior(hlin,'DataCursor'),'Enable',0)
    end
    reMapYTickLabels = false;
    setappdata(hFig,'msheatmapType','categorical')
    setappdata(hFig,'msheatmapCategoricalOriginalIDs',hID)
elseif isLCMS
    ylabel('Retention Time')
    setappdata(hFig,'msheatmapType','lcms')
    setappdata(hFig,'msheatmapRetentionTimeVector',T)
    setappdata(hFig,'msheatmapMZRange',[min(newMZ),max(newMZ)])
else
    ylabel('Spectrogram Indices')
    setappdata(hFig,'msheatmapType','indexed')
end
xlabel('Mass/Charge (M/Z)')

% draw markers
for i = 1:numel(B)
  hmar = plot(B(i),min(ylim)-diff(ylim)/100,'Vr','clipping','off',...
       'markersize',5,'linewidth',2,'Tag','msheatmapMarkers');
  setappdata(hmar,'legend_hgbehavior',0)
  set(hggetbehavior(hmar,'DataCursor'),'Enable',0)
end

if isLCMS
    specidx = num2cell(T);
elseif isnumeric(specidx)
    specidx = num2cell(specidx);
end
setappdata(hFig,'MyYTickLabels',specidx)
setupmsheatmapListeners(reMapYTickLabels,hFig,hAxes); % listeners to update markers when zooming and panning

hold off 

colormap(privateColorMap(cmmp,colormapSize)); % set colormap
hc = colorbar('Ytick',[]);
ylabel(hc,'Relative Intensity')   

% We intercept the callbacks because the colorbar is specially formatted
% for this figure
set(findall(hFig,'type','uimenu','tag','figMenuInsertColorbar'),'Callback',@myColorbarCallback)
set(findall(hFig,'Tag','Annotation.InsertColorbar'),'ClickedCallback',@myColorbarCallback)

setAllowAxesRotate(rotate3d(hFig),hAxes,false)
msDataCursor(hFig)

if nargout==0
    clear hFig
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function myColorbarCallback(~,~)
     insertmenufcn(gcbf,'Colorbar')
     hc = findall(gcbf,'type','axes','tag','Colorbar');
     if ~isempty(hc)
         set(hc,'Ytick',[]);
         ylabel(hc,'Relative Intensity')
     end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function limitsListener(hSrc,event,hf,ha) %#ok
% Update the markers on the top edge
xlim = get(ha,'Xlim');
ylim = get(ha,'Ylim');
hMarkers = findobj(ha,'Tag','msheatmapMarkers');
for i = 1:numel(hMarkers)
    x = get(hMarkers(i),'Xdata');
    if x>=min(xlim) && x<=max(xlim)
        set(hMarkers(i),'Visible','on')
        set(hMarkers(i),'Ydata',min(ylim)-diff(ylim)/100)
    else
        set(hMarkers(i),'Visible','off')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ytickListener(hSrc,event,hf,ha) %#ok
% Update YTickLabels
specidx = getappdata(hf,'MyYTickLabels');
ticks = get(ha,'YTick');
[h,g] = ismember(ticks,1:numel(specidx));
newLab = cell(numel(ticks),1);
newLab(h) = specidx(g(h));
set(ha,'YTickLabel',newLab);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setupmsheatmapListeners(reMapYTickLabels,hFig,hAxes)  
% helper function to setsup listeners for the markers, so we can detect if
% we would need to change the relocate

% listens when the Ylim of axes has changed
YLimListener = addlistener(hAxes,'YLim',...
               'PostSet',@(hSrc,event)limitsListener(hSrc,event,hFig,hAxes));
% listens when the Xlim of axes has changed
XLimListener = addlistener(hAxes,'XLim',...
               'PostSet',@(hSrc,event)limitsListener(hSrc,event,hFig,hAxes)); 
if reMapYTickLabels           
  % listens when the Ytick of axes has changed
  YTickListener = addlistener(hAxes,'YTick',...
               'PostSet',@(hSrc,event)ytickListener(hSrc,event,hFig,hAxes));           
else
    YTickListener = [];
end
% store the listeners
setappdata(hFig,'msheatmapListeners',[YLimListener, XLimListener, YTickListener]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pcmap = privateColorMap(mp,colormapSize)
%PRIVATECOLORMAP returns a custom color map colors with the midpoint
%(white) at the specified location

mp = double(mp); % can't set a colormap with singles 

% set proportions for every region
prop = [.15 .25 .05 .05 .05 .05 .25 .15];
% adjust prop according to the midpoint
prop = 2*prop .* [mp([1 1 1 1]) 1-mp([1 1 1 1])];
% set control points for colors
pts = [0  0 .5;
       0  0  1;
       0 .8 .8;
      .4  1 .4;
      .9  1 .9;
       1  1  0;
       1 .5  0;
       1  0  0;
      .4  0  0];
% calculate colormap
pcmap = interp1(cumsum([0 prop])*(colormapSize-1)+1,pts,1:colormapSize,'pchip');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function th = quickQuantile(P,quan_int,res_int)

% computes a quick quantile using accumarray. This avoids extra copies of
% the data that occur in quantile.m
if nargin <3
    res_int = max(round(numel(P)/1000),1000);
end

% We could use histc here but this seems to be quicker.
imi = min(P);
ima = max(P);
in2idx =  @(x,r) round((x - imi) / (ima-imi) * (r-1) + .5);
inva = accumarray(in2idx(P(~isnan(P)),res_int),1,[res_int 1]);

idx2in = @(x,r) (x-.5) / (r-1) * (ima-imi) + imi;
th = idx2in(interp1q(cumsum(inva)/sum(inva),(1:res_int)',quan_int),res_int);
