function hp = msdotplot(fh,p,varargin)
%MSDOTPLOT generates a dotplot of a LCMS (GCMS) dataset
%
%  MSDOTPLOT(P,T) plots the multiple scans of a LCMS (GCMS) dataset. P is a
%  cell array of peak lists, where each element is a two-column matrix with
%  m/z values in the first column and ion intensity values in the second
%  column. Each element corresponds to a spectrum or retention time. T is a
%  vector of retention times associated with an LCMS or GCMS data set. The
%  number of elements in T equals the number of elements in the cell array
%  P. MSDOTPLOT uses the current figure if it is a heatmap generated with
%  MSHEATMAP, otherwise it creates a new figure.
%
%  MSDOTPLOT(FH,P,T) plots into the axes cointained by the figure with
%  handle FH.
%
%  MSDOTPLOT(...,'QUANTILEVALUE',Q) determines which peaks are displayed in
%  the dotplot. Choices are a scalar between 0 and 1. Default is 0. For
%  example, when Q = 0 all the peaks are plotted, and when Q = 0.5 only the
%  approximately half of the peaks are plotted, those with the larger ion
%  intensity value.
%
%  PH = MSDOTPLOT(...) returns a handle to the lineseries object. You can
%  use this handle as input to the get function to display a list of the
%  plot's properties. You can use this handle as input to the set function
%  to change the plot's properties, including showing and hiding points.
%
%  Examples:
%
%    load lcmsdata
%
%    % Create a dot plot with only the 5% largest peaks:
%    msdotplot(ms_peaks,ret_time,'Quantile',0.95)
%
%    % Draw all the peak information over a heatmap of the LCMS dataset:
%    [MZ,Y] = msppresample(ms_peaks,5000);
%    msheatmap(MZ,ret_time,log(Y))
%    msdotplot(ms_peaks,ret_time)
%    axis([480 532 375 485]) % zoom in to see the detail
%
%  See also DIFFPROTDEMO, LCMSDEMO, MSHEATMAP, MSPALIGN, MSPPRESAMPLE,
%  MZXML2PEAKS, MZXMLREAD.

%   Copyright 2006-2008 The MathWorks, Inc.


% check inputs
bioinfochecknargin(nargin,2,mfilename);

if numel(fh)>1 % MSDOTPLOT(P,T)
    t = p(:);
    p = fh;
    if isempty(get(0,'children'))
        createNewFigure = true;
    elseif isempty(findobj(gcf,'Tag','MSHeatMap'))
        createNewFigure = true;
    else
        createNewFigure = false;
        fh = gcf;
    end
else % MSDOTPLOT(FH,p,t)
    t = varargin{1}(:);
    varargin(1)=[];
    createNewFigure = false;
    if isempty(findobj(fh,'Tag','MSHeatMap'))
        error(message('bioinfo:msdotplot:invalidHeatMap'))
    end
end

%set defaults
Q = 0;
useMemSafeQuantile = false;

if  numel(varargin) > 0
    if rem(numel(varargin),2) == 1
        error(message('bioinfo:msdotplot:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'quantilevalue','memsafequantile'};
    for j=1:2:numel(varargin)
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:msdotplot:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:msdotplot:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % quantile
                    if ~isnumeric(pval) || ~isscalar(pval) || pval<0 || pval>=1
                        error(message('bioinfo:msdotplot:badLevels'))
                    end
                    Q = pval;
                case 2 % use memory safe quantile estimate
                    useMemSafeQuantile = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            end
        end
    end
end

if createNewFigure
    figure;
    if ~useMemSafeQuantile
        try
            le = cumsum(cellfun(@(x) size(x,1),p));
            h = zeros(le(end),1);h(le(1:end-1)+1)=1;h = cumsum(h)+1;
            K = single([t(h) ,cell2mat(p)]);
            K(K(:,3)<=quantile(K(:,3),Q),:)=[];
        catch outerException
            try % try the memory safe estimate
                K = memSafeQuantile(p,Q,t);
            catch innerException
                rethrow(outerException)
            end
        end
    else
        K = memSafeQuantile(p,Q,t);
    end
    hp = plot3(K(:,2),K(:,1),K(:,3),'k.','DisplayName','Centroids','Tag','msdotplot');
    set(hp,'MarkerSize',1)
    xlabel('Mass/Charge (M/Z)')
    ylabel('Retention Time')
    axis([min(K(:,2)) max(K(:,2)) min(K(:,1)) max(K(:,1))])
else % work over the existing heatmap
    figure(fh)
    
    if ~ishold
        hold on
        holdStatus = 'off';
    else
        holdStatus = 'on';
    end
    
    le = cumsum(cellfun(@(x) size(x,1),p));
    h = zeros(le(end),1);h(le(1:end-1)+1)=1;h = cumsum(h)+1;
    
    if numel(getappdata(fh,'msheatmapListeners'))==3
        t = (1:numel(t))';
    end
    
    T = t(h);
    K = single([T ,cell2mat(p)]);
    K(K(:,3)<=quantile(K(:,3),Q),:)=[];
    
    % get the range of the heatmap image for lcms data and remove anything
    % outside of that range. Note that K(:,2) are the MZ values
    figAppdata = getappdata(fh);
    if isfield(figAppdata,'msheatmapMZRange')
        K((K(:,2) < figAppdata.msheatmapMZRange(1)),:) = [];
        K((K(:,2) > figAppdata.msheatmapMZRange(2)),:) = [];
    end
    
    hp = plot3(K(:,2),K(:,1),K(:,3),'k.','DisplayName','Centroids','Tag','msdotplot');
    
    set(hp,'MarkerSize',1)
    hold(holdStatus)
end

setappdata(hp,'legend_hgbehavior',0)
setAllowAxesRotate(rotate3d(gcf),gca,false)
msDataCursor(gcf)

if nargout==0
    clear hp
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vals = memSafeQuantile(p,Q,t)
% MEMSAFEQUANTILE small memory footprint guess at quantile

% find the range of the data
le = (cellfun(@(x) size(x,1),p));
maxVals = cellfun(@(x) max(x(:,2)),p);
minVals = cellfun(@(x) min(x(:,2)),p);
totalPoints = sum(le);

totalMax = max(maxVals(~isinf(maxVals)));
totalMin = min(minVals(~isinf(minVals)));

% deal with the case where we still have Inf values
if isempty(totalMax)
    totalMax = -inf;
    for i=1:numel(p)
        infMask = isinf(p{i}(:,2));
        localMax = max(p{i}(~infMask,2));
        if localMax>totalMax
            totalMax = localMax;
        end
        
    end
end

if isempty(totalMin)
    totalMin = inf;
    for i=1:numel(p)
        infMask = isinf(p{i}(:,2));
        localMin = min(p{i}(~infMask,2));
        if localMin<totalMin
            totalMin = localMin;
        end
    end
end

% figure out how many points we are looking for

qPoints = ceil((1-Q)*totalPoints);

% guess an initial value -- midpoint of the range is probably good enough
fraction = .5;
guess = totalMin + fraction*(totalMax-totalMin);
count = totalPoints;
oldCount = 0;
iter = 0;
% now count how many points are greater than the guess and then bisect the
% range until we converge.
while (count~=qPoints) && (count~=oldCount && (totalMax-totalMin)>1) && iter<1000
    oldCount = count;
    count = sum(sum(cellfun(@(x) sum(x(:,2) >= guess),p)));
    if count == qPoints
        break
    elseif count<qPoints
        totalMax = guess;
        guess = totalMin + fraction*(guess-totalMin);
    else
        totalMin = guess;
        guess = totalMax - fraction*(totalMax-guess);
    end
    % provide an escape just in case we get stuck
    iter = iter+1;
end

% now create the output structure
best = guess;
vals = zeros(count,3,class(p{1}));
currentPos = 1;
for i=1:numel(p)
    ndx = p{i}(:,2)>= best;
    numFound = sum(ndx);
    vals(currentPos:currentPos+numFound-1,2) = p{i}(ndx,1);
    vals(currentPos:currentPos+numFound-1,3) = p{i}(ndx,2);
    vals(currentPos:currentPos+numFound-1,1) = t(i);
    currentPos = currentPos+numFound;
end
vals(currentPos:end,:)=[];
