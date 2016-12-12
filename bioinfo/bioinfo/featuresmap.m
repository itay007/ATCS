function [legH,outFeatList] = featuresmap(s,featList,c,varargin)
%FEATURESMAP draws the map with the features of a GenBank structure.
%
% FEATURESMAP(S) creates a graphical map showing all the features in the
% GenBank structure S.
%
% FEATURESMAP(S,FEATLIST) indicates the features to plot. FEATLIST is a
% cell array of strings. Features not found in the GenBank structure are
% ignored. FEATLIST defaults to a list of all the features in the GenBank
% structure. If FEATLIST includes '-' as the first string in the cell
% array, then the remaining strings (features) are not mapped.
%
% FEATURESMAP(S,FEATLIST,LEVELS) or FEATURESMAP(S,LEVELS) indicate in which
% level of the map the features are drawn. Defaults to [1:numel(FEATLIST)].
% The total number of levels for the map is max(LEVELS).
%
% [H,OUTFEATLIST] = FEATURESMAP(...) returns a list of handles and the
% respective feature labels. Use them for putting a legend.
%
% FEATURESMAP(...,'FONTSIZE',FS) sets the font size (points) of the
% annotations of the features. FS defaults to 9 points.
%
% FEATURESMAP(...,'COLORMAP',CM) specifies a list of colors to be used for
% the features. CM may have any number of rows, but it must have exactly 3
% columns.  Each row is interpreted as a color, with the first element
% specifying the intensity of red light, the second green, and the third
% blue.  Color intensity can be specified on the interval 0.0 to 1.0.
%
% FEATURESMAP(...,'SHOWPOSITIONS',true) includes the sequence position of
% the features in the annotations. Default is false.
%
% FEATURESMAP(...,'QUALIFIERS',QUALLIST) sets an ordered list of feature
% qualifiers to search for and use as annotations. For each feature, the
% first matching qualifier found from QUALLIST is used for its annotation.
% If a feature does not include any of the qualifiers, no annotation
% displays for that feature. QUALLIST is a cell array of  strings. QUALLIST
% defaults to {'gene','product','locus_tag','note','db_xref','protein_id'}.
% Use an empty cell to indicate that no labels should be drawn.
%
% Examples:
%
%     % Create a circular map of five different features mapped on three
%     % levels and use the outputs from the featuresmap function as inputs
%     % to the legend function;
%     gbk = getgenbank('J01415')
%     [h,f] = featuresmap(gbk,{'CDS','D_loop','mRNA','tRNA','rRNA'},[1 2 2 2 3])
%     legend(h,f,'interpreter','none','location','bestoutside')
%     title('Human mitochondrion, complete genome')
%
%     % Create a linear map showing only the "gene" features, change the
%     % font of the labels to seven points and include the sequence
%     % position in the labels:
%     herpes = getgenbank('NC_001348')
%     featuresmap(herpes,{'gene'},'showpositions',true,'fontsize',7);
%     title('Genes in Human herpesvirus 3 (strain Dumas)')
%
% See also FEATURESPARSE, GENBANKREAD, GETGENBANK, SEQSTATSDEMO, SEQVIEWER.

%   Copyright 2006-2012 The MathWorks, Inc.



% =========================================================================
% == Hard-coded constants for circular map
numPoints = 2000;         % number of points to draw the complete circle
outestRingRadius = 1.25;  % width of the ring will be 0.25
imaginaryCircleForLabels = 1.5; % imaginary circle to place text


% == Hard-coded constants for both types of map
textSeparationFactor = 1;     % to adjust text overlaps, set to 1 to use
% the default extent of the text height
horizontalLineLength = 0.2;   % extension horizontal line

if ~isfield(s,'LocusSequenceLength')
    error(message('bioinfo:featuresmap:notGenbankStructure'))
end

% =========================================================================
% == Grab data from genbank structure
if ~isfield(s,'LocusTopology')
    warning(message('bioinfo:featuresmap:noTopologyFound'))
    isCircularTopology = false;
else
    isCircularTopology = strcmp(s.LocusTopology,'circular');
end
le = str2double(s.LocusSequenceLength);
f = featuresparse(s); % get features from genbank structure
clear s; % no longer needed

% =========================================================================
% == Read in madatory inputs
if  nargin == 1
    c = [];        % will force c to defaults
    featList = {}; % will force featlist to defaults
end
if nargin == 2
    if iscellstr(featList)
        c = [];    % will force c to defaults
    elseif isnumeric(featList)
        c = featList;
        featList = {};
    end
end
% check if the second input was indeed an optional input argument
if nargin>1 && ischar(featList)
    if nargin==2
        varargin = {featList};    % will error later, no paired
    else
        varargin = {featList c varargin{:}}; %#ok<CCAT> % featlist and c were optional inputs
    end
    c = [];        % will force c to defaults
    featList = {}; % will force featlist to defaults
end
% check if the third input was indeed an optional input argument
if nargin>2 && ischar(c)
    varargin = {c varargin{:}}; %#ok<CCAT>
    c = []; % will force c to defaults
end
% at this point featList must be a cell array of strings and c a numeric
% array:
if ~iscellstr(featList)
    error(message('bioinfo:featuresmap:invalidFeatlist'))
end
if ~isnumeric(c)
    error(message('bioinfo:featuresmap:invalidLevels'))
end
% check if user wants to remove features from the default feature list
if ~isempty(featList) && strcmp(featList{1},'-')
    removeFields = featList(2:end);
else
    removeFields = {};
end
% need to get the default feature list
if isempty(featList) || strcmp(featList{1},'-')
    featList = fieldnames(f);
end
% update feature list
[~,h] = setdiff(featList,removeFields);
featList = featList(sort(h));
% need to set the default levels
if isempty(c)
    c = 1:numel(featList);
end
m = max(c); % num of levels
c = c(:);
if numel(c) < numel(featList);
    l1 = m+1;
    l2 = m+numel(featList)-numel(c);
    warning(message('bioinfo:featuresmap:tooFewLevels', numel( c ), numel( featList ), l1, l2))
    c = [c;(l1:l2)'];
    m = l2;
end

% =========================================================================
% Parse optional inputs
[labelFontSize, colors, showPositionsInLabels, quallist, nquallist] = parse_inputs(varargin{:});

% replicate the number of colors to have enough for the number of features
kCol = ceil(numel(featList)/size(colors,1));
if kCol>1
    colors = repmat(colors,kCol,1);
    warning(message('bioinfo:featuresmap:tooFewColors'))
end

% create figure and axes
hf = figure('color','w','tag','Bioinfo:featuresmap');
ha = axes('Position',[0 0 1 .92],'xtick',[],'ytick',[]);
hold on; axis square; axis off;

% =========================================================================
% == Draw map

if isCircularTopology
    x = @(x,l,n) -[(1+(l-1)*(outestRingRadius-1)/n)*sin([x(1):(sign(diff(x)).*le/numPoints):x(2) x(2)]*2*pi./le) (1+(l)*(outestRingRadius-1)/n)*sin([x(2):-(sign(diff(x)).*le/numPoints):x(1) x(1)]*2*pi./le) (1+(l-1)*(outestRingRadius-1)/n)*sin(x(1)*2*pi./le)];
    y = @(x,l,n) [(1+(l-1)*(outestRingRadius-1)/n)*cos([x(1):(sign(diff(x)).*le/numPoints):x(2) x(2)]*2*pi./le) (1+(l)*(outestRingRadius-1)/n)*cos([x(2):-(sign(diff(x)).*le/numPoints):x(1) x(1)]*2*pi./le) (1+(l-1)*(outestRingRadius-1)/n)*cos(x(1)*2*pi./le)];
else
    x = @(x,l,n) (l-[1 0 0 1 1])/n;
    y = @(x,l,n) [x(1) x(1) x(2) x(2) x(1)]/le;
end

if max(c)==1
    hOutline = plot(x([0 le],25,49),y([0 le],25,49),'k','linewidth',2);
else
    hOutline = plot(x([0 le],1,1),y([0 le],1,1),'k');
end
% Disable legend for the outline
hasbehavior(hOutline,'legend',false);
j = 1;
legH = cell(numel(featList),1);
estMaxj = sum(structfun(@length,f));
labels = cell(estMaxj,1);
labelPos = zeros(1,estMaxj);
labelCol = zeros(estMaxj,3);
labelPatch = zeros(estMaxj,1);

if numel(featList)==1
    l = 1;
    legH{1} = drawPatches(featList{l},colors,c(1));
else
    for l = 1:numel(featList)
        legH{l} = drawPatches(featList{l},colors(l,:),c(l));
    end
end

lilegH = false(numel(legH),1);
for ilegH = 1:numel(legH)
    if ishghandle(legH{ilegH})
        lilegH(ilegH) = true;
    end
end
        
% outFeatList = featList(legH>0);
% legH = legH(legH>0);

outFeatList = featList(lilegH);
legH = [legH{lilegH}];

labels = labels(1:j-1);
labelPos = labelPos(1:j-1);
labelCol = labelCol(1:j-1,:);
labelPatch = labelPatch(1:j-1);

% == figure out the vertical extension of the labels
htext = text(0,0,'test','FontSize',labelFontSize);
ex = get(htext,'extent')*[0 0 0 1]';
vertext = textSeparationFactor * ex;
delete(htext)

if isCircularTopology
    alpha = labelPos.*2*pi/le;
    margin = false; % not implemented
    done = false;
    while ~done
        [hl,ht] = circularannotations(outestRingRadius,imaginaryCircleForLabels,...
            vertext,labelFontSize,horizontalLineLength,...
            margin,alpha,labels,labelPatch);
        if numel(hl)==1
            d = get(hl,'Ydata');
        else
            d = cell2mat(get(hl,'Ydata'));
        end
        if max(-min(d(:)),max(d(:))) > imaginaryCircleForLabels*1.5
            estimatedReductionFactor = imaginaryCircleForLabels*1.5./max(-min(d(:)),max(d(:)));
            vertext = vertext * min(0.9,estimatedReductionFactor);
            delete(hl);
            delete(ht);
        else
            done = true;
        end
    end
    if numel(d)
        y1 = min(d(:))-ex/2;
        y2 = max(d(:))+ex/2;
        x = (y2-y1)*560/420/0.92/2;
        axis([-x x y1 y2])
    end
    set(ha,'DataAspectRatio',[1 1 1])
    % place seqlength
    text(-.05,.8,'1 bp','horizontal','right','vertical','top','fontsize',labelFontSize,'clipping','on')
    text(.05,.8,sprintf('%d bp',le),'horizontal','left','vertical','top','fontsize',labelFontSize,'clipping','on')
else
    pos = labelPos/le;
    done = false;
    while ~done
        [hl,ht] = linearannotations(1,2,vertext,labelFontSize,...
            horizontalLineLength,pos,labels,labelPatch);
        if numel(hl)==1
            d = get(hl,'Ydata');
        else
            d = cell2mat(get(hl,'Ydata'));
        end
        if max(-min(d(:)),max(d(:)))>1.2
            estimatedReductionFactor = 1.2./max(-min(d(:)),max(d(:)));
            vertext = vertext * min(0.9,estimatedReductionFactor);
            delete(hl);
            delete(ht);
        else
            done = true;
        end
    end
    % and adjust axes
    if numel(d)
        axis([-3 9 min(-.1,min(d(:))-ex/2) max(1.1,max(d(:))+ex/2)])
    else
        axis([-6 6 -.1 1.1])
    end
    % place seqlength
    text(-.5,0,'1 bp','horizontal','right','vertical','bottom','fontsize',labelFontSize,'clipping','on')
    text(-.5,1,sprintf('%d bp',le),'horizontal','right','vertical','middle','fontsize',labelFontSize,'clipping','on')
end

%coloring lines & link to respective patches
for count = 1:numel(hl)
    set(hl(count),'color',labelCol(count,:),'UserData',labelPatch(count))
    set(ht(count),'UserData',labelPatch(count))
    % Disable legend for the annotation lines
    hasbehavior(hl(count),'legend',false);
end

set([ht;hl],'ButtonDownFcn', @(h,e) placeDataTip(get(h,'UserData')))

for jLeg = 1:numel(legH)
    set(legH(jLeg),'DisplayName',outFeatList{jLeg});
end

if nargout == 0
    clear legH;
end

% fix to make sure zoom works appropriately, we could remove it after 2006b
set(ha,'PlotBoxAspectRatioMode','Auto')

% prepare data tip and its callbacks
hTip = text(0,0,'','BackgroundColor', [1 1 0.933333],'Color', [0 0 0],...
    'EdgeColor', [0.8 0.8 0.8],'Visible','off',...
    'Tag','Bioinfo:featuresmap:dataTip','interpreter','none',...
    'fontsize',9,'fontname','Helvetica','horizontal','left');
set(hf,'WindowButtonUpFcn',@(h,e) deleteDataTip);
set(findobj(hf,'Type','Patch'),'ButtonDownFcn', @(h,e) placeDataTip(h))




% =========================================================================
% == Nested Functions
    function legH = drawPatches(field,color,level)
        legH = NaN;
        if (level>0) && isfield(f,field)
            feats = f.(field);
            fnames = fieldnames(feats);
            for i = 1:numel(feats)
                for k = 1:2:numel(feats(i).Indices)
                    patchColor = color(rem(j,size(color,1))+1,:);
                    legH = patch(x(feats(i).Indices(k+[0 1])+[0 1],level,m),y(feats(i).Indices(k+[0 1])+[0 1],level,m),patchColor,'UserData',[l,i,level]);
                    % disable legend. Later re-enable for the final patch
                    hasbehavior(legH,'legend',false);
                    % try to find a qualifier to use it as label for this feature
                    done = false; p = 0; testLabel = '';
                    while ~done && p<nquallist
                        p = p+1;
                        annotFieldInd = strmatch(quallist{p},fnames,'exact');
                        if ~isempty(annotFieldInd)
                            testLabel = feats(i).(fnames{annotFieldInd});
                            if ~isempty(testLabel)
                                % mark segmented features
                                if numel(feats(i).Indices)>2
                                    testLabel = [testLabel ' (' mat2str((k+1)/2) ')']; %#ok<AGROW>
                                end
                                done = true;
                            end
                        end
                    end
                    if showPositionsInLabels
                        testLabel = [testLabel ' [' mat2str(feats(i).Indices(k)) ':' mat2str(feats(i).Indices(k+1)) ']']; %#ok<AGROW>
                    end
                    if ~isempty(testLabel)
                        if iscell(testLabel)
                            labels{j} = testLabel{1};
                        else
                            labels{j} = testLabel;
                        end
                        labelPos(j) = mean(feats(i).Indices(k+[0 1]));
                        labelCol(j,:) = patchColor;
                        labelPatch(j) = legH;
                        j = j + 1;
                    end
                end
            end
        end
        % Re-enable legend for the final patch
        if ishghandle(legH)
           hasbehavior(legH,'legend',true);
        end
    end % drawPatches nested function
% =========================================================================
    function placeDataTip(hp)
        % hp - handle to patch
        in = get(hp,'UserData'); % which feature and which index (pre-stored in userdata)
        set(hTip,'String',featurestruct2str(f.(featList{in(1)})(in(2)),featList{in(1)},50),'visible','on')
        set(hf,'WindowButtonMotionFcn',@(h,e) moveDataTip)
        moveDataTip
    end % placeDataTip nested function
% =========================================================================
    function moveDataTip
        set(hTip,'Position',[1 1]*get(ha,'CurrentPoint')/2)
    end % moveDataTip nested function
% =========================================================================
    function deleteDataTip
        set(hf,'WindowButtonMotionFcn',[])
        set(hTip,'String','','visible','off')
    end % deleteDataTip nested function
% =========================================================================
end % featuresmap main function

% =========================================================================
% == Sub-Functions
function [hl,ht] = linearannotations(l1,l2,vertext,labelFontSize,horline,pos,labels,hp) %#ok
% hl,ht   - handles to the lines and text objects
% l1,l2   - internal and external limits
% vertext - vertical extent of the labels
% labelFontSize - font size
% horline - length of the horizontal part of the line to the annotation
% pos     - origin points for every label
% label   - cell string with the annotations
% hp      - handles to respective patches

% calculate the default position for every label
lx = repmat(l2,size(pos));
ly = pos;

% adjust the vertical position of the labels
[u,h] = sort(ly);
ly(h) = spaceAdjust(u,vertext);

% place annotations
ht = text(lx+(lx>=0)*horline*2-horline,ly,labels,'fontsize',labelFontSize,...
    'interpreter','none','Vertical','middle','Horizontal','left',...
    'clipping','on');

% draw lines to annotations
hl = line([lx*0+l1;lx;lx+(lx>=0)*horline*2-horline],[pos;ly;ly]);

end % linearannotations function
% =========================================================================
function [hl,ht] = circularannotations(r1,r2,vertext,labelFontSize,horline,margin,alpha,labels,hp) %#ok
% hl,ht   - handles to the lines and text objects
% r1,r2   - internal and external radious
% vertext - vertical extent of the labels
% labelFontSize - font size
% horline - length of the horizontal part of the line to the annotation
% margin  - extend horizontal line to margins ? (not implemented yet)
% alpha   - origin angles for every label
% label   - cell string with the annotations
% hp      - handles to respective patches

% some hard coded constants
adjustTextAtBottomTopOfCircularMaps = 0.5 * r1;

% calculate the default position for every label
lx = -r2*sin(alpha);
ly = r2*cos(alpha);

% adjust the vertical position of the right labels
h1 = find(lx>=0);
[u,h2] = sort(ly(h1));
ly(h1(h2)) = spaceAdjust(u,vertext);

% adjust the vertical position of the left labels
h1=find(lx<0);
[u,h2]=sort(ly(h1));
ly(h1(h2))=spaceAdjust(u,vertext);

% shifts the horizontal position to r2 at least
tmp = sqrt(r2.^2 - ly.^2);
tmp((tmp.*tmp)<0) = eps * 1000;
lx = sign(lx).* tmp;

% tries to give some horizontal shifts to the labels at the bottom and top
% to avoid the overlapping of lines
h = lx>=0 & ly<0;
h2 = lx<adjustTextAtBottomTopOfCircularMaps & h;
if any(h2)
    h3 = find(h2);
    [dummy,h4] = sort(ly(h3));
    p = min(lx(h&~h2))/sum(h2); if isempty(p); p = r1; end
    lx(h3(h4)) = (0:sum(h2)-1)*p;
end
h = lx<0 & ly<0;
h2 = lx>-adjustTextAtBottomTopOfCircularMaps & h;
if any(h2)
    h3 = find(h2);
    [dummy,h4] = sort(ly(h3)); 
    p = max(lx(h&~h2))/sum(h2); if isempty(p); p = -r1; end
    lx(h3(h4)) = -eps * 1000 + (0:sum(h2)-1)*p;
end
h = lx>=0 & ly>=0;
h2 = lx<adjustTextAtBottomTopOfCircularMaps & h;
if any(h2)
    h3 = find(h2);
    [dummy,h4] = sort(ly(h3),'descend'); 
    p = min(lx(h&~h2))/sum(h2); if isempty(p); p = r1; end
    lx(h3(h4)) = (0:sum(h2)-1)*p;
end
h = lx<0 & ly>=0;
h2 = lx>-adjustTextAtBottomTopOfCircularMaps & h;
if any(h2)
    h3 = find(h2);
    [dummy,h4] = sort(ly(h3),'descend'); 
    p = max(lx(h&~h2))/sum(h2); if isempty(p); p = -r1; end
    lx(h3(h4)) = -eps * 1000 + (0:sum(h2)-1)*p;
end

lxmz = lx>=0;
% a final check to avoid reversed angles of the annotation lines
lx(lxmz) = max(lx(lxmz),-r1*sin(alpha(lxmz)));
lx(~lxmz) = min(lx(~lxmz),-r1*sin(alpha(~lxmz)));

% place annotations
lxt = lx+(lxmz)*horline*2-horline;
ht = text(lxt,ly,labels,'fontsize',labelFontSize,'interpreter','none',...
    'Vertical','middle','Horizontal','left','clipping','on');
set(ht(~lxmz),'Horizontal','right')

% draw lines to annotations
hl = line([-r1*sin(alpha);lx;lxt],[r1*cos(alpha);ly;ly]);

end % circularannotations function
% =========================================================================
function w = spaceAdjust(u,v)
% u - goal (should be ordered)
% v - minimum separation between output points
% w - output

eps1000 = 1000*eps;
n = numel(u);

% initialize output
w = u;
for i = 2:n
    w(i) = max(w(i-1)+v,u(i));
end

moving = true;
while moving % main loop
    moving = false;
    i=1;
    while i<=n % move block by block
        b = find(abs(diff(w(i:n))-v)>eps1000,1)-1; % finds next block
        if isempty(b); b = n-i; end   % gets last block
        sh = mean(u(i:i+b)-w(i:i+b)); % estimate shift
        if abs(sh)>eps1000
            if i==1;   leftLim  =-inf; else leftLim  = w(i-1)+v;   end
            if i+b==n; rightLim = inf; else rightLim = w(i+b+1)-v; end
            if w(i)+sh   < leftLim;  sh = leftLim-w(i);    end
            if w(i+b)+sh > rightLim; sh = rightLim-w(i+b); end
            w(i:i+b) = w(i:i+b) + sh;
            moving = true;
        end
        i = i + b + 1; % next block
    end
    i=1;
    while ~isempty(i) % move singles
        k0 = abs(diff(w)-v)>eps1000;
        for i = find(([1 k0] & w>u) | ([k0 1] & w<u))
            if i==1; leftLim  =-inf; else leftLim  = w(i-1)+v; end
            if i==n; rightLim = inf; else rightLim = w(i+1)-v; end
            w(i) = max(min(u(i),rightLim),leftLim);
            moving = true;
        end
    end
end % moving
end % spaceAdjust function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [labelFontSize, colors, showPositionsInLabels, quallist, nquallist] = parse_inputs(varargin)
% Parse the varargin parameter/value inputs
% Check that we have the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:featuresmap:IncorrectNumberOfArguments', mfilename));
end

% The allowed inputs

okargs = {'fontsize','colormap','showpositions','qualifiers'};
% set defaults of optional input arguments
quallist = {'gene','product','locus_tag','note','db_xref','protein_id'};
nquallist = numel(quallist);
colors = [0 0 1;0 .7 0;1 0 0;0 1 1;1 0 1;1 1 0;.6 0 0;.4 1 .4;1 .5 0;.7 0 1;.6 .6 .4;.8 .8 .8];
showPositionsInLabels = false;
labelFontSize = 9;

for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1 % 'fontsize'
            if isscalar(pval) && isnumeric(pval) && pval>0
                labelFontSize = pval;
            else
                error(message('bioinfo:featuresmap:InvalidFontSize'))
            end
        case 2 % 'colormap'
            if size(pval,2)==3 && size(pval,1)>0 && max(pval(:))<=1 && min(pval(:))>=0
                colors = pval;
            else
                error(message('bioinfo:featuresmap:InvalidColorMap'))
            end
        case 3 % 'showpositions'
            showPositionsInLabels = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 4 % 'qualifiers'
            quallist = pval(:);
            if ~iscellstr(quallist)
                error(message('bioinfo:featuresmap:invalidQuallist'))
            end
            nquallist = numel(quallist);
    end
end
end
