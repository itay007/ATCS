function [count,p] = seqdotplot(seq1,seq2,varargin)
% SEQDOTPLOT generates a dotplot of sequence matches.
%
%   SEQDOTPLOT(S,T) plots the sequence matches of sequences S and T.
%
%   SEQDOTPLOT(S,T,WINDOW,NUM) plots sequence matches when there are at
%   least NUM matches in a window of size WINDOW. For nucleotide sequences
%   the literature recommends a WINDOW of 11 and NUM of 7.
%
%   MATCHES = SEQDOTPLOT(...) returns the number of dots in the dotplot
%   matrix.
%
%   [MATCHES, MATRIX] = SEQDOTPLOT(...) returns the dotplot as a sparse
%   matrix.
%
%   SEQDOTPLOT(...,'DISPLAY',true) turns the verbose mode on. Useful for
%   large queries where the algorithm takes large amounts of memory and
%   time.
%
%   Example:
%
%       moufflon = getgenbank('AB060288','sequence',true)
%       takin = getgenbank('AB060290','sequence',true)
%       seqdotplot(moufflon,takin,11,7)
%
%   This shows the similarities between prion protein (PrP) nucleotide
%   sequences of two ruminants, the moufflon and the golden takin.
%
%   See also ALIGNDEMO, NWALIGN, SWALIGN.

% 	Example reference:
% 	Comparative analysis of the prion protein open reading frame nucleotide
% 	sequences of two wild ruminants, the moufflon and golden takin.
%
% 	Seo SW, Hara K, Kubosaki A, Nasu Y, Nishimura T, Saeki K, Matsumoto Y,
% 	Endo H, Onodera T.

% Copyright 2002-2012 The MathWorks, Inc.

%   SEQDOTPLOT(...,'MEMORY_BLOCKSIZE',N) Set the maximum block of memory
%   used by the function.
%
%   SEQDOTPLOT(...,'SCREENSIZE',N) Overrides the screensize setting -- used
%   for testing.

%%% setting some constants
memoryBlock    = 2^25;
limitForSparse = 2^27;  % (in bytes, i.e. can contain limitForSparse/4 matches)

%%% setting defaults for parameters
verbosity = false;
window = 1;
stringency =1;

%%% get space limit for image in pixels
rootUnits = get(0, 'Units');
set(0, 'Units', 'pixels');
screnSize = get(0,'ScreenSize')*[0 0 1 0;0 0 0 1]';
set(0, 'Units', rootUnits);

%%% processing inputs
switch nargin
    case 0,
        seq1 = randseq(15);
        seq2 = seq1;
        fprintf('Creating some random data... \n')
    case 1,
        seq2 = seq1;
    otherwise
        if nargin > 2
            val = varargin{1};
            if isreal(val) && numel(val)==1 && isnumeric(val)
                window = val;
                varargin(1) = [];
            end
        end
        if nargin > 3
            val = varargin{1};
            if isreal(val) && numel(val)==1 && isnumeric(val)
                stringency = val;
                varargin(1) = [];
            else
                stringency = window;
            end
        end
end

% processing remaining varargins
if  numel(varargin)
    if rem(numel(varargin),2) == 1
        error(message('bioinfo:seqdotplot:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'display','memory_blocksize','screensize'};
    for j=1:2:numel(varargin)
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:seqdotplot:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:seqdotplot:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % verbosity on/off
                    verbosity = pval == true;
                case 2  % memoryblock
                    memoryBlock = pval;
                case 3  % screensize
                    screnSize = pval;
            end
        end
    end
end

% figure out the window limits
windowLimits = screnSize - [20 110];
windowLimits = round( windowLimits * 0.85);

% If the input is a structure then extract the Sequence data.
if isstruct(seq1)
    seq1 = bioinfoprivate.seqfromstruct(seq1);
end
if isstruct(seq2)
    seq2 = bioinfoprivate.seqfromstruct(seq2);
end

% testing valid conditions of input parameters
if stringency > window
    error(message('bioinfo:seqdotplot:DotplotTooStringent'));
end

% convert to integers and check that both sequences are of the same type
if (ischar(seq1) && ischar(seq2))
    if (bioinfoprivate.isnt(seq1) && bioinfoprivate.isnt(seq2))
        iSeq1 = nt2int(seq1);
        iSeq2 = nt2int(seq2);
    elseif (bioinfoprivate.isaa(seq1) && bioinfoprivate.isaa(seq2))
        iSeq1 = aa2int(seq1);
        iSeq2 = aa2int(seq2);
    else
        error(message('bioinfo:seqdotplot:DotplotMismatchCharSeq'));
    end
elseif (isnumeric(seq1) && isnumeric(seq2))
    if (bioinfoprivate.isnt(seq1) && bioinfoprivate.isnt(seq2))
        iSeq1 = seq1;
        iSeq2 = seq2;
    elseif (bioinfoprivate.isaa(seq1) && bioinfoprivate.isaa(seq2))
        iSeq1 = seq1;
        iSeq2 = seq2;
    else
        error(message('bioinfo:seqdotplot:DotplotMismatchNumSeq'));
    end
else
    error(message('bioinfo:seqdotplot:DotplotMismatchType'));
end

% length of input sequences
le1=length(iSeq1);
le2=length(iSeq2);

if le1>le2
    switchedSequences = true;
    tem = iSeq1;
    iSeq1 = iSeq2;
    iSeq2 = tem;
    le1=length(iSeq1);
    le2=length(iSeq2);
else
    switchedSequences = false;
end

% pad sequences to contain window at the end
padSeq = zeros(1,window-1);
iSeq1 = [iSeq1 padSeq];
iSeq2 = [iSeq2 padSeq];

mapSize = double(max([iSeq1(:);iSeq2(:)]));

% we cannot index 0's so we re-map to the last position of map
iSeq1(iSeq1 == 0) = mapSize+1;
iSeq2(iSeq2 == 0) = mapSize+1;

% create map, we do not want zeros to match zeros (last position in map)
map = diag(uint8([true(1,mapSize),false]));

columnBlock = floor(memoryBlock/le1); %including window columns

p=uint32([]);
count = 0;

firstCols =  1:columnBlock-window+1:le2;
numFirstCols = numel(firstCols);

if (numFirstCols > 1)
    % Allocate a big block of memory
    s = zeros(le1,columnBlock-window+1,'uint8');
    if verbosity
        t0 = clock;
        fprintf('Search is divided in %d memory blocks \n', numel(firstCols))
    end
else
    % allocate space for le1xle2
    s = zeros(le1,le2,'uint8');
end

for i = 1:numFirstCols
    j = firstCols(i);
    % compute some constants for this block
    lastCol = min( le2+window-1 , j+columnBlock-1);
    numCols = lastCol - j + 1 - window + 1;
    colInd  = j:lastCol;
    
    mask = map(iSeq1,iSeq2(colInd));
    
    if (numCols+window-1 == columnBlock)
        s(:)=uint8(0);
        for k = 1:window
            s = s + mask(k:le1+k-1,k:numCols+k-1);
        end
        s = s .* mask(1:le1,1:numCols);
        h = find(s >= stringency);
        count = count + numel(h);
        if (4*count) > limitForSparse
            error(message('bioinfo:seqdotplot:LimitExceeded'));
        end
        p = [p;uint32((j-1)*le1+h)]; %#ok<AGROW>
    else % when last block is incomplete
        if numCols ~= le2
            s = zeros(le1,numCols,'uint8');
        end
        for k = 1:window
            s = s + mask(k:le1+k-1,k:numCols+k-1);
        end
        s = s .* mask(1:le1,1:numCols);
        h = find(s >= stringency);
        count = count + numel(h);
        if (4*count) > limitForSparse
            error(message('bioinfo:seqdotplot:LimitExceededLast'));
        end
        p = [p;uint32((j-1)*le1+h)]; %#ok<AGROW>
    end
    if (numFirstCols > 1) && verbosity
        fprintf('Memory block: %d            Elapsed time: %f seconds \n',...
            i,etime(clock,t0))
        fprintf('Match matrix size: %d  Memory used: %d bytes\n',...
            count,count*4)
    end
end

% match matrix info should be in sparse, as it is needed fot the output,
% but also for easy downsampling computation
p=reshape(sparse(double(p),1,true,le1*le2,1),le1,le2);

% computing the scaling factors in the case the match matrix does not fit in
% the screen
downSampleX = ceil( le2 / windowLimits(1) );
downSampleY = ceil( le2 / windowLimits(2) );

if (downSampleX > 1) || (downSampleY > 1)
    
    warning(message('bioinfo:seqdotplot:imageTooBigForScreen', downSampleX, downSampleY))
    
    
    imageSizePixelsX = ceil(size(p,2)/downSampleX);
    imageSizePixelsY = ceil(size(p,1)/downSampleY);
    
    % setting the size to a multiple of downSample for easy downsampling,
    % (like padding with zeros a full matrix)
    padCordinatesX = imageSizePixelsX*downSampleX;
    padCordinatesY = imageSizePixelsY*downSampleY;
    
    if ~isequal(size(p),[padCordinatesY,padCordinatesX]) % do not pad if p
        % is already that size
        p(padCordinatesY,padCordinatesX)=false;
    end
    
    % setting the size and class of the x-downsampled matrix
    tp = sparse([],[],false,padCordinatesY,imageSizePixelsX);
    
    % x-downsampling in the sparse domain
    for k = 1:downSampleX
        tp = tp | p(:,k:downSampleX:padCordinatesX);
    end
    
    % setting the size and class of the output Image
    I = false(imageSizePixelsY,imageSizePixelsX);
    
    % y-downsampling
    I(ceil(find(tp(:))/downSampleY))=true;
    
    % resize the padded p so the output has the sequence lengths
    p = p(1:le1,1:le2);
    
else % the image fits in the screen, no need to downsample
    I =full(p);
end

if switchedSequences
    xLabel = 'Sequence 1';
    yLabel = 'Sequence 2';
else
    xLabel = 'Sequence 2';
    yLabel = 'Sequence 1';
end

hFig = figure('visible','off','tag','seqDotPlot',...
    'nextplot','replacechildren');

xData = (0:size(I,2)-1)*downSampleX+1;
yData = (0:size(I,1)-1)*downSampleY+1;
imagesc(I,'xData',xData,'yData',yData);

colormap(1-gray);

hAxis = gca;
if (downSampleX==1 && downSampleY==1)
    axis(hAxis,'image');
end
set(hAxis,'XAxisLocation','top','Units','pixels');
xlabel(xLabel,'fontsize',12);
ylabel(yLabel,'fontsize',12);

windowSize = ceil(size(I)/.85);

if any(windowSize>get(hAxis,'Position')*[0 0 0 1;0 0 1 0]')
    figPos([1 3]) = [ceil((screnSize(1)-windowSize(2))/2) windowSize(2)];
    figPos([2 4]) = [screnSize(2) - 90 - windowSize(1) windowSize(1)];
    set(hFig,'Position',figPos)
    
    axesPos([2 1]) = ceil(windowSize.*[0.05 0.1])+.5 ;
    axesPos([4 3]) = size(I);
    axesPos(2) = axesPos(2) + max(0,floor((windowSize(1)*0.85 - size(I,1))/2));
    set(hAxis,'Position',axesPos)
    set(hAxis,'Units','Normalized')
else % image fits do not need to resize
    set(hAxis,'Units','Normalized')
    pos = get(hAxis,'position');
    set(hAxis,'position',pos - [0 .05 0 0]);
    % for short sequences, force the tick marks to be integer values.
    if le1<6
        set(hAxis,'ytick',1:le1);
    end
    if le2<6
        set(hAxis,'xtick',1:le2);
    end
end

set(hFig,'Visible','on')

% with output args force everything to be doubles, clear not used outputs
switch nargout
    case 2;
        if switchedSequences
            p = double(p');
        else
            p = double(p);
        end
    case 1; clear p;
    otherwise; clear p; clear count;
end
