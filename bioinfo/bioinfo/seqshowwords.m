function wordOutput =seqshowwords(seq,word,varargin)
%SEQSHOWWORDS displays a sequence with selected words highlighted.
%
%   SEQSHOWWORDS(SEQ,WORD) displays the sequence, SEQ, in the MATLAB figure
%   window, with all occurrences of the word WORD highlighted. SEQSHOWWORDS
%   returns a structure of the start and stop positions of all occurrences
%   of WORD in SEQ. WORD can be a regular expression. WORD can also be a
%   cell array of words.
%
%   SEQSHOWWORDS(...,'COLOR',COLOR) selects the color used to highlight the
%   words in the output display. The default color is red. COLOR can be a
%   1x3 RGB vector whose elements specify the intensities of the red, green
%   and blue component of the color; the intensities can be in the range [0
%   1]. COLOR can also be a character from the following list:
%
%            'b'     blue
%            'g'     green
%            'r'     red
%            'c'     cyan
%            'm'     magenta
%            'y'     yellow
%
%   SEQSHOWWORDS(...,'COLUMNS',COLS) specifies how many columns per line
%   to use in the output. The default is 64.
%
%   SEQSHOWWORDS(...,'ALPHABET',A) specifies that SEQ and WORD are
%   amino acids ('AA') or nucleotides ('NT'). The default is NT.
%
%   SEQSHOWWORDS(...,'EXACT',true) only finds exact matches for ambiguous
%   or extended nucleotide or amino acid symbols. For example K will only
%   match K and not match G and T.
%
%   Notes: 
%       If the search word or words contain amino acid or nucleotide
%       symbols that represent multiple symbols, then seqshowwords shows
%       all possible matches. For example, the symbol R represents either G
%       or A (purines). If Word is 'ART', then seqshowwords shows
%       occurrences of both 'AAT' and 'AGT'.
%
%       SEQSHOWWORDS does not highlight overlapping patterns multiple times.
%
%   Examples:
%
%       seqshowwords('GCTATAACGTATATATATA','TATA');
%
%       % This highlights two places, the first occurrence of 'TATA' and
%       % the 'TATATATA' immediately after 'CG'. The final 'TA' is not
%       % highlighted because the preceding 'TA' is part of an already
%       % matched pattern. To highlight all multiple repeats of TA, use the
%       % regular expression 'TA(TA)*TA'.
%
%       seqshowwords('GCTATAACGTATATATATA','TA(TA)*TA');
% 
%       % Show multiple words
%       seqshowwords('GCTATAACGTATATATATA',{'CG', 'GC'});
%
%   See also CLEAVE, PALINDROMES, REGEXP, RESTRICT, SEQDISP, SEQMATCH,
%   SEQSHOWORFS, SEQVIEWER, SEQWORDCOUNT, STRFIND.

%   Copyright 2002-2012 The MathWorks, Inc.

% If the input is a structure then extract the Sequence data.
if isstruct(seq)
    seq = bioinfoprivate.seqfromstruct(seq);
end

wrap = 64;
origword = word;
useRegexp = true;
alphabet = {};

color = 'FF0000';
noDisplay = false;
if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:seqshowwords:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'color','columns','exact','nodisplay','alphabet'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname))); 
        if isempty(k)
            error(message('bioinfo:seqshowwords:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:seqshowwords:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 %color
                    color = setcolorpref(pval, color);
                case 2% wrap
                    wrap = pval;
                case 3% no regexp
                    useRegexp = ~bioinfoprivate.opttf(pval);
                    if isempty(useRegexp)
                        error(message('bioinfo:seqshowwords:InputOptionNotLogical', okargs{ k }));
                    end
                case 4% noDisplay
                    noDisplay = pval;
                case 5 % alphabet
                    alphabet = pval;
                    if strcmpi(pval,'aa')
                        alphabet = 'aa';
                    end
            end
        end
    end
end

if (~iscell(word))
    if useRegexp && ((bioinfoprivate.isnt(word) && bioinfoprivate.isnt(seq)) || (bioinfoprivate.isaa(word) &&bioinfoprivate.isaa(seq)))
        if isempty(alphabet)
            word = seq2regexp(word);
        else
            word = seq2regexp(word,'ALPHABET',alphabet);
        end
    end
else
%     Remove empty cells
    empidx = arrayfun(@(x) isempty(x{:}), word);
    word(empidx) = [];
    for i = 1:length(word)
        if useRegexp && ((bioinfoprivate.isnt(word{i}) && bioinfoprivate.isnt(seq)) || (bioinfoprivate.isaa(word{i}) &&bioinfoprivate.isaa(seq)))
            if isempty(alphabet)
                word{i} = seq2regexp(word{i});
            else
                word{i} = seq2regexp(word{i},'ALPHABET',alphabet);
            end
        end
    end
end

seqLen = length(seq);
if seqLen  > 500000  % large display will mess up the help browser's memory
    error(message('bioinfo:seqshowwords:ShowWordsLimit'))
end

% regexp struggles with more than 4 or 5 [ ] groups

%numGroups = sum(word == '[');
%if seqLen > 20000 && numGroups > 4
%    warning('bioinfo:seqshowwords:RegexpMayTakeTime','This search may take a long time...');
%end

% % [starts,stops] = regexpi(seq,word);
if iscellstr(word)
    [~, idxW] = unique(word, 'first');
    [starts,stops] = regexpi(seq,word(sort(idxW)));
    
    starts2 = [];
    stops2 = [];
    for i = 1:length(starts)
        starts2 = [starts2, starts{i}]; %#ok
        stops2 = [stops2, stops{i}]; %#ok
    end
    starts = starts2;
    stops = stops2;
    
    [starts, sidx] = sort(starts);
    stops = stops(sidx);
    
    zeroidx = find(diff(starts) == 0);
    while ~isempty(zeroidx)
        diff_stops = diff(stops);
        negidx = find(diff_stops(zeroidx)<= 0);
        if ~isempty(negidx)
            starts(zeroidx(negidx)+1) = [];
            stops(zeroidx(negidx)+1) = [];
        end
        zeroidx = find(diff(starts) == 0);
        
        diff_stops = diff(stops);
        posidx = find(diff_stops(zeroidx)> 0);
        if ~isempty(posidx)
            stops(zeroidx(posidx)) = stops(zeroidx(posidx)+1);
            starts(zeroidx(posidx)+1) = [];
            stops(zeroidx(posidx)+1) = [];
        end
        
        zeroidx = find(diff(starts) == 0);
    end
    
    w = 1;
    while w <= numel(starts)-1
        i = w+1;
        while i <= numel(starts)
            del_s = stops(w) - starts(i);
            del_e = stops(w) - stops(i);
            
            if del_s >= 0
                if del_e <= 0
                    stops(w) = stops(i);
                end
                starts(i) = [];
                stops(i) = [];
                w = 1;
                i = w+1;
            else
                i = i+1;
            end
        end
        w = w + 1;
    end
else
    [starts,stops] = regexpi(seq,word);
end

% save the output
wordOutput.Start = starts;
wordOutput.Stop = stops;

if ~noDisplay
    showwords(seq, starts, stops, origword, wrap, color)
end
% end % end of main function

function showwords(str, startidx, stoptidx, realword, width, colour)
% Display words in sequence in figure window

import com.mathworks.mwswing.MJScrollPane;
import java.awt.Color;
import java.awt.Dimension;
import com.mathworks.toolbox.bioinfo.sequence.*;

% Create Java color
c = Color(hex2dec(colour(1:2))/255, hex2dec(colour(3:4))/255, hex2dec(colour(5:6))/255);

% Create the viewer
b = SequenceShowWords(c);

if isempty(startidx)
    startidx=length(str)+1;
end

startidx = startidx-1; % Java start with 0
% show
b.displayWords(str, startidx, stoptidx, width);

% Setup a scrollpane and put b into the scrollPane
scrollpanel = MJScrollPane(b, MJScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, MJScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);

% Create a figure
if ischar(realword)
    fTitle = sprintf('Occurrences of %s', realword);
elseif iscell(realword)
    tstr = bioma.util.catCellStrToStr(realword);
    if length(tstr) > 64
        tstr = [tstr(1:64) '...'];
    end
    fTitle = sprintf('Occurrences of the words: %s', tstr);
end
hFigure = figure( ...
    'WindowStyle', 'normal', ...
    'Menubar', 'none', ...
    'Toolbar', 'none', ...
    'NumberTitle','off',...
    'Tag', 'seqshowwords',...
    'Resize', 'on', ...
    'HandleVisibility', 'Callback',...
    'Name', fTitle,...
    'DeleteFcn',{@deleteView, b});
% Set the figure widow size to fit the scrollPane
d = awtinvoke(scrollpanel, 'getPreferredSize()');
pos = getpixelposition(hFigure);
pos(3) = d.getWidth;
pos(4) = d.getHeight;      
setpixelposition(hFigure,pos);

figurePosition = get(hFigure, 'Position');
[viewP, viewC] = javacomponent( scrollpanel, ...
    [0, 0, figurePosition(3), figurePosition(4)], ...
    hFigure);

set(viewC, 'units', 'normalized');
 set(hFigure, 'userdata', viewC);
 
%  Add custom toolbar
% Get toolbox/matlab/icon path
iconPath = fullfile(toolboxdir('matlab'),'icons');

tb = uitoolbar(hFigure);
cicon= load(fullfile(iconPath,'printdoc.mat')); % load cdata of print icon from toolbox/matlab/icon/printdoc.mat
a1=uipushtool(tb, ...
    'CData', cicon.cdata,...
    'TooltipString', 'Print',...
    'ClickedCallback', {@print_cb, b}); %#ok

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function print_cb(hSrv, event, b) %#ok
awtinvoke(b, 'printWords()');

function deleteView(hfig, event, b)%#ok
awtinvoke(b, 'clearPrintView()');
viewC = get(hfig, 'userdata');
if ~isempty(viewC)
    delete(viewC)
end


