function showalignment(alignment,varargin)
%SHOWALIGNMENT colored display of an alignment.
%
%   SHOWALIGNMENT(ALIGNMENT) displays an alignment in the MATLAB figure
%   window. For pairwise alignments matches and similar residues are
%   highlighted and ALIGNMENT is the output from one of the functions
%   NWALIGN or SWALIGN. For multiple sequence alignment highly conserved
%   columns are highlighted and ALIGNMENT is the output from the function
%   MULTIALIGN.
%
%   SHOWALIGNMENT(...,'MATCHCOLOR',COLOR) selects the color used to
%   highlight the matches or the totally conserved columns in the output
%   display. The default color is red.
%
%   SHOWALIGNMENT(...,'SIMILARCOLOR',COLOR) selects the color used to
%   highlight similar residues that are not exact matches or highly
%   conserved columns. The default color is magenta.
%
%   COLOR can be a 1x3 RGB vector whose elements specify the intensities of
%   the red, green and blue component of the color; the intensities can be
%   in the range [0 1]. COLOR can also be a character from the following
%   list:
%
%            'b'     blue
%            'g'     green
%            'r'     red
%            'c'     cyan
%            'm'     magenta
%            'y'     yellow
%
%   The following options are only available when showing pairwise
%   alignments:
%
%   SHOWALIGNMENT(...,'STARTPOINTERS',STARTAT) specifies the starting
%   indices in the original sequences of a local alignment. STARTAT is the
%   two element vector returned as the third output of the SWALIGN
%   function.
%
%   SHOWALIGNMENT(...,'COLUMNS',COLS) specifies how many columns per line to
%   use in the output. The default is 64.
%
%   SHOWALIGNMENT(...,'TERMINALGAP',TF) excludes terminal gaps from the
%   count of matches and similar residues when TF is FALSE. Default is
%   to count terminal gaps (TRUE). This option is valid for pairwise
%   alignments only.
%
%       Examples:
%
%       % Display an alignment between two amino acid sequences:
%       [score, alignment] = nwalign('VSPAGMASGYD','IPGKASYD');
%       showalignment(alignment)
%
%       % Display a multiple alignment:
%       gag = multialignread('aagag.aln');
%       showalignment(gag)
%
%       % Note: To view a multiple sequence alignment and also interact
%       % with it, use the SEQALIGNVIEWER function.
%
%   See also ALIGNDEMO, MULTIALIGN, NWALIGN, SEQALIGNVIEWER, SWALIGN.

%   Copyright 2002-2012 The MathWorks, Inc.

bioinfochecknargin(nargin,1,mfilename);

%== figure out which type of alignment we have:
% 1. a pairwise alignment can be a char array with two or three rows, if
%    three rows then the middle row must contain pipes, spaces and colons
% 2. a multiple alignment can be a char array of 3 or more rows or a
%    structure with the field Sequence. Field Header or Name are optional,
%    if found they are used to label the sequences.
seqNames = []; % by default no sequence names are given
if isstruct(alignment) && isfield(alignment,'Sequence')
    %try to get headers
    if isfield(alignment,'Header')
        seqNames = char(alignment(:).Header);
    elseif isfield(alignment,'Name')
        seqNames = char(alignment(:).Name);
    elseif isfield(alignment,'LocusName')
        seqNames = char(alignment(:).LocusName);
    end
    alignment = {alignment(:).Sequence};
    if numel(unique(cellfun('length',alignment)))==1
        alignment = char(alignment);
    else
        error(message('bioinfo:showalignment:SequencesNotAligned'))
    end
end

if ischar(alignment)
    [numRows,alignmentLen] = size(alignment);
    alignment = upper(alignment);
    
    if numRows == 2
        % create our own match string
        matchString = blanks(alignmentLen);
        matches = (alignment(1,:) == alignment(2,:));
        similar = matches & false;
        matchString(matches) = '|';
        alignment = [alignment(1,:); matchString; alignment(2,:)];
        isMultipleAlignment = false;
        
    elseif numRows == 3
        
        if all(ismember(alignment(2,:),'|: '))
            matches = ( alignment(1,:) == alignment(3,:));
            if ~all((alignment(2,:) == '|') == matches)
                warning(message('bioinfo:showalignment:InconsistentAlignment'));
            end
            similar = (alignment(2,:) == ':');
            isMultipleAlignment = false;
        else
            isMultipleAlignment = true;
        end
    else
        isMultipleAlignment = true;
    end
else
    error(message('bioinfo:showalignment:BadAlignmentFormat'))
end

if numRows < 2
    error(message('bioinfo:showalignment:TooFewSequences'))
end

% in case of multiple alignment we still need to figure out the highly
% conserved columns and if header exists we need to add the to the
% displayed char array
if isMultipleAlignment
    % try to guess which type of alphabet we have
    if bioinfoprivate.isnt(alignment)
        alpha = 'nt';
    elseif bioinfoprivate.isaa(alignment)
        alpha = 'aa';
    else
        error(message('bioinfo:showalignment:InvalidSymbolsInSequences'))
    end
    [con,score] = seqconsensus(alignment,'ambiguous','count','alpha',alpha,'gaps','all');
    alignment = [seqNames repmat(' ',numRows,1) alignment];
    con = [repmat(' ',1,size(seqNames,2)+1) con];
    score = [repmat(inf,1,size(seqNames,2)+1) score];
    matches = ~score;
    switch alpha
        case 'aa'
            similar = (score < 7) & ~matches & con~='-';
        case 'nt'
            similar = (score < 5) & ~matches & con~='-';
    end
end

% now deal with any options
if isMultipleAlignment
    wrap = size(alignment,2);
else
    wrap = 64;
end

color =  'ff0000';
simcolor = 'ff00ff';
startat = [1;1];

noDisplay = false;  % use this for testing
terminalgap = true;

if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:showalignment:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'matchcolor','columns','nodisplay',...
        'similarcolor','startpointers','terminalgap'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:showalignment:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:showalignment:AmbiguousParameterName', pname));
        else
            switch(k)
                
                case 1 %match color
                    color = setcolorpref(pval, color);
                case 2% wrap
                    if ~isMultipleAlignment
                        wrap = pval;
                    else
                        warning(message('bioinfo:showalignment:MultiAlignColumns'));
                    end
                case 3% noDisplay
                    noDisplay = pval;
                case 4 %similar color
                    simcolor = setcolorpref(pval, simcolor);
                case 5 %startat
                    startat = pval(:);
                    if ~isnumeric(startat) || numel(startat) > 2
                        error(message('bioinfo:showalignment:BadStartPointers'));
                    elseif numel(startat) == 1
                        startat = [startat;startat]; %#ok
                    end
                    if isMultipleAlignment
                        warning(message('bioinfo:showalignment:MultiAlignPointers'));
                    end
                case 6 %terminalgap
                    terminalgap = bioinfoprivate.opttf(pval);
            end
        end
    end
end

if ~noDisplay
    import com.mathworks.mwswing.MJScrollPane;
    import java.awt.Color;
    import java.awt.Dimension;
    import com.mathworks.toolbox.bioinfo.sequence.*;
    
    % Create the viewer
    b = awtcreate('com.mathworks.toolbox.bioinfo.sequence.ShowLocalAlignment');
    
    % Set Java color for match and similar
    b.changeMatchColor(Color(hex2dec(color(1:2))/255, hex2dec(color(3:4))/255, hex2dec(color(5:6))/255));
    b.changeSimilarColor(Color(hex2dec(simcolor(1:2))/255, hex2dec(simcolor(3:4))/255, hex2dec(simcolor(5:6))/255));
    
    % show
    if isMultipleAlignment
        b.displayAlignment(alignment, matches, similar, wrap);
    else
        if terminalgap
            count = numel(matches);
        else
            mask = all(~(alignment([1 3],:)=='-' | alignment([1 3],:)==' '));
            count = find(mask,1,'last')-find(mask,1)+1;
        end
        b.displayAlignment(alignment, matches, similar, count, startat, wrap);
    end
    
    % Setup a scrollpane and put b into the scrollPane
    scrollpanel = MJScrollPane(b, MJScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, MJScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
    
    % Create a figure
    fTitle = sprintf('Aligned Sequences');
    hFigure = figure( ...
        'WindowStyle', 'normal', ...
        'Menubar', 'none', ...
        'Toolbar', 'none', ...
        'NumberTitle','off',...
        'Tag', 'seqshoworfs',...
        'Resize', 'on', ...
        'Name', fTitle,...
        'HandleVisibility', 'Callback',...
        'DeleteFcn',{@deleteView, b});
    % Set the figure widow size to fit the scrollPane
    d = awtinvoke(scrollpanel, 'getPreferredSize()');
    pos = getpixelposition(hFigure);
    if ~isMultipleAlignment
        pos(3) = d.getWidth;
    end
    pos(4) = d.getHeight;
    setpixelposition(hFigure,pos);
    
    figurePosition = get(hFigure, 'Position');
    [viewP, viewC] = javacomponent( scrollpanel, ...
        [0, 0, figurePosition(3), figurePosition(4)], ...
        hFigure);
    set(viewC, 'units', 'normalized');
    set(hFigure, 'userdata', viewC);
    % Get toolbox/matlab/icon path
    iconPath = fullfile(toolboxdir('matlab'),'icons');
    
    tb = uitoolbar(hFigure);
    cicon= load(fullfile(iconPath,'printdoc.mat')); % load cdata of print icon from toolbox/matlab/icon/printdoc.mat
    a1=uipushtool(tb, ...
        'CData', cicon.cdata, ...
        'TooltipString', 'Print',...
        'ClickedCallback', {@print_cb, b}); %#ok
end

%----------------------------------------------------
function print_cb(hSrv, event, b) %#ok
awtinvoke(b, 'printAlignment()');

%----------------------------------------------
function deleteView(hfig, event, b)%#ok
if ~isempty(b)
    awtinvoke(b, 'clearPrintView()');
end
viewC = get(hfig, 'userdata');
if ~isempty(viewC)
    delete(viewC)
end
