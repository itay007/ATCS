function [prof,h1,h2] = profalign(prof1,prof2,varargin)
%PROFALIGN performs Needleman-Wunsch global alignment of two profiles.
%
%  PROF = PROFALIGN(PROF1,PROF2) returns a new profile PROF for the optimal
%  global alignment of two profiles. The profiles PROF1 and PROF2 are
%  numeric arrays of size [(4 or 5 or 20 or 21) x Profile Length] with
%  counts or weighted profiles. Weighted profiles are used to down-weight 
%  similar sequences and up-weight divergent sequences. The output profile
%  will be a numeric matrix of size [(5 or 21) x New Profile Length], the
%  last row represents gaps. Original gaps in the input profiles are
%  preserved. The output profile is the result of adding the aligned
%  columns of the input profiles. 
%
%  [PROF, H1, H2] = PROFALIGN(PROF1,PROF2) returns pointers that indicate
%  how to rearrange the columns of the original profiles into the new
%  profile. 
%
%  PROFALIGN(...,'SCORINGMATRIX',SM) defines the scoring matrix SM to be
%  used for the alignment. The default is BLOSUM50 for AA or NUC44 for NT.
%
%  PROFALIGN(...,'GAPOPEN',{G1,G2}) sets the penalties for opening a gap in
%  the first and second profiles respectively. G1 and G2 can be either 
%  scalars or vectors. When a vector is used the number of elements is one
%  more than the length of the input profile; every element indicates the
%  position specific penalty for opening a gap between two consecutive
%  symbols in the sequence, being the first and the last element the gap
%  penalties employed at the flanks of the sequence. The default gap open
%  penalties are {10,10}.
%
%  PROFALIGN(...,'EXTENDGAP',{E1,E2}) sets the penalties for extending a
%  gap in the first and second profile respectively. E1 and E2 can be
%  either scalars or vectors. When a vector is used the number of elements 
%  is one more than the length of the input profile; every element
%  indicates the position specific penalty for extending a gap between two
%  consecutive symbols in the sequence, being the first and the last
%  element the gap penalties employed at the flanks of the sequence. If
%  EXTENDGAP is not specified, then extensions to gaps are scored with the
%  same value as GAPOPEN. 
%
%  PROFALIGN(...,'EXISTINGGAPADJUST',false) turns off the automatic
%  adjustment, based on existing gaps, of the position specific penalties
%  for opening a gap. By default (true), for every profile position
%  PROFALIGN will lower proportionally the penalty for opening a gap
%  towards the penalty of extending a gap based on the proportion of gaps
%  found in the contiguous symbols and on the weight of the input profile.
%
%  PROFALIGN(...,'TERMINALGAPADJUST',true) adjusts the penalty for opening
%  a gap at the ends equal to the penalty for extending a gap. Default is
%  false.
%
%  PROFALIGN(...,'SHOWSCORE',true) displays the scoring space and the
%  winning path.
%
%  Examples: 
%
%      ma1 = ['RGTANCDMQDA';'RGTAHCDMQDA';'RRRAPCDL-DA'];
%      ma2 = ['RGTHCDLADAT';'RGTACDMADAA'];
%      p1  = seqprofile(ma1,'gaps','all','counts',true);
%      p2  = seqprofile(ma2,'counts',true);
%
%      % Merge two profiles into a single one by aligning them.
%      p = profalign(p1,p2);
%      seqlogo(p)
%
%      % Using the output pointers to generate the multiple alignment.
%      [p, h1, h2] = profalign(p1,p2);
%      ma = repmat('-',5,12);
%      ma(1:3,h1) = ma1;
%      ma(4:5,h2) = ma2;
%      disp(ma)
%
%      % Up-weighting the gap penalty before Cysteine in the second profile.
%      gapVec = 10 + [p2(aa2int('C'),:) 0] * 10
%      p3 = profalign(p1,p2,'gapopen',{10,gapVec});
%      seqlogo(p3)
% 
%      % Add a new sequence to a profile without inserting new gaps into
%      % the profile.
%      gapVec = [0 inf(1,11) 0];
%      p4 = profalign(p3,seqprofile('PLHFMSVLWDVQQWP'),'gapopen',{gapVec,10});
%      seqlogo(p4)
%             
%  See also  HMMPROFALIGN, MULTIALIGN, NWALIGN, SEQPROFILE, SEQCONSENSUS.

% References:
%   J.D. Thompson, D.G. Higgins, and T.J. Gibson. Nucleic Acids Res. (1994)
%   22(22):4673-4680.
%   R. Durbin, S. Eddy, A. Krogh, and G. Mitchison. Biological Sequence
%   Analysis. Cambridge UP, 1998.
%   Needleman, S. B., Wunsch, C. D., J. Mol. Biol. (1970) 48:443-453.
%
% Undocumented constants to score the existing gap-gap and gap-residue
% cases that exist already in the profiles:
%
%  PROFALIGN(...,'GAPGAPSCORE',GGS) defines the score for matching an
%  existing gap in one profile to an existing gap in the other profile. GGS
%  can be a scalar or a function specified using @. PROFALIGN passes four
%  values to the function: the average score for two matched residues (sm),
%  the average score for two mismatched residues (sx), and, the length of
%  both profiles or sequences (len1 and len2). GGS defaults to
%  @(sm,sx,len1,len2) 0.1*sm. 
%
%  PROFALIGN(...,'GAPRESSCORE',GRS) defines the score for matching an
%  existing gap in one profile to a residue in the other profile. GRS can
%  be a scalar or a function specified using @. PROFALIGN passes four
%  values to the function: the average score for two matched residues (sm), 
%  the average score for two mismatched residues (sx), and, the length of
%  both profiles or sequences (len1 and len2). GRS defaults to
%  @(sm,sx,len1,len2) 0.1*sx. 
%

%   Copyright 2005-2012 The MathWorks, Inc.

% Defaults
setGapExtend = false;
oldGapAdjust = true;
endGapAdjust = false;
showscore = false;

% Undocumented constants to score the existing gap-gap and gap-residue
% cases that exist already in the profiles, sm and sx are the mean of the
% scores for matches and mismatches respectively in the score matrix
matchGapGapFunction = @(sm,sx,len1,len2) 0.1*sm;
matchGapGapIsFunctionHandle = true;
matchGapResFunction = @(sm,sx,len1,len2) 0.1*sx;
matchGapResIsFunctionHandle = true;

[numSym1,len1] = size(prof1);
[numSym2,len2] = size(prof2);

% Default gap penalty vectors
go1(1:len1+1) = -8;
go2(1:len2+1) = -8;

% Obtain profile weights
wep1 = mean(sum(prof1));
wep2 = mean(sum(prof2));

if ~len1||~len2
    error(message('bioinfo:profalign:InvalidProfileLengths'));
end

% add row for gaps (if necessary)
if any(numSym1 == [4 20])
    prof1(end+1,:) = 0;
    numSym1 = numSym1+1;
elseif  all(numSym1 ~= [5 21])
     error(message('bioinfo:profalign:invalidProfile1'))
end
if any(numSym2 == [4 20])
    prof2(end+1,:) = 0;
    numSym2 = numSym2+1;
elseif  all(numSym2 ~= [5 21])
     error(message('bioinfo:profalign:invalidProfile2'))
end
% check size of input profiles
if numSym1==numSym2
    numSym = numSym1;
    isAminoAcid = (numSym==21);
else
    error(message('bioinfo:profalign:differentProfileType'))
end

% Check input arguments
if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:profalign:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'scoringmatrix','gapopen','extendgap','showscore',...
              'existinggapadjust','terminalgapadjust','gapgapscore',...
              'gapresscore'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:profalign:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:profalign:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % scoring matrix
                    if isnumeric(pval)
                        ScoringMatrix = pval;
                    else
                        if ischar(pval)
                            pval = lower(pval);
                        end
                        try
                            ScoringMatrix = feval(pval);
                        catch allExceptions
                            error(message('bioinfo:profalign:InvalidScoringMatrix'));
                        end
                    end
                case 2 %gap open penalty
                    if iscell(pval)
                        go1 = -pval{1};
                        go2 = -pval{2};
                        if isscalar(go1)
                            go1(1:len1+1) =  go1;
                        else
                            go1 = go1(:)';
                        end
                        if isscalar(go2)
                            go2(1:len2+1) =  go2;
                        else
                            go2 = go2(:)';
                        end
                    elseif isscalar(pval)
                        go1(1:len1+1) = -pval; 
                        go2(1:len2+1) = -pval;
                    else
                        error(message('bioinfo:profalign:InvalidGapOpen'))
                    end
                    if length(go1)~=(len1+1) || length(go2)~=(len2+1)
                        error(message('bioinfo:profalign:InvalidLengthGapOpen'))
                    end
                case 3 %gap extend penalty
                    if iscell(pval)
                        ge1 = -pval{1};
                        ge2 = -pval{2};
                        if isscalar(ge1)
                            ge1(1:len1+1) =  ge1;
                        else
                            ge1 = ge1(:)';
                        end
                        if isscalar(ge2)
                            ge2(1:len2+1) =  ge2;
                        else
                            ge2 = ge2(:)';
                        end
                    elseif isscalar(pval)
                        ge1(1:len1+1) = -pval; 
                        ge2(1:len2+1) = -pval;
                    else
                        error(message('bioinfo:profalign:InvalidExtendGap'))
                    end
                    if length(ge1)~=(len1+1) || length(ge2)~=(len2+1)
                        error(message('bioinfo:profalign:InvalidLengthExtendGap'))
                    end
                    setGapExtend = true;
                case 4 % showscore
                    showscore = bioinfoprivate.opttf(pval);
                    if isempty(showscore)
                        error(message('bioinfo:profalign:showscoreInputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 5 % oldgapadjust
                    oldGapAdjust = bioinfoprivate.opttf(pval);
                    if isempty(oldGapAdjust)
                        error(message('bioinfo:profalign:oldgapadjustInputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 6 % terminalgapadjust  
                    endGapAdjust = bioinfoprivate.opttf(pval);
                    if isempty(endGapAdjust)
                        error(message('bioinfo:profalign:endgapadjustInputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                % Undocumented constants to score the existing gap-gap and gap-residue
                % cases that exist already in the profiles:
                case 7 % gapgapscore
                    if isscalar(pval) && isnumeric(pval)
                        matchGapGapIsFunctionHandle = false;
                        matchGapGapConstant = pval;
                    elseif isa (pval, 'function_handle')
                        matchGapGapFunction = pval;
                    else
                        error(message('bioinfo:profalign:InvalidGapGapScore'))
                    end
                case 8 % gapresscore
                    if isscalar(pval) && isnumeric(pval)
                        matchGapResIsFunctionHandle = false;
                        matchGapResConstant = pval;
                    elseif isa (pval, 'function_handle')
                        matchGapResFunction = pval;
                    else
                        error(message('bioinfo:profalign:InvalidGapResScore'))
                    end
            end
        end
    end
end

% setting the default scoring matrix
if ~exist('ScoringMatrix','var')
    if isAminoAcid
        ScoringMatrix = blosum50;
    else
        ScoringMatrix = nuc44;
    end
end

[smR,smC] = size(ScoringMatrix);

% introduce the gap-gap and gap-residue into the scoring matrix
if isAminoAcid
    if smR<20 || smC<20 || smR~=smC
        error(message('bioinfo:profalign:invalidAAScoringMatrix'))
    else
        SM = zeros(21,21);
        SM(1:20,1:20) = ScoringMatrix(1:20,1:20);
        if matchGapGapIsFunctionHandle || matchGapResIsFunctionHandle
            sm = mean(diag(SM));
            sx = sum(sum(SM-diag(diag(SM))))/380;
            if matchGapResIsFunctionHandle
                matchGapResConstant = matchGapResFunction(sm,sx,len1,len2);
            end
            if matchGapGapIsFunctionHandle
                matchGapGapConstant = matchGapGapFunction(sm,sx,len1,len2);
            end
        end
        SM(21,:) = matchGapResConstant;
        SM(:,21) = matchGapResConstant;
        SM(21,21) = matchGapGapConstant;
    end
else % for nucleotides
    if smR<4 || smC<4 || smR~=smC
        error(message('bioinfo:profalign:invalidNTScoringMatrix'))
    else
        SM = zeros(5,5);
        SM(1:4,1:4) = ScoringMatrix(1:4,1:4);
        if matchGapGapIsFunctionHandle || matchGapResIsFunctionHandle
            sm = mean(diag(SM));
            sx = sum(sum(SM-diag(diag(SM))))/12;
            if matchGapResIsFunctionHandle
               matchGapResConstant = matchGapResFunction(sm,sx,len1,len2);
            end
            if matchGapGapIsFunctionHandle
               matchGapGapConstant = matchGapGapFunction(sm,sx,len1,len2);
            end
        end
        SM(5,:) = matchGapResConstant;
        SM(:,5) = matchGapResConstant;
        SM(5,5) = matchGapGapConstant;
    end
end

% rescale gap penalties to account for sequence weights or counts
if oldGapAdjust
    wg1 = sum(prof1(1:end-1,:));
    wg2 = sum(prof2(1:end-1,:));
    go1 = go1 .* min([wg1 inf;inf wg1]);
    go2 = go2 .* min([wg2 inf;inf wg2]);
else
    wg1 = repmat(wep1,1,len1);    
    wg2 = repmat(wep2,1,len2);    
    go1 = go1 .* wep1;
    go2 = go2 .* wep2;
end

if setGapExtend
    ge1 = ge1 * wep1;
    ge2 = ge2 * wep2;
    % make sure that gap open penalties were not scaled to less that the
    % gap extension penalties otherwise the resulting alignment is
    % unbalanced
    if any(go1>ge1) 
        go1 = min(go1,ge1);
    end
    if any(go2>ge2)
        go2 = min(go2,ge2);
    end
end

if endGapAdjust % adjust terminal penalties
    go1(1)=ge1(1);
    go2(1)=ge2(2);
    go1(end)=ge1(end);
    go2(end)=ge2(end);
end

% do the alignment
if setGapExtend  
    [F, pointer] = affinegap(prof1,len1,prof2,len2,SM,go1,go2,ge1,ge2,wg1,wg2);
else
    [F, pointer] = simplegap(prof1,len1,prof2,len2,SM,go1,go2,wg1,wg2);
end

% trace back through the pointer matrix
half2=ceil((len2+1)/2);
half1=ceil((len1+1)/2);

i = len2+1; j = len1+1;
path = zeros(len2+len1,2);
step = 1;

score =  max(F(len2+1,len1+1,:));

if setGapExtend

    if F(len2+1,len1+1,3)==score      % favor with left-gap
        laststate=3;
    elseif F(len2+1,len1+1,2)==score  % then with up-gap
        laststate=2;
    else                              % at last with match
        laststate=1;
    end

    while (i>half2 && j>half1) % in the rigth half favor gaps when several
        % paths lead to the highest score
        state=laststate;
        if bitget(pointer(i,j,state),3)
            laststate=3;
        elseif bitget(pointer(i,j,state),2)
            laststate=2;
        else
            laststate=1;
        end

        switch state
            case 1 % is diagonal
                j = j - 1;
                i = i - 1;
                path(step,:) = [j,i];
            case 2 % is up
                i = i - 1;
                path(step,2) = i;
            case 3 % is left
                j = j - 1;
                path(step,1) = j;
        end
        step = step + 1;
    end

    while (i>1 || j>1)          % in the rigth half favor matchs when several
        % paths lead to the highest score
        state=laststate;
        if bitget(pointer(i,j,state),1)
            laststate=1;
        elseif bitget(pointer(i,j,state),2)
            laststate=2;
        else
            laststate=3;
        end

        switch state
            case 1 % is diagonal
                j = j - 1;
                i = i - 1;
                path(step,:) = [j,i];
            case 2 % is up
                i = i - 1;
                path(step,2) = i;
            case 3 % is left
                j = j - 1;
                path(step,1) = j;
        end
        step = step + 1;
    end

else % ~setGapExtend
    while (i > 1 || j > 1)
        switch pointer(i,j)
            case 1 % diagonal only
                j = j - 1;
                i = i - 1;
                path(step,:) = [j,i];
            case 2 % up only
                i = i - 1;
                path(step,2) = i;
            case 4 % left only
                j = j - 1;
                path(step,1) = j;
            case 6 % up or left --> up (favors gaps in seq2)
                j = j - 1;
                path(step,1) = j;
            otherwise %3 diagonal or up         --> diagonal (favors no gaps)
                %4 diagonal or left       --> diagonal (favors no gaps)
                %7 diagonal or left or up --> diagonal (favors no gaps)
                j = j - 1;
                i = i - 1;
                path(step,:) = [j,i];
        end
        step = step + 1;
    end
end % if setGapExtend

path(step:end,:) = []; % step-1 is the length of the new profile
path = flipud(path);

%setting the size of the output profile
prof = zeros(numSym,step-1);

mask1 = path(:,1)>0;
mask2 = path(:,2)>0;

% adding the aligned profiles
prof(:,mask1) = prof1;
prof(:,mask2) = prof(:,mask2)+prof2;

% updating to gaps 
prof(numSym,~mask1) = prof(numSym,~mask1)+wep1;
prof(numSym,~mask2) = prof(numSym,~mask2)+wep2;

h1=find(mask1);
h2=find(mask2);
    
if showscore
    figure
    F=max(F(2:end,2:end,:),[],3);
    clim=max(max(max(abs(F(~isinf(F))))),eps);
    imagesc(F,[-clim clim]);
    colormap(privateColorMap(1));
    set(colorbar,'YLim',[min([F(:);-eps]) max([F(:);eps])])
    title('Score for best path')
    xlabel('Profile 1')
    ylabel('Profile 2')
    hold on
    plot(path(all(path>0,2),1),path(all(path>0,2),2),'k.')
end

%-------------------------------------------------------------------------%
function [F, pointer] = simplegap(p1,m,p2,n,SM,g1,g2,wg1,wg2)
% Standard Needleman-Wunsch algorithm

% p1,p2 input profiles
% m,n profile lengths
% SM scoring matrix
% g1,g2 gap penalties
% wg1,wg2 weights for gaps

% set up storage for dynamic programming matrix
F = zeros(n+1,m+1);
F(2:end,1) = g1(1) * wg2(1) * (1:n)';
F(1,2:end) = g2(1) * wg1(1) * (1:m);

% and for the back tracing matrix
pointer= repmat(uint8(4),n+1,m+1);
pointer(:,1) = 2;  % up
pointer(1,1) = 1;  

% initialize buffers to the first column
ptr = pointer(:,2); % ptr(1) is always 4
currentFColumn = F(:,1);

% Score matches
SC = p2' * SM * p1;

% main loop runs through the matrix looking for maximal scores
for outer = 2:m+1

    % grab the data from the matrices and initialize some values
    lastFColumn    = currentFColumn;
    currentFColumn = F(:,outer);
    best = currentFColumn(1);
    
    for inner = 2:n+1
        % score the three options
        up       = best + g1(outer)*wg2(inner-1);
        left     = lastFColumn(inner) + g2(inner)*wg1(outer-1);
        diagonal = lastFColumn(inner-1) +SC(inner-1,outer-1);

        % max could be used here but it is quicker to use if statements
        if up > left
            best = up;
            pos = 2;
        else
            best = left;
            pos = 4;
        end

        if diagonal >= best
            best = diagonal;
            ptr(inner) = 1;
        else
            ptr(inner) = pos;
        end
        currentFColumn(inner) = best;

    end % inner
    % put back updated columns
    F(:,outer)   = currentFColumn;
    % save columns of pointers
    pointer(:,outer)  = ptr;
end % outer

%-------------------------------------------------------------------------%
function [F,pointer] = affinegap(p1,m,p2,n,SM,go1,go2,ge1,ge2,wg1,wg2)
% Needleman-Wunsch algorithm modified to handle affine gaps

% p1,p2 input profiles
% m,n profile lengths
% SM scoring matrix
% go1,go2 gap open penalties
% ge1,ge2 gap extension penalties
% wg1,wg2 weights for gaps

% Set states
inAlign =   1;
inGapUp =   2;
inGapLeft = 3;
numStates = 3;

% Set up storage for dynamic programming matrix:
% for keeping the maximum scores for every state
F =  zeros(n+1,m+1,numStates);
F(:,1,:) = -inf;
F(1,:,:) = -inf;
F(1,1,inAlign) = 0;

F(2:end,1,inGapUp)   = (go1(1) + ge1(1) * (0:n-1)') * wg2(1);
F(1,2:end,inGapLeft) = (go2(1) + ge2(1) * (0:m-1)) * wg1(1);

% and for the back tracing pointers
pointer(n+1,m+1,numStates) = uint8(0);
pointer(2:end,1,inGapUp)   = 2;  % up
pointer(1,2:end,inGapLeft) = 4;  % left

% initialize buffers to the first column
ptrA = pointer(:,1,inAlign);
ptrU = pointer(:,1,inGapLeft);
ptrL = pointer(:,1,inGapUp);

currentFColumnA = F(:,1,inAlign);
currentFColumnU = F(:,1,inGapUp);
currentFColumnL = F(:,1,inGapLeft);

% Score matches
SC = p2' * SM * p1;

% main loop runs through the matrix looking for maximal scores
for outer = 2:m+1
    
    % grab the data from the matrices and initialize some values for the
    % first row the most orderly possible
    lastFColumnA    = currentFColumnA;
    currentFColumnA = F(:,outer,inAlign);
    bestA           = currentFColumnA(1);
    currentinA      = lastFColumnA(1);
    
    lastFColumnU    = currentFColumnU;
    currentFColumnU = F(:,outer,inGapUp);
    bestU           = currentFColumnU(1);
    
    lastFColumnL    = currentFColumnL;
    currentFColumnL = F(:,outer,inGapLeft);
    currentinGL     = lastFColumnL(1);
    
    for inner = 2:n+1
        
        % grab the data from the columns the most orderly possible
        upOpen      = bestA + go1(outer)*wg2(inner-1); 
        inA         = currentinA;
        currentinA  = lastFColumnA(inner);
        leftOpen    = currentinA + go2(inner)*wg1(outer-1);
        
        inGL        = currentinGL;
        currentinGL = lastFColumnL(inner);
        leftExtend  = currentinGL + ge2(inner)*wg1(outer-1);
        
        upExtend = bestU + ge1(outer)*wg2(inner-1); 
        inGU     = lastFColumnU(inner-1);
        
         % operate state 'inGapUp'
            
        if upOpen > upExtend
            bestU = upOpen; ptr = 1;   % diagonal
        elseif upOpen < upExtend
            bestU = upExtend; ptr = 2; % up
        else % upOpen == upExtend
            bestU = upOpen; ptr = 3;   % diagonal and up
        end
        currentFColumnU(inner)=bestU;
        ptrU(inner)=ptr;

        % operate state 'inGapLeft'
        
        if leftOpen > leftExtend
            bestL = leftOpen; ptr = 1;   % diagonal
        elseif leftOpen < leftExtend
            bestL = leftExtend; ptr = 4; % left
        else % leftOpen == leftExtend
            bestL = leftOpen; ptr = 5;   % diagonal and left
        end
        currentFColumnL(inner) = bestL;
        ptrL(inner) = ptr;

        % operate state 'inAlign'
                
        if  inA > inGU
            if inA > inGL
                bestA = inA; ptr = 1;  % diagonal
            elseif inGL > inA
                bestA = inGL; ptr = 4; % left
            else
                bestA = inA; ptr = 5;  % diagonal and left
            end
        elseif inGU > inA
            if inGU > inGL
                bestA = inGU; ptr = 2; % up
            elseif inGL > inGU
                bestA = inGL; ptr = 4; % left
            else
                bestA = inGU; ptr = 6; % up & left
            end
        else
            if inA > inGL
                bestA = inA; ptr = 3;  % diagonal & up
            elseif inGL > inA
                bestA = inGL; ptr = 4; % left
            else
                bestA = inA; ptr = 7;  % all
            end
        end
        bestA = SC(inner-1,outer-1) + bestA;
        currentFColumnA(inner) = bestA;
        ptrA(inner) = ptr;
                
    end %inner
    
    % put back updated columns
    F(:,outer,inGapLeft) = currentFColumnL;
    F(:,outer,inGapUp)   = currentFColumnU;
    F(:,outer,inAlign)   = currentFColumnA;
    % save columns of pointers
    pointer(:,outer,inAlign)  = ptrA;
    pointer(:,outer,inGapUp)  = ptrU;
    pointer(:,outer,inGapLeft)= ptrL;
end %outer

%-------------------------------------------------------------------------%
function pcmap = privateColorMap(selection)
%PRIVATECOLORMAP returns a custom color map
switch selection
    case 1, pts = [0 0 .3 20;
            0 .1 .8 25;
            0 .9 .5 15;
            .9 1 .9 8;
            1 1 0 26;
            1 0 0 26;
            .4 0 0 0];
    otherwise, pts = [0 0 0 128; 1 1 1 0];
end
xcl(cumsum(pts(1:end-1,4))+1) = 1;
pcmap = interp1(pts(:,1:3),cumsum([1;1./pts(cumsum(xcl(1:end-1))+1,4)]));
