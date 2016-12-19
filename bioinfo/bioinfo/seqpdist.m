function D = seqpdist(seqs,varargin)
%SEQPDIST computes pairwise distance between sequences
%
%  D = SEQPDIST(SEQS) returns a vector D containing biological distances
%  between each pair of sequences stored in the M elements of the cell
%  SEQS. D is an 1-by-(M*(M-1)/2) vector, corresponding to the M*(M-1)/2
%  pairs of sequences in SEQS. The output D is arranged in the order of
%  ((2,1),(3,1),..., (M,1),(3,2),...(M,2),.....(M,M-1), i.e., the lower
%  left triangle of the full M-by-M distance matrix.  To get the distance
%  between the Ith and Jth sequences (I > J), use the formula
%  D((J-1)*(M-J/2)+I-J). SEQS can also be a vector of structures with the
%  field 'Sequence' or a matrix of chars.
%
%  SEQPDIST(...,'METHOD',method) selects the method used to compute the
%  distances between every pair of sequences. Choices are: 
%
%  Distances defined for both nucleotides and amino acids:
%
%     'p-distance'      -  proportion of sites at which the two sequences
%                          are different. p --> 1 for poorly related and 
%                          p --> 0 for similar sequences.
%     'Jukes-Cantor'    -  maximum likelihood estimate of the number of
%     (default)            of substitutions between two sequences:
%                          for NT  d = -3/4 * log(1 - p * 4/3)
%                          for AA  d = -19/20 * log(1 - p * 20/19)
%     'alignment-score' -  distance (d) between two sequences (1 and 2) is
%                          computed from the pairwise alignment score (s)
%                          as follows: 
%                          d(1,2) = (1-s(1,2)/s(1,1)) * (1-s(1,2)/s(2,2))
%                          This option does not imply that pre-aligned
%                          input sequences will be realigned, it only
%                          scores them. Use with care; this distance method
%                          does not comply with the ultrametric condition.
%                          * In the rare case s(x,y)>s(x,x) then d(x,y)=0
%
%  Distances defined only for nucleotides and no scoring of gaps:
%
%     'Tajima-Nei'      -  maximum likelihood estimate considering the
%                          background nucleotide frequencies. It can be
%                          computed from the input sequences or given by
%                          setting 'OPTARGS' to [gA gC gG gT].
%     'Kimura'          -  considers separately the transitional and
%                          transversional nucleotide substitutions.
%     'Tamura'          -  considers separately the transitional and
%                          transversional nucleotide substitutions and the
%                          GC content. GC content can be computed from the
%                          input sequences or given by setting 'OPTARGS'.
%     'Hasegawa'        -  considers separately the transitional and
%                          transversional nucleotide substitutions and the
%                          background nucleotide frequencies.
%                          Background frequencies can be computed from the
%                          input sequences or given by setting 'OPTARGS' to
%                          [gA gC gG gT].
%     'Nei-Tamura'      -  considers separately the transitional
%                          substitutions between purines, the transitional
%                          substitutions between pyrimidines and the
%                          transversional substitutions and the background
%                          nucleotide frequencies. Background frequencies
%                          can be computed from the input sequences or
%                          given by setting 'OPTARGS' to [gA gC gG gT].
%
%  Distances defined only for amino acids and no scoring of gaps:
%    
%     'Poisson'         -  assumes that the number of amino acid
%                          substitutions at each site has a Poisson
%                          distribution.
%     'Gamma'           -  assumes that the number of amino acid
%                          substitutions at each site has a Gamma
%                          distribution with parameter 'a'. 'a' can be set
%                          by 'OPTARGS'. Default value is 2.
%
%  A user defined distance function may also be specified using @. For
%  example @DISTFUN, the distance function must be of the form:
%
%  function D = DISTFUN(S1, S2, OPTARGS)
%
%  taking as arguments two same length sequences (NT or AA) plus zero or
%  more additional problem-dependent arguments in OPTARGS, and returning a
%  scalar which represents the distance between S1 and S2.
%
%  SEQPDIST(...,'ALPHABET',alpha) specifies whether the sequences are
%  amino acids ('AA') or nucleotides ('NT'). The default is AA.
%
%  SEQPDIST(...,'INDELS',indels) indicates how to treat sites with gaps.
%  Options are: 
%
%     'score'           -  scores these sites either as a point mutation or
%     (default)            with the alignment parameters depending on the
%                          method selected. 
%     'pairwise-delete' -  for every pairwise comparison it ignores the
%                          sites with gaps.
%     'complete-delete' -  ignores all the columns in the multiple
%                          alignment which contain a gap, this option is
%                          available only if a multiple alignment was
%                          provided at the input SEQS.
%
%  SEQPDIST(...,'OPTARGS',optargs) some distance methods require or accept
%  optional arguments. Use a cell array to pass more than one input argument
%  (e.g., The nucleotide frequencies in Tajima-Nei distance function can be
%  specified instead of computing them from the input sequences).
%
%  SEQPDIST(...,'PAIRWISEALIGNMENT',true) forces pairwise alignment,
%  ignoring the multiple alignment of the input sequences (if any). If
%  input sequences do not have the same length and the flag is not set, a
%  pairwise alignment is performed by default. Pairwise alignment may be
%  slow when working with a large number of sequences. The default is false.
%
%  SEQPDIST(...,'USEPARALLEL',true) computes the pairwise alignments in
%  parallel on separate MATLAB workers. In order to execute the alignments
%  in parallel you must have the Parallel Computing Toolbox installed and a
%  Parallel Pool available. If a Parallel Pool does not exist then one is
%  automatically created when the auto-creation option is enabled in your
%  parallel preferences. Computation runs in serial mode without issuing
%  any warning otherwise. Defaults to false, or serial computation mode.
%
%  SEQPDIST(...,'SQUAREFORM',true) converts the output into a square format,
%  so that D(I,J) denotes the distance between the Ith and the Jth
%  sequence. The output matrix is symmetric and has a zero diagonal.
%
%  The remaining input parameters are analogous to NWALIGN and are used
%  when PAIRWISEALIGNMENT is true or METHOD is 'alignment-score'. Get help
%  on NWALIGN for more info.
%
%   SEQPDIST(...,'SCORINGMATRIX',matrix) 
%   SEQPDIST(...,'SCALE',scale) 
%   SEQPDIST(...,'GAPOPEN',penalty)
%   SEQPDIST(...,'EXTENDGAP',penalty) 
%
%  Examples:
%
%     % Load a multiple alignment of amino acids:
%     seqs = fastaread('pf00002.fa');
%
%     % For every possible pair of sequences in the multiple alignment
%     % removes sites with gaps and scores with the substitution matrix
%     % PAM250:
%
%     dist = seqpdist(seqs,'method','alignment-score',...
%                  'indels','pairwise-delete','scoringmatrix',@pam250)
%
%     % To force the re-alignment of every pair of sequences ignoring the
%     % provided multiple alignment:
%
%     dist = seqpdist(seqs,'method','alignment-score',...
%                   'indels','pairwise-delete','scoringmatrix',@pam250,...
%                   'pairwisealignment',true)
%
%     % To measure the 'Jukes-Cantor' pairwise distances after re-aligning
%     % every pair of sequences, counting the gaps as point mutations:
%
%     dist = seqpdist(seqs,'method','jukes-cantor','indels','score',...
%                  'scoringmatrix',@pam250,'pairwisealignment',true)
%
%  See also BIRDFLUDEMO, DNDS, DNDSML, FASTAREAD, PHYTREE, PHYTREE/PDIST,
%  SEQLINKAGE, SEQNEIGHJOIN.

% Copyright 2003-2012 The MathWorks, Inc.


%-% set defaults
pairWiseAlignment = false;
useParallel = false;
squareFormOut = false;
showDots = false;
isAmino = true;
distMethod = 'j';
indelMethod = 's';
nwalignArgs = {};
optargs = {};
gapopen = 8;
gapextend = 8;
distanceFunction = [];
%setGapExtend = false;
%scale=1;

%-% check input sequences 
% they can be an array of chars, a vector of string cells or a vector array
% of structures with the field 'Sequence'

if ischar(seqs)
    seqs(seqs=='.')='-';
    if ~any(seqs(:)==' ')  %no padded sequences (what happens after strvcat of sequences)
        sameLength = true;
        seqs=mat2cell(seqs,ones(1,size(seqs,1)),size(seqs,2));
    else
        sameLength = false;
        pairWiseAlignment = true; % need to align, force pairwise alignment
        seqs = mat2cell(seqs,ones(1,size(seqs,1)),size(seqs,2));
        seqs = regexprep(seqs,'[ .-]',''); % remove padding and gap chars
    end
elseif iscell(seqs) || isfield(seqs,'Sequence')
    if isfield(seqs,'Sequence') % if struct put them in a cell
        seqs = {seqs(:).Sequence};
    end
    seqs = seqs(:);
    L = zeros(length(seqs),1);   % to keep track of the lengths
    seqs = strrep(seqs,' ','');  % padding spaces are not considered gap chars
    seqs = strrep(seqs,'.','-'); % gap chars can only be '-'
    
    % check if sequences have same length
    for i = 1:length(seqs)
        L(i) = length(seqs{i});
    end
    if all(~diff(L)) 
        sameLength = true;
    else
        sameLength = false;
        pairWiseAlignment = true; % need to align, force pairwise alignment
        seqs = strrep(seqs,'-','');  % remove potential gap chars
     end
else
    error(message('bioinfo:seqpdist:IncorrectInputType'))
end

%-% identify input arguments

if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:seqpdist:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'method','indels','pairwisealignment','alphabet',...
              'scoringmatrix','scale','gapopen','extendgap',...
              'squareform','optargs','showdots','useparallel'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs, length(pname)));
        if isempty(k)
            error(message('bioinfo:seqpdist:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:seqpdist:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % select method to compute distances
                    if isa(pval,'function_handle')
                        distMethod = 'u';
                        distanceFunction = pval;
                    else                      
                        distMethods = {'p-distance','jukes-cantor',...
                          'alignment-score','tajima-nei','kimura','tamura',...
                          'hasegawa','nei-tamura','poisson','gamma'};
                        distMethod = strmatch(lower(pval),distMethods); 
                        if isempty(distMethod) 
                            error(message('bioinfo:seqpdist:NotValidMethod'))
                        elseif numel(distMethod)>1
                            error(message('bioinfo:seqpdist:AmbiguousMethod'))
                        else
                            codes = 'pjatkmhnog';
                            distMethod = codes(distMethod);
                        end
                    end
                case 2 % select method to handle gaps
                    indelMethods = {'score','pairwise-delete',...
                                    'complete-delete'};
                    indelMethod = strmatch(lower(pval),indelMethods); 
                    if isempty(indelMethod) 
                        error(message('bioinfo:seqpdist:NotValidIndelMethod'))
                    end
                    indelMethod = indelMethods{indelMethod}(1); % 's' 'p' or 'c'
                case 3 % force pairwise alignment
                    pairWiseAlignment = (pval == true);
                    if ~pairWiseAlignment && ~sameLength
                        error(message('bioinfo:seqpdist:PairwiseAlignmentNeeded'))
                    end
                case 4 % alphabet
                    if strcmpi(pval,'nt')
                        isAmino = false;
                        nwalignArgs = {nwalignArgs{:},'alphabet','NT'}; %#ok<*CCAT>
                    end
                case 5 % scoring matrix
                    if isnumeric(pval)
                        ScoringMatrix = pval;
                    else
                        if ischar(pval)
                            pval = lower(pval);
                        end
                        try
                            [ScoringMatrix,ScoringMatrixInfo] = feval(pval); %#ok
                        catch allExceptions
                            error(message('bioinfo:seqpdist:InvalidScoringMatrix'));
                        end
                    end
                case 6 % scale
                    scale = pval;
                    nwalignArgs = {nwalignArgs{:},'scale',scale};
                case 7 % gap open penalty
                    gapopen = pval;
                    nwalignArgs = {nwalignArgs{:},'gapopen',gapopen};
                case 8 % gap extend penalty
                    gapextend = pval;
                    setGapExtend = true; %#ok 
                    nwalignArgs = {nwalignArgs{:},'extend',gapextend};
                case 9 % square form output
                    squareFormOut = (pval == true);
                case 10 % optional arguments for the distance function
                    if iscell(pval)
                        optargs = pval;
                    else
                        optargs = {pval};
                    end
                case 11 % show dots flag
                    showDots = (pval==true);                
                case 12 % Use parallel
                    useParallel = bioinfoprivate.opttf(pval);
                    if isempty(useParallel)
                        error(message('bioinfo:seqpdist:useParallelNotLogical', upper( char( okargs( k ) ) )));
                    end
            end
        end
    end
end

% If no pairwisealignment is required, make sure it runs in serial mode
useParallel = useParallel && pairWiseAlignment;

% setting the default scoring matrix 
if ~exist('ScoringMatrix','var')
    if isAmino
        [ScoringMatrix,ScoringMatrixInfo] = blosum50; %#ok
    else
        [ScoringMatrix,ScoringMatrixInfo] = nuc44; %#ok
    end
end

% check that the distance method is appropriate for the sequence type
if isAmino 
    if any(distMethod == 'tkmhn')  
        error(message('bioinfo:seqpdist:NotValidDistanceMethodAA'))
    end
elseif any(distMethod == 'og') 
    error(message('bioinfo:seqpdist:NotValidDistanceMethodNT'))
end

% check that some methods cannot score gaps
if (indelMethod == 's') && any(distMethod == 'tkmhnog')
    error(message('bioinfo:seqpdist:NotGapsIndelMethod'))
end

M=length(seqs);           % number of sequences

% compute the overall nucleotide ratio (when needed)
if any(distMethod == 'thn')
    if numel(optargs);
        g = optargs{1};
        g=g(:);
    else
        map('ACGT- ') = 1:6;
        temp = upper(char(seqs));
        g = accumarray(map(temp(:))',1);
        g = g(1:4) / sum(g(1:4));
    end
else
    g = [.25 .25 .25 .25]; % Any value so parfor can execute
end

% compute the GC ratio (when needed)
if distMethod == 'm'
    if numel(optargs);
        gc = optargs{1};
    else
        map('ACGT- ') = 1:6;
        temp = char(seqs);
        gc = accumarray(map(temp(:))',1);
        gc = sum(gc(2:3)) / sum(gc(1:4));
    end
else
    gc = [.5 .5]; % Any value so parfor can execute
end

% Gamma parameter
if distMethod == 'g'
    if numel(optargs);
        a = optargs{1};
    else
        a = 2;
    end
else
    a = 2; % Any value so parfor can execute
end

%-% algorithm starts here

% delete entire columns with gaps if required
if indelMethod == 'c'
    if sameLength 
        h = all(cell2mat(seqs)~='-');
        if sum(h)
            for i = 1:M
                seqs{i} = seqs{i}(h);
            end
        else
            error(message('bioinfo:seqpdist:NoColumnsWithoutGaps'))
        end
    else
        error(message('bioinfo:seqpdist:DeleteColumnOnUnalignedSeqs'))
    end
end 
% delete alignment information if we realign
if pairWiseAlignment && sameLength
    seqs = strrep(seqs,'-','');
end



if useParallel  % main loop implemented with PARFOR
        
    D = cell(1,M-1);
    parfor i1 = 1:M-1
        [vsc,val] = nwalign(seqs(i1),seqs(i1+1:M),nwalignArgs{:},'scoringmatrix',ScoringMatrix);
        if showDots 
            fprintf('.'); 
        end
        j = 1;
        d = zeros(1,M-i1);
        for i2 = i1+1:M
            sc = vsc(j);
            s1 = val{j}(1,:);
            s2 = val{j}(3,:);
        
            % h will contain the locs that need to be scored in the pairwise
            % alignment
            if indelMethod == 's'  % score gaps
                h = (s1~='-') | (s2~='-'); % only avoid sites with double gaps
            else                   % pairwise deleted gaps
                h = (s1~='-') & (s2~='-'); % avoid any site with gaps
            end    
            
            switch distMethod
                case 'p' % p-distance
                    d(j) = pdistance(s1,s2,h);
                case 'j' % Jukes-Cantor
                    d(j) = jukescantor(s1,s2,h,isAmino);
                case 'a' % alignment-score
                    if indelMethod == 's' 
                        d(j) = alignmentscore(s1,s2,h,isAmino,pairWiseAlignment,...
                                           sc,ScoringMatrix,-gapopen,-gapextend);
                    else                  
                        d(j) = alignmentscore(s1,s2,h,isAmino,false,...
                                           0,ScoringMatrix,0,0);
                    end
                case 't' % Tajima-Nei
                    d(j) = tajimanei(s1,s2,h,g);
                case 'k' % Kimura
                    d(j) = kimura(s1,s2,h);
                case 'm' % Tamura
                    d(j) = tamura(s1,s2,h,gc);
                case 'h' % Hasegawa
                    d(j) = hasegawa(s1,s2,h,g);
                case 'n' % Tamura-Nei
                    d(j) = tamuranei(s1,s2,h,g);
                case 'o' % Poisson
                    d(j) = poissondistance(s1,s2,h);
                case 'g' % Gamma
                    d(j) = gammadistance(s1,s2,h,a);
                case 'u' % User defined
                    try
                        d(j) = distanceFunction(s1(h),s2(h),optargs{:});
                    catch theException
                        error(message('bioinfo:seqpdist:ErrorInUserDistanceFunction', func2str( distanceFunction ), theException.message));
                    end
            end %switch    
            j = j + 1;
        end % for i2 = 
        D{i1} = d;
    end % for i1 = 
    D = cell2mat(D); 
    
else % main loop implemented in serial mode
    
    D = zeros(1,(M*(M-1)/2)); % pairwise distances (allocating space)
    id = 1;                   % index to D
    
    for i1 = 1:M-1
        if showDots
            fprintf('.');
        end
        if pairWiseAlignment
            [vsc,val] = nwalign(seqs(i1),seqs(i1+1:M),nwalignArgs{:},'scoringmatrix',ScoringMatrix);
        end
        for i2 = i1+1:M
            if pairWiseAlignment
                sc = vsc(i2-i1);
                s1 = val{i2-i1}(1,:);
                s2 = val{i2-i1}(3,:);
            else
                sc = 0;
                s1 = seqs{i1};
                s2 = seqs{i2};
            end
            
            % h will contain the locs that need to be scored in the pairwise
            % alignment
            if indelMethod == 's'  % score gaps
                h = (s1~='-') | (s2~='-'); % only avoid sites with double gaps
            else                   % pairwise deleted gaps
                h = (s1~='-') & (s2~='-'); % avoid any site with gaps
            end
            
            switch distMethod
                case 'p' % p-distance
                    D(id) = pdistance(s1,s2,h);
                case 'j' % Jukes-Cantor
                    D(id) = jukescantor(s1,s2,h,isAmino);
                case 'a' % alignment-score
                    if indelMethod == 's'
                        D(id) = alignmentscore(s1,s2,h,isAmino,pairWiseAlignment,...
                            sc,ScoringMatrix,-gapopen,-gapextend);
                    else
                        D(id) = alignmentscore(s1,s2,h,isAmino,false,...
                            0,ScoringMatrix,0,0);
                    end
                case 't' % Tajima-Nei
                    D(id) = tajimanei(s1,s2,h,g);
                case 'k' % Kimura
                    D(id) = kimura(s1,s2,h);
                case 'm' % Tamura
                    D(id) = tamura(s1,s2,h,gc);
                case 'h' % Hasegawa
                    D(id) = hasegawa(s1,s2,h,g);
                case 'n' % Tamura-Nei
                    D(id) = tamuranei(s1,s2,h,g);
                case 'o' % Poisson
                    D(id) = poissondistance(s1,s2,h);
                case 'g' % Gamma
                    D(id) = gammadistance(s1,s2,h,a);
                case 'u' % User defined
                    try
                        D(id) = distanceFunction(s1(h),s2(h),optargs{:});
                    catch theException
                        error(message('bioinfo:seqpdist:ErrorInUserDistanceFunction', func2str( distanceFunction ), theException.message));
                    end
            end %switch
            id = id + 1;
        end
    end
    
end % isDistributedComputing || useParallel || serial mode    
    
if showDots 
    fprintf('\n'); 
end  

if squareFormOut
   D = squareform(D); % use stats function
end

%-% p-distance
function d = pdistance(s1,s2,h)
% p-distance
d = sum( (s1~=s2) & h ) / sum(h);

%-% Jukes-Cantor distance
function d = jukescantor(s1,s2,h,isAmino)
% Jukes-Cantor distance
%
% Reference:  Jukes, T. H. and C. R. Cantor. (1969) Mammalian Protein
%            Metabolism Acad. Press
f = sum( (s1~=s2) & h ) / sum(h);
if isAmino       
   d = -19/20 * log(max(eps,1-20*f/19)); % max prevents negative log
else
   d = -3/4 * log(max(eps,1-4*f/3));
end

%-% alignment-score distance
function d = alignmentscore(s1,s2,h,isAmino,pwa,sc,SM,go,ge)
% Alignment-score distance
sms = uint16(size(SM,1));
if isAmino                    
    sint1 = uint16(aa2int(s1(h)));
    sint2 = uint16(aa2int(s2(h)));
    gaps1 = sint1==25; gaps2 = sint2==25;
else
    sint1 = uint16(nt2int(s1(h)));
    sint2 = uint16(nt2int(s2(h)));
    gaps1 = sint1==16; gaps2 = sint2==16;
end                  
if ~pwa % compute the score of the alignment if seqs have not been aligned
    gaps      = [gaps1;gaps2];
    firstGaps = any(diff([[0;0] gaps],[],2)==1);
    otherGaps = any(gaps) & ~firstGaps;
    matchLocs = ~any(gaps);
    sc = sum(SM((sint1(matchLocs)-1) * sms + sint2(matchLocs))) + ...
         sum(firstGaps)*go + sum(otherGaps)*ge;
end
sc1 = sum(SM((sint1(~gaps1)-1) * sms + sint1(~gaps1)));
sc2 = sum(SM((sint2(~gaps2)-1) * sms + sint2(~gaps2)));
d =  max(0, sc1 - sc) * max(0, sc2 - sc) /sc1/sc2;
     % max(0,x) prevents the rare case when sc > scx
                
%-% Tajima-Nei distance
function d = tajimanei(s1,s2,h,g)
% Tajima-Nei distance
% 
% Let X be the nucleotide substitution count (4x4 matrix), and g the
% nucleotide frequencies over all the sequences:
% 
% fX = (X-diag(diag(X)))/sum(X(:)); %frequency of substitutions
% p  = 1-sum(diag(X))/sum(X(:));    % proportion of substitutions
% gi = (1./g);
% c  = gi' * (fX.*fX) * gi /2 ;
% b  = (1-g'*g+(p^2)/c) /2;
% d  = - b * log( 1 - p/b);
%
% g is computed from the input sequences or may be given by the user by
% setting 'OPTARGS' to [gA gC gG gT].
%
% Reference: Tajima and Nei (1984) Molecular Biology and Evolution

X  = accumarray(double([nt2int(s1(h))',nt2int(s2(h))']),1,[4 4]);
fX = (X-diag(diag(X)))/sum(X(:)); %frequency of substitutions
p  = 1-sum(diag(X))/sum(X(:));    % proportion of substitutions
gi = (1./g);
c  = gi' * (fX.*fX) * gi /2 ;
b  = (1-g'*g+(p^2)/c) /2;
d  = - b * log( 1 - p/b);

%-% Kimura distance
function d = kimura(s1,s2,h)
% Kimura distance
%
% Let P and Q be the transition and transversional difference proportions:
% d = - log(1 -2*P -Q)/2 - log(1 -2*Q)/4;
%
% Reference: Kimura, M. (1980)  Journal of Molecular Evolution
X  = accumarray(double([nt2int(s1(h))',nt2int(s2(h))']),1,[4 4]);
numPairs = sum(X(:));
numTransitions = sum(X([3 8 9 14]));
numTransversions = numPairs - sum(diag(X)) - numTransitions;
P = numTransitions / numPairs;
Q = numTransversions / numPairs;
d = - log(1 -2*P -Q)/2 - log(1 -2*Q)/4;

%-% Tamura distance
function d = tamura(s1,s2,h,gc)
% Tamura distance
%
% Let P and Q be the transition and transversional difference proportions
% and gc the average GC content over all the sequences:
% k = 2 * gc * (1-gc);
% d = - k * log(1-P/k-Q) - (1-k) * log(1-2*Q)/2;
%
% gc is computed from the input sequences or may be given by the user by
% setting 'OPTARGS' to the GC proportion.
%
% Reference: Tamura, K. (1992) Molecular Biology and Evolution
X  = accumarray(double([nt2int(s1(h))',nt2int(s2(h))']),1,[4 4]);
numPairs = sum(X(:));
numTransitions = sum(X([3 8 9 14]));
numTransversions = numPairs - sum(diag(X)) - numTransitions;
P = numTransitions / numPairs;
Q = numTransversions / numPairs;
k = 2 * gc * (1-gc);
d = - k * log(1-P/k-Q) - (1-k) * log(1-2*Q)/2;

%-% Hasegawa distance
function d = hasegawa(s1,s2,h,g)
% Hasegawa distance
%
% Let P and Q be the transition and transversional difference proportions
% and g the nucleotide frequencies over all the sequences:
%
% gr = g(1)+g(3);  % purines frequency
% gy = g(2)+g(4);  % pyrimidines frequency
% k1 = 2*g(1)*g(3) / gr;
% k2 = 2*g(2)*g(4) / gy;
% d = - k1 * log(1-P/2/k1-Q/gr) - k2 * log(1-P/2/k2-Q/gy) ...
%     - (2*gr*gy-gy*k1-gr*k2) * log(1-Q/gr/gy/2);
% 
% g is computed from the input sequences or may be given by the user by
% setting 'OPTARGS' to [gA gC gG gT].
%
% Reference: Tamura and  Nei (1993) Molecular Biology and Evolution
X  = accumarray(double([nt2int(s1(h))',nt2int(s2(h))']),1,[4 4]);
numPairs = sum(X(:));
numAGTransitions = sum(X([3 9]));
numCTTransitions = sum(X([8 14]));
numTransitions = numCTTransitions + numAGTransitions;
numTransversions = numPairs - sum(diag(X)) - numTransitions;
P = numTransitions / numPairs;
Q = numTransversions / numPairs;
gr = g(1)+g(3); % purines frequency
gy = g(2)+g(4); % pyrimidines frequency
k1 = 2*g(1)*g(3) / gr;
k2 = 2*g(2)*g(4) / gy;
d = - k1 * log(1-P/2/k1-Q/gr) - k2 * log(1-P/2/k2-Q/gy) ...
    - (2*gr*gy-gy*k1-gr*k2) * log(1-Q/gr/gy/2);

%-% Tamura-Nei distance
function d = tamuranei(s1,s2,h,g)
% Tamura-Nei distance
%
% Let P1, P2, and Q be the purines transition, pyrimidines transition and
% the transversional difference proportions, and g the nucleotide
% frequencies over all the sequences:
%
% gr = g(1)+g(3);  % purines frequency
% gy = g(2)+g(4);  % pyrimidines frequency
% k1 = 2*g(1)*g(3) / gr;
% k2 = 2*g(2)*g(4) / gy;
% d = - k1 * log(1-P1/k1-Q/gr) - k2 * log(1-P2/k2-Q/gy) ...
%     - (2*gr*gy-gy*k1-gr*k2) * log(1-Q/gr/gy/2);
% 
% g is computed from the input sequences or may be given by the user by
% setting 'OPTARGS' to [gA gC gG gT].
%
% Reference: Tamura and  Nei (1993) Molecular Biology and Evolution
X  = accumarray(double([nt2int(s1(h))',nt2int(s2(h))']),1,[4 4]);
numPairs = sum(X(:));
numAGTransitions = sum(X([3 9]));
numCTTransitions = sum(X([8 14]));
numTransitions = numCTTransitions + numAGTransitions;
numTransversions = numPairs - sum(diag(X)) - numTransitions;
P1 = numAGTransitions / numPairs;
P2 = numCTTransitions / numPairs;
% P = numTransitions / numPairs;
Q = numTransversions / numPairs;
gr = g(1)+g(3); % purines frequency
gy = g(2)+g(4); % pyrimidines frequency
k1 = 2*g(1)*g(3) / gr;
k2 = 2*g(2)*g(4) / gy;
d = - k1 * log(1-P1/k1-Q/gr) - k2 * log(1-P2/k2-Q/gy) ...
    - (2*gr*gy-gy*k1-gr*k2) * log(1-Q/gr/gy/2);

%-% Poisson distance
function d = poissondistance(s1,s2,h)
% Poisson distance
% 
% Let p be the proportion of substitutions:
% d = -log(1-p);
f = sum( (s1~=s2) & h ) / sum(h);
d = -log(1-f);
 
%-% Gamma distance
function d = gammadistance(s1,s2,h,a)
% Gamma distance
% Let p be the proportion of substitutions and 'a' the Gamma parameter:
% d = a * ( (1-p)^(-1/a) - 1);
%
% default of 'a' is 2 unless another value is given by the user by setting
% 'OPTARGS'
f = sum( (s1~=s2) & h ) / sum(h);
d = a * ( (1-f)^(-1/a) - 1);

% OBSOLETED
% %-%-----------------------------------------------------------------------
% function C = splittasks(N,W)
% % Helper function to calculate the best task partition that balances the
% % distributed computation of the pairwise alignments. 
% %
% % N - number of sequences
% % W - number of available workers
% % C - Columns of the half pairwise matrix assigned to every task
% 
% if W==0 % if no workers available return a default partition
%     C = ones(N-1,1);
%     return
% end
% 
% % M - number of pairwise comparisons
% M = (N-1)*N/2;
% % I - ideal partitions
% Ka = 0.9;
% I = M/W*[Ka*(1:W)  W*Ka+(1-Ka)*(1:W)];
% % P - Possible partitions
% P = cumsum(N-1:-1:1);
% 
% % balancing the load for every task
% C = zeros(2*W,1);
% c = 1;
% for i = 1:2*W
%    while P(c)<I(i)
%        c = c+1;
%    end
%    C(i) = c;
% end
% C = diff([0;C]);
% 
% % remove tasks with no work
% C(C==0)=[];

