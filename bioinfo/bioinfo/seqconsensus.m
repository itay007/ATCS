function [C,S] = seqconsensus(P,varargin)
%SEQCONSENSUS computes the consensus sequence for a set of sequences
%
%  C = SEQCONSENSUS(SEQS) returns a string with the consensus sequence for
%  SEQS, an aligned set of sequences. SEQS can be a char array, a cell
%  array of strings, or an array of structures with the field 'Sequence'.
%
%  [C,S] = SEQCONSENSUS(SEQS) returns the conservation score of the
%  consensus sequence. Scores are computed with BLOSUM50 for AA or NUC44
%  for NT. Scores are the average Euclidean distance between the scored
%  symbol and the M-dimensional consensus value. M is the size of the
%  alphabet. The consensus value is the profile weighted by the scoring
%  matrix.
%
%  C = SEQCONSENSUS(P) P is a sequence profile such as the one returned by
%  SEQPROFILE. P is a matrix of size [20 (or 4) x seq Length] with the
%  frequency/count of amino acids (or nucleotides) for every position. P
%  can also have 21 (or 5) rows if gaps are included in the consensus. 
%
%  SEQCONSENSUS(...,'SCORINGMATRIX',SM) specifies the scoring matrix. The
%  default is BLOSUM50 for AA or NUC44 for NT. SM can also be a 21-by-21,
%  5-by-5, 20-by-20, or 4-by-4 numeric array. For the last two cases, gap
%  scores (last row/column) are set to mean(diag(SM)) for a gap matching
%  another gap, and to mean(nodiag(SM)) for a gap matching another symbol.
%
%  The following input parameters are analogous to SEQPROFILE; however,
%  alphabet is restricted to 'AA' or 'NT' only. Type help SEQPROFILE for
%  more info:
%
%  SEQCONSENSUS(...,'ALPHABET',A)
%  SEQCONSENSUS(...,'GAPS',G)
%  SEQCONSENSUS(...,'AMBIGUOUS',C)
%  SEQCONSENSUS(...,'LIMITS',L)
%
%  Example:
%
%      seqs = fastaread('pf00002.fa');
%      [C,S] = seqconsensus(seqs,'limits',[50 60],'gaps','all')
%
%  See also FASTAREAD, MULTIALIGNREAD, PROFALIGN, SEQDISP, SEQLOGO,
%   SEQPROFILE.

% Copyright 2003-2008 The MathWorks, Inc.


% set defaults
isAAalpha = true;
alpha = 'aa';
gaps = 'none';
ambiguous = 'ignore';
limits = [0 inf];
numericSMProvided = false;
predefSMProvided  = false;
gapScoresIncluded = false;

nvarargin = numel(varargin);
if nvarargin
    if rem(nvarargin,2)
        error(message('bioinfo:seqconsensus:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'alphabet','gaps','ambiguous','limits','scoringmatrix'};
    for j=1:2:nvarargin
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname))); 
        if isempty(k)
            error(message('bioinfo:seqconsensus:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:seqconsensus:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % alphabet
                    alphaOptions = {'aa','nt'};
                    alpha = strmatch(lower(pval),alphaOptions); 
                    if isempty(alpha)
                        error(message('bioinfo:seqconsensus:NotValidAlphabet'))
                    end
                    alpha = alphaOptions{alpha};
                    isAAalpha = isequal(alpha,'aa');
                case 2 % gaps
                    gapsOptions = {'all','noflanks','none'};
                    gaps = strmatch(lower(pval),gapsOptions); 
                    if isempty(gaps)
                        error(message('bioinfo:seqconsensus:NotValidGapsOptions'))
                    end
                    gaps = gapsOptions{gaps};
                case 3 % ambiguous
                    ambiguousOptions = {'count','ignore'};
                    ambiguous =  strmatch(lower(pval),ambiguousOptions); 
                    if isempty(ambiguous)
                        error(message('bioinfo:seqconsensus:NotValidAmbiguousOptions'))
                    end
                    ambiguous = ambiguousOptions{ambiguous};
                case 4 % limits
                    limits = round(pval(:)');
                    if numel(limits)~=2 || diff(limits)<=0
                        error(message('bioinfo:seqconsensus:badLimits'))
                    end
                case 5 % scoring matrix
                    if isnumeric(pval)
                        ScoringMatrix = pval;
                        numericSMProvided = true;
                    else
                        predefSMProvided = true;
                        if ischar(pval)
                            pval = lower(pval);
                        end
                        try
                            ScoringMatrix = feval(pval);
                        catch allExceptions
                            error(message('bioinfo:seqconsensus:InvalidScoringMatrix'));
                        end
                    end
            end
        end
    end
end

% setting the default scoring matrix
if ~numericSMProvided && ~predefSMProvided
    if isAAalpha
        ScoringMatrix = blosum50;
    else
        ScoringMatrix = nuc44;
    end
end

% check sizes and type of scoring matrix
if numericSMProvided
    if all(size(ScoringMatrix)==21) ||  all(size(ScoringMatrix)==5)
        gapScoresIncluded = true;
    elseif any(size(ScoringMatrix)~=20) &&  any(size(ScoringMatrix)~=4)
        error(message('bioinfo:seqconsensus:InvalidNumericScoringMatrix'));
    end
else
    if isAAalpha
        ScoringMatrix = ScoringMatrix(1:20,1:20);
    else
        ScoringMatrix = ScoringMatrix(1:4,1:4);
    end
end

if isnumeric(P) %P has the profile
    if isAAalpha && all(size(P,1)~=[20 21])
        error(message('bioinfo:seqconsensus:IncorrectProfileAA'))
    elseif  ~isAAalpha && all(size(P,1)~=[4 5])
        error(message('bioinfo:seqconsensus:IncorrectProfileNT'))
    end
else %P has the aligned sequences, call SEQPROF internally to get the profile
    % check input sequences
    % they can be an array of chars, a vector of string cells or a vetor
    % array  of structures with the field 'Sequence'

    if iscell(P) || isfield(P,'Sequence')
        if isfield(P,'Sequence') % if struct put them in a cell
            P = {P(:).Sequence};
        end
        P = P(:);
        P = strrep(P,' ',''); % padding spaces are not considered 'align' chars
        P = char(P); % now seqs must be a char array
    end

    if ~ischar(P)
        error(message('bioinfo:seqconsensus:IncorrectInputType'))
    end
    SEQS = P;

    P = seqprofile(SEQS,'alpha',alpha,'amb',ambiguous,'lim',limits,'gaps',gaps);
end

% append the gap scores to SM when necessary
SMs = size(ScoringMatrix,1);
if ~gapScoresIncluded && any(size(P,1)==[5 21])
    temp = (sum(ScoringMatrix,2)-diag(ScoringMatrix))/(SMs-1);
    ScoringMatrix(SMs+1,SMs+1) = mean(diag(ScoringMatrix));
    ScoringMatrix(1:SMs,SMs+1) = temp;
    ScoringMatrix(SMs+1,1:SMs) = temp';
    SMs=SMs+1;
end

% make sure the profile is as frequency (and not as a count), so for the
% computations we use product vector instead of means.
P = P./repmat(sum(P),size(P,1),1);

% consensus values
X = ScoringMatrix*P;

% find out the consensus symbols
% note: take care of floating point rounding error carried during matrix
%       multiplication
[mX,idx] = max(repmat(max(X),SMs,1) < X + SMs*eps ); 

if isAAalpha
    C = int2aa(idx);
    C(idx==21)='-';
else
    C = int2nt(idx);
    C(idx==5)='-';
end

if nargout>1
    % Euclidean distances to consensus and then dot mult by the profile
    S = zeros(1,size(P,2));
    for j=1:size(P,2)
        S(j) = sqrt(sum((repmat(X(:,j),1,size(P,1))-ScoringMatrix).^2))*P(:,j);
    end
end

