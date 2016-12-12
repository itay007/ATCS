function [P,S] = seqprofile(seqs,varargin)
%SEQPROFILE computes the sequence profile of a multiple alignment
%
%  P = SEQPROFILE(SEQS) returns a matrix P of size [20 (or 4) x seq Length]
%  with the frequency of amino acids (or nucleotides) for every column in
%  the multiple alignment. The order of the rows is given by AA2INT (or
%  NT2INT). SEQS can be a char array, a cell array of strings, or an array
%  of structures with the field 'Sequence'.
%
%  [P,S] = SEQPROFILE(SEQS) returns a unique symbol list. Every symbol in S
%  corresponds to a row in P.
%
%  SEQPROFILE(...,'ALPHABET',A) specifies whether the sequences are amino
%  acids ('AA') or nucleotides ('NT'). When A is 'none', SEQPROFILE creates
%  a symbol list based on the observed symbols. Every character can be a
%  symbol except '-' and '.', which are reserved for gaps. The default
%  alphabet is 'AA'.
%
%  SEQPROFILE(...,'COUNTS',true) returns the counts instead of the
%  frequency.
%
%  SEQPROFILE(...,'GAPS',G) appends a row to the bottom of P with the count
%  for gaps. Options are 'all' (counts all gaps), 'noflanks' (counts all
%  gaps except those at the flanks of every sequence), and 'none'. The
%  default is 'none'.
%
%  SEQPROFILE(...,'AMBIGUOUS','COUNT') counts ambiguous symbols. For
%  example, in calculating amino acid frequencies, an X in a sequence adds
%  0.05 (1/20) to every row, while a B adds 0.50 (1/2) to the D and N rows.
%  The default is 'IGNORE', i.e., ambiguous symbols are ignored.
%
%  SEQPROFILE(...,'LIMITS',L) sets the start and end positions for the
%  profile relative to the indices of the multiple alignment. L should be a
%  1-by-2 vector and defaults to [1,inf].
%
%  Example:
%
%      seqs = fastaread('pf00002.fa');
%      [P,S] = seqprofile(seqs,'limits',[50 60],'gaps','all')
%
%  See also FASTAREAD, MULTIALIGNREAD, PROFALIGN, SEQCONSENSUS, SEQDISP,
%   SEQLOGO.

% Copyright 2003-2012 The MathWorks, Inc.


% set defaults
alpha = 'aa';
gaps = 'none';
counts = false;
ambiguous = false;
limits = [0 inf];

nvarargin = numel(varargin);
if nvarargin
    if rem(nvarargin,2)
        error(message('bioinfo:seqprofile:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'alphabet','counts','gaps','ambiguous','limits'};
    for j=1:2:nvarargin
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname))); 
        if isempty(k)
            error(message('bioinfo:seqprofile:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:seqprofile:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % alphabet
                    alphaOptions = {'aa','nt','none'};
                    alpha = strmatch(lower(pval),alphaOptions); 
                    if isempty(alpha)
                        error(message('bioinfo:seqprofile:NotValidAlphabet'))
                    end
                    alpha = alphaOptions{alpha};
                case 2 % counts
                    counts = bioinfoprivate.opttf(pval);
                case 3 % gaps
                    gapsOptions = {'all','noflanks','none'};
                    gaps = strmatch(lower(pval),gapsOptions); 
                    if isempty(gaps)
                        error(message('bioinfo:seqprofile:NotValidGapsOptions'))
                    end
                    gaps = gapsOptions{gaps};
                case 4 % ambiguous
                    if ischar(pval)
                        ambiguous = ~isempty(strmatch(lower(pval),'count')); 
                    else
                        error(message('bioinfo:seqprofile:NotValidAmbiguousOptions'))
                    end
                case 5 % limits
                    limits = round(pval(:)');
                    if numel(limits)~=2 || diff(limits)<=0
                        error(message('bioinfo:seqprofile:badLimits'))
                    end
            end
        end
    end
end

% check input sequences
% they can be an array of chars, a vector of string cells or an array of
% structures with the field 'Sequence'

if iscell(seqs) || isfield(seqs,'Sequence')
    if isfield(seqs,'Sequence') % if struct put them in a cell
        seqs = {seqs(:).Sequence};
    end
    seqs = seqs(:);
    seqs = strrep(seqs,' ',''); % padding spaces are not considered 'align' chars
    seqs = strvcat(seqs); %#ok<VCAT> % now seqs must be a char array
end

if ~ischar(seqs)
    error(message('bioinfo:seqprofile:IncorrectInputType'))
end

% ... at this point seqs can only be a char array.

seqs = upper(seqs); % work in uppercase
seqs(seqs=='.')='-'; % align symbols can only be '-'
if any(seqs(:)==' ')  %aligned ?
    error(message('bioinfo:seqprofile:IncorrectMultipleAlignment'))
end

% Find the gaps and eliminate gapflanks is required
[numSeqs,numCols] = size(seqs);
if isequal(gaps,'noflanks')
    for i = 1:numSeqs
        seqs(i,1:find(seqs(i,:)~='-',1,'first')-1)=' ';
        seqs(i,find(seqs(i,:)~='-',1,'last')+1:end)=' ';
    end
end

% Trim arrays if necessary
limits(1) = max(1,limits(1));
limits(2) = min(numCols,limits(2));
if any(limits~=[1,numCols])
    seqs = seqs(:,limits(1):limits(2));
    numCols = diff(limits)+1;
end

switch alpha
    case 'aa'
        C = zeros(26,numCols);
        % reshape necessary as XX2int flips column vectors
        intSeqs = reshape(aa2int(seqs),size(seqs));
        intSeqs(~intSeqs) = 26; %aa2int marks unknowns as 0
        for i = 1:numCols
            C(:,i)=accumarray(double(intSeqs(:,i)),1,[26,1]);
        end
        if isequal(gaps,'none')
            P = C(1:20,:);
            S = int2aa(1:20);
        else
            P = C([1:20 25],:);
            S = int2aa([1:20 25]);
        end
        if ambiguous
            P(1:20,:)  = P(1:20,:) +repmat(C(23,:)/20,20,1);  % X -> any
            P([6 7],:) = P([6 7],:)+repmat(C(22,:)/2,2,1);    % Z -> EQ
            P([3 4],:) = P([3 4],:)+repmat(C(21,:)/2,2,1);    % B -> ND
        end
    case 'nt'
        C = zeros(17,numCols);
        seqs(seqs==' ')='*'; % spaces cannot be counted by nt2int
        % reshape necessary as XX2int flips column vectors
        intSeqs = reshape(nt2int(seqs,'unknown',17),size(seqs));
        for i = 1:numCols
            C(:,i)=accumarray(double(intSeqs(:,i)),1,[17,1]);
        end
        if isequal(gaps,'none')
            P = C(1:4,:);
            S = int2nt(1:4);
        else
            P = C([1:4 16],:);
            S = int2nt([1:4 16]);
        end
        if ambiguous
            P(1:4,:)  = P(1:4,:) +repmat(C(15,:)/4,4,1);       % N -> ACGT
            P([1 2 3],:) = P([1 2 3],:)+repmat(C(14,:)/3,3,1); % V -> ACG
            P([1 2 4],:) = P([1 2 4],:)+repmat(C(13,:)/3,3,1); % H -> ACT
            P([1 3 4],:) = P([1 3 4],:)+repmat(C(12,:)/3,3,1); % D -> AGT
            P([2 3 4],:) = P([2 3 4],:)+repmat(C(11,:)/3,3,1); % B -> CGT
            P([1 4],:) = P([1 4],:)+repmat(C(10,:)/2,2,1);     % W -> AT
            P([2 3],:) = P([2 3],:)+repmat(C(9,:)/2,2,1);      % S -> CG
            P([1 2],:) = P([1 2],:)+repmat(C(8,:)/2,2,1);      % M -> AC
            P([3 4],:) = P([3 4],:)+repmat(C(7,:)/2,2,1);      % K -> GT
            P([2 4],:) = P([2 4],:)+repmat(C(6,:)/2,2,1);      % Y -> CT
            P([1 3],:) = P([1 3],:)+repmat(C(5,:)/2,2,1);      % R -> AG
        end
    otherwise % makes its own alphabet up (for seqlogo)
        C = zeros(64,numCols);
        intSeqs = min(64,max(0,seqs-32)); % symbols from 0 to 64
        intSeqs(~intSeqs) = 64; % will be the unknown (not counted)
        for i = 1:numCols
            C(:,i)=accumarray(double(intSeqs(:,i)),1,[64,1]);
        end
        uniqueRows = setdiff(unique(intSeqs),[13 64]); %no gaps and unknowns
        if ~isequal(gaps,'none')
            uniqueRows = [uniqueRows 13];
        end
        P = C(uniqueRows,:);
        S = char(uniqueRows+32);
end

if ~counts  %compute frequency
    Q = sum(P);
    Q(Q==0) = 1;
    P = P./repmat(Q,size(P,1),1);
end
