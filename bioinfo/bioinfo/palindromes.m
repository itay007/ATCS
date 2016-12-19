function [pos,plen,sequences,gps] = palindromes(seq,varargin)
%PALINDROMES finds palindromes in a sequence.
%
%   [POS, PLENGTH] = PALINDROMES(SEQ) finds any palindromes of length
%   greater than or equal to 6 in sequence SEQ and returns POS, the
%   starting indices, and PLENGTH, the lengths of the palindromes.
%
%   [POS,PLENGTH,PAL] = PALINDROMES(SEQ) also returns a cell array, PAL, of
%   the palindromes.
%
%   PALINDROMES(...,'LENGTH',LEN) finds all palindromes of length greater
%   than or equal to LEN. The default value is 6.
%
%   PALINDROMES(...,'GAP',GAPLEN) allows gaps of up to length GAPLEN in the
%   middle of the palindromes. When GAPLEN is greater than 1, gaps are
%   shown in lower case and the palindrome is shown in upper case.
%
%   PALINDROMES(...,'COMPLEMENT',true) finds complementary palindromes,
%   that is, where the elements match their complementary pairs A-T(or U)
%   and C-G instead of an exact nucleotide match.
%
%   Examples:
%
%       [p,l,s] = palindromes('GCTAGTAACGTATATATAAT')
%
%       [pc,lc,sc] = palindromes('GCTAGTAACGTATATATAAT','complement',true)
%
%   See also REGEXP, SEQRCOMPLEMENT, SEQSHOWWORDS, STRFIND.

%   Copyright 2002-2012 The MathWorks, Inc.


% If the input is a structure then extract the Sequence data.
if isstruct(seq)
    seq = bioinfoprivate.seqfromstruct(seq);
end

minPal = 6;
useComplement = false;
gapLen = 1;
if nargin > 1
    if rem(nargin,2)== 0
        error(message('bioinfo:palindromes:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'length','complement','gap'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:palindromes:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:palindromes:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1   % length
                    minPal = pval;
                    if minPal < 2
                        warning(message('bioinfo:palindromes:LengthTooShort'));
                        minPal = 2;
                    end
                case 2 %  complement
                    useComplement = bioinfoprivate.opttf(pval);
                    if isempty(useComplement)
                        error(message('bioinfo:palindromes:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 3   % gap
                    gapLen = max(1,pval);
            end
        end
    end
end

% Save original sequence for output purposes
origseq = seq;
% work in lower case
seq = lower(seq);

% make a copy of sequence -- use complement if requested
cseq = seq;
if useComplement
    cseq = seqcomplement(seq);
    if bioinfoprivate.isrna(seq) && ~bioinfoprivate.isrna(cseq)
        cseq = dna2rna(cseq);
    end
end

% make the hits vector to store lengths of palindromes
seqLength = numel(seq);
hits = zeros(seqLength,1);
gapHits = hits;
gapStep = floor(gapLen/2);

% for each position in the sequence
for outer = floor(minPal/2):(seqLength-floor(minPal/2))
    % look for palindromes centred at the position with gap at most 2*gap+gapOdd
    for gap = 0:gapStep
        for gapOdd = 0:1
            len = 2*gap+gapOdd;
            if len <= gapLen
                start = 0;
                % figure out how far we can look to the left or to the
                % right
                for  inner = 1:min(outer-gapOdd-gap,seqLength-outer-gap)
                    % if we have a match, add two to the count and note
                    % starting point
                    if (seq(outer-(inner-1+gapOdd)-gap) == cseq(outer+inner+gap))
                        len = len + 2;
                        start = outer-(inner-1+gapOdd)-gap;
                    else
                        % if not, jump out of inner loop
                        break;
                    end
                end
                % if we had a match start > 0
                if (start > 0) && (len > hits(start))
                    % if len > hits(start) then we have found a new best
                    % length palindrome for that starting point.
                    hits(start) = len;
                    gapHits(start) = 2*gap+gapOdd;
                    if gapHits(start) == 1
                        gapHits(start) = 0;
                    end
                end
            end
        end
    end
end


pos = find(hits>=minPal);
gps = gapHits(pos);
% now extract the palindrome for the hits > length
if nargout > 1
    % get the palindrome length
    plen = hits(pos);
    if nargout > 2 % extract the actual sequences
        numPalindromes = length(pos);
        sequences = cell(numPalindromes,1);
        for count = 1:numPalindromes
            sequences{count} = origseq(pos(count):pos(count)+hits(pos(count))-1);
            % when we are working with gaps, show the gap in lower case and
            % the palindrome in upper case.
            if gapLen > 1
                thePal = upper(sequences{count});
                if gps(count) > 1
                    startGap = (floor(length(thePal)/2)-floor(gps(count)/2))+1;
                    stopGap = startGap +gps(count)-1 ;
                    thePal(startGap:stopGap) = lower(thePal(startGap:stopGap));
                end
                sequences{count} = thePal;
            end
        end
    end
end

