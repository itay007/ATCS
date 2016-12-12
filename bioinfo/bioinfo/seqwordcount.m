function count = seqwordcount(seq,word,varargin)
%SEQWORDCOUNT reports the number of occurrences of a word in a sequence.
%
%   SEQWORDCOUNT(SEQ,WORD) counts the number of times that the word WORD
%   appears in the sequence SEQ.
%
%   If WORD contains nucleotide or amino acid symbols that represent
%   multiple possible symbols, then SEQWORDCOUNT counts all matches. For
%   example, the symbol R represents either G or A (purines). If WORD
%   is 'ART', then SEQWORDCOUNT counts occurrences of both 'AAT' and 'AGT'.
%
%   SEQWORDCOUNT(...,'ALPHABET',A) specifies that SEQ and WORD are
%   amino acids ('AA') or nucleotides ('NT'). The default is NT.
%
%   Example:
%       seqwordcount('GCTAGTAACGTATATATAAT','BART');
%       reports 2 matches  ('TAGT' and 'TAAT')
%
%   SEQWORDCOUNT does not count overlapping patterns multiple times.
%
%   Examples:
%
%       seqwordcount('GCTAGTAACGTATATATAAT','BART')                     
%       % This example reports 2 matches  ('TAGT' and 'TAAT')                            
%
%       % SEQWORDCOUNT does not count overlapping patterns multiple times.
%       seqwordcount('GCTATAACGTATATATAT','TATA')
%       % This reports 3 matches. TATATATA is counted as two distinct
%       % matches, not three overlapping occurrences.
%
%   See also CODONCOUNT, SEQ2REGEXP, SEQMATCH, SEQSHOWORFS, SEQSHOWWORDS,
%   SEQVIEWER, STRFIND.

%   Copyright 2002-2012 The MathWorks, Inc.

% If we are sure that both sequence and word are both DNA/RNA or both AA
% use seq2regexp to deal with extended symbols.

% If the input is a structure then extract the Sequence data.
if isstruct(seq)
    seq = bioinfoprivate.seqfromstruct(seq);
end

alphabet = '';
if (bioinfoprivate.isnt(seq) && bioinfoprivate.isnt(word))
    alphabet = 'nt';
elseif (bioinfoprivate.isaa(seq) && bioinfoprivate.isaa(seq))
    alphabet = 'aa';
end;

if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:seqwordcount:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'alphabet'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname))); 
        if isempty(k)
            error(message('bioinfo:seqwordcount:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:seqwordcount:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1
                    alphabet = pval;
                    if strcmpi(pval,'aa')
                        alphabet = 'aa';
                    end
            end
        end
    end
end

if ~isempty(alphabet)
    word = seq2regexp(word,'ALPHABET',alphabet);
end

count = length(regexpi(seq,word));
