function out=seq2regexp(seq,varargin)
%SEQ2REGEXP converts extended NT or AA symbols into a regular expression.
%
%   SEQ2REGEXP(SEQUENCE) converts extended nucleotide or amino acid symbols
%   in SEQUENCE into regular expression format.
%
%   SEQ2REGEXP(...,'ALPHABET',type) specifies whether the sequence is amino
%   acids ('AA') or nucleotides ('NT'). The default is NT.
%
%   SEQ2REGEXP(...,'AMBIGUOUS',false) removes the ambiguous characters from
%   the output regular expressions. This was the default behavior in older
%   versions of the Bioinformatics Toolbox.
%
%   IUB/IUPAC nucleic acid code conversions:
%
%   A --> A                   M --> [AC] (amino)
%   C --> C                   S --> [GC] (strong)
%   G --> G                   W --> [AT] (weak)
%   T --> T                   B --> [GTC]
%   U --> U                   D --> [GAT]
%   R --> [GA] (purine)       H --> [ACT]
%   Y --> [TC] (pyrimidine)   V --> [GCA]
%   K --> [GT] (keto)         N --> [AGCT] (any)
%
%   Amino acid conversions:
%
%   B --> [DN] 	aspartic acid or asparagine
%   Z --> [EQ]	glutamic acid or glutamine
%   X --> [ARNDCQEGHILKMFPSTWYV]
%
%   Example:
%
%      r = seq2regexp('ACWTMAN')
%
%   See also REGEXP, REGEXPI, RESTRICT, SEQWORDCOUNT.

%   Copyright 2002-2012 The MathWorks, Inc.

% If the input is a structure then extract the Sequence data.
if isstruct(seq)
    seq = bioinfoprivate.seqfromstruct(seq);
end

isAminoAcid = false;
isNucleotide = false;
useAmbiguous = true;
if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:seq2regexp:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'alphabet','ambiguous'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:seq2regexp:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:seq2regexp:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1
                    if strcmpi(pval,'aa')
                        isAminoAcid = true;
                    elseif  strcmpi(pval,'nt')
                        isNucleotide = true;
                    end
                case 2  % old behaviour ACGT Only
                    useAmbiguous = bioinfoprivate.opttf(pval);
                    if isempty(useAmbiguous)
                        error(message('bioinfo:seq2regexp:InputOptionNotLogical', okargs{ k }));
                    end
            end
        end
    end
end


%out=upper(seq);
[m,p] = regexp(seq, '(?<!\\)\w', 'match'); % find non-metacharachers
seq(p) = upper([m{:}]); % change non-metacharacters to upper case

% Some databases use parentheses to demark quantifiers. This regexprep
% corrects those statements to the standard curly braces without making a
% global replacement of all parentheses.

out=regexprep(seq,'\((\d+,?\d*)\)','{$1}');

% guess the alphabet independently of the metacharacters
if ~isAminoAcid && (isNucleotide || bioinfoprivate.isnt(regexprep(seq, '[{(][\d,]+[)}]|\\\w', '')))
    iub = cell(11,1);
    if useAmbiguous
        iub{1,1}={'R','[AGR]'};
        iub{2,1}={'Y','[CTY]'};
        iub{3,1}={'M','[ACM]'};
        iub{4,1}={'K','[GTK]'};
        iub{5,1}={'S','[CGS]'};
        iub{6,1}={'W','[ATW]'};
        iub{7,1}={'B','[CGTYKSB]'};
        iub{8,1}={'D','[AGTRKWD]'};
        iub{9,1}={'H','[ACTYMWH]'};
        iub{10,1}={'V','[ACGRMSV]'};
        iub{11,1}={'N','[ACGTRYKMSWBDHVN]'};
    else
        iub{1,1}={'R','[AG]'};
        iub{2,1}={'Y','[CT]'};
        iub{3,1}={'M','[AC]'};
        iub{4,1}={'K','[GT]'};
        iub{5,1}={'S','[CG]'};
        iub{6,1}={'W','[AT]'};
        iub{7,1}={'B','[CGT]'};
        iub{8,1}={'D','[AGT]'};
        iub{9,1}={'H','[ACT]'};
        iub{10,1}={'V','[ACG]'};
        iub{11,1}={'N','[ACGT]'};
    end
    % leave gaps as -
    %iub{14,1}={'-','*'};

    for j=1:size(iub,1)
        patt = ['(?<!\\)' iub{j,1}{1,1}]; % look for non-metacharacters
        out = regexprep(out, patt, iub{j,1}{1,2});
    end
elseif isAminoAcid || bioinfoprivate.isaa(regexprep(seq, '[{(][\d,]+[)}]', ''))
    iub = cell(3,1);
    if useAmbiguous
        iub{1,1}={'B','[NDB]'};
        iub{2,1}={'Z','[QEZ]'};
        iub{3,1}={'X','[ARNDCQEGHILKMFPSTWYVBZX]'};
    else
        iub{1,1}={'B','[ND]'};
        iub{2,1}={'Z','[QE]'};
        iub{3,1}={'X','[ARNDCQEGHILKMFPSTWYV]'};
    end
    for j=1:size(iub,1)
        patt = ['(?<!\\)' iub{j,1}{1,1}]; % look for non-metacharacters
        out = regexprep(out, patt, iub{j,1}{1,2});
    end


end
