function nt = aa2nt(aa,varargin)
%AA2NT converts an amino acid sequence to a sequence of nucleotides.
%
%   AA2NT(SEQ) converts amino acid sequence SEQ to nucleotides using
%   the Standard genetic code. In general, the mapping from amino acid to
%   nucleotide codon is not a one-to-one mapping so the function chooses a
%   random codon corresponding to a particular amino acid.
%
%   AA2NT(...,'GENETICCODE',CODE) converts a nucleotide sequence to an
%   amino acid sequence using the genetic code CODE. CODE can be a string
%   or an ID number from the list below, or a structure created using the
%   function REVGENETICCODE. When a text string name is used, it can be
%   truncated to the first two characters of the name. The Standard genetic
%   code (ID = 1) is used by default.
%
%   ID  Name
%
% 	1	Standard
% 	2	Vertebrate Mitochondrial
% 	3	Yeast Mitochondrial
% 	4	Mold, Protozoan, and Coelenterate Mitochondrial and Mycoplasma/Spiroplasma
% 	5	Invertebrate Mitochondrial
% 	6	Ciliate, Dasycladacean, and Hexamita Nuclear
% 	9	Echinoderm Mitochondrial
% 	10	Euplotid Nuclear
% 	11	Bacterial and Plant Plastid
% 	12	Alternative Yeast Nuclear
% 	13	Ascidian Mitochondrial
% 	14	Flatworm Mitochondrial
% 	15	Blepharisma Nuclear
% 	16	Chlorophycean Mitochondrial
% 	21	Trematode Mitochondrial
% 	22	Scenedesmus Obliquus Mitochondrial
% 	23	Thraustochytrium Mitochondrial
%
%   AA2NT(...,'ALPHABET',A) defines which nucleotide alphabet to use. The
%   default value is 'DNA' which uses the symbols A,C,T,G. If ALPHABET is
%   set to 'RNA', then A,C,U,G are used instead.
%
%   Example:
%
%           aa2nt('MATLAP')
%
%   See also GENETICCODE, NT2AA, RAND, REVGENETICCODE, SEQVIEWER.

%   Copyright 2002-2012 The MathWorks, Inc.



% If the input is a structure then extract the Sequence data.
if isstruct(aa)
    aa = bioinfoprivate.seqfromstruct(aa);
end

code = 1;
rna = false;
if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:aa2nt:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'alphabet','geneticcode'};
    for j=1:2:nargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:aa2nt:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:aa2nt:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % alphabet
                    [isAminoAcid,rna] = bioinfoprivate.optAlphabet(pval,okargs{k}, mfilename);
                    if isAminoAcid
                        error(message('bioinfo:aa2nt:InvalidAlphabet'));
                    end
                case 2  % geneticcode
                    code = pval;
                    if isstruct(code)
                        % allow geneticcode created struct
                        numFields = numel(fieldnames(code));
                        if numFields == 66 && isfield(code, 'Name')
                            code = revgeneticcode(code.Name);
                            numFields = numel(fieldnames(code));
                        end
                        if numFields ~= 23 || ~isfield(code, 'Name')
                            error(message('bioinfo:aa2nt:BadCodeStruct'));
                        end
                    end
            end
        end
    end

end

% send individual rows for char array
if size(aa,1) > 1
    for i = size(aa,1):-1:1 %preallocates the resulting array
        nt(i,:) = aa2nt(aa(i,:),varargin{:});
    end
else
    if ~ischar(aa)
        intFlag = true;
        aa = int2aa(aa);
    else
        aa = upper(aa);
        intFlag = false;
    end
    aa = strrep(aa,'*','s');
    aa = strrep(aa,'?','u');
    aa = strrep(aa,'-','g');
    
    % B	Asp or Asn	D|N
    % Z	Glu or Gln	E|Q
    % X	   ANY		ANY
    ambiguous = regexp(aa,'[BZX]');
    if ~isempty(ambiguous)
        randoms = rand(size(ambiguous));
        warning(message('bioinfo:aa2nt:AmbiguousCharacters'));
        for count = 1:length(ambiguous)
            switch aa(ambiguous(count))
                case 'B'
                    if randoms(count) < 0.5
                        aa(ambiguous(count)) = 'D';
                    else
                        aa(ambiguous(count)) = 'N';
                    end
                case 'Z'
                    if randoms(count) < 0.5
                        aa(ambiguous(count)) = 'E';
                    else
                        aa(ambiguous(count)) = 'Q';
                    end
                case 'X'
                    aa(ambiguous(count)) = upper(int2aa(ceil(20*randoms(count))));
            end
        end
    end

    if ~isstruct(code)
        code = revgeneticcode(code);
    end
    code.s = code.Stops;
    code.u = {'???'};
    code.g = {'---'};
    
    numAAs = length(aa);
    numNts = 3 * numAAs;
    nt = blanks(numNts);
    randoms = rand(1,numAAs);
    place = 1;
    for count = 1:numAAs
        nt(place:place+2) = code.(aa(count)){ceil(randoms(count)*size(code.(aa(count)),2))};
        place = place + 3;
    end

    if rna == true
        nt = strrep(nt,'T','U');
    end
    if intFlag == true
        nt = nt2int(nt);
    end
end
