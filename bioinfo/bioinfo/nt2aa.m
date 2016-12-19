function aa = nt2aa(nt,varargin)
% NT2AA converts a nucleotide sequence to a sequence of amino acids.
%
%   NT2AA(SEQ) converts nucleotide sequence SEQ to an amino acid sequence
%   using the Standard genetic code.
%
%   NT2AA(...,'FRAME',RF) converts a nucleotide sequence for the reading
%   frame RF to an amino acid sequence. Set RF to 'All' to convert all
%   three reading frames. In this case, the output is a 3x1 cell array. The
%   default FRAME is 1.
%
%   NT2AA(...,'ALTERNATIVESTARTCODONS',TF) is used to control the use of
%   alternative start codons. By default, if the first codon of a sequence
%   corresponds to a known alternative start codon, the codon is translated
%   to methionine. If this option is set to false, then alternative start
%   codons at the start of a sequence are translated to their corresponding
%   amino acids for the genetic code that is being used, which might not
%   necessarily be methionine.
%
%   See http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=t#SG1
%   for more details of alternative start codons.
%
%   NT2AA(...,'GENETICCODE',CODE) converts a nucleotide sequence to an
%   amino acid sequence using the genetic code CODE. CODE can be a string
%   or an ID number from the list below, or a structure created using the
%   function GENETICCODE. When a text string name is used, it can be
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
%   NT2AA(...,'ACGTOnly',TF) allows you to specify the behavior of
%   ambiguous nucleotide characters (R, Y, K, M, S, W, B, D, H, V and N)
%   and unknown characters. If TF is true (default) then the function will
%   error if any of these characters are present. If TF is false then the
%   function will attempt to resolve ambiguities. If it cannot it will
%   return X for the affected codon.
%
%   Note that gaps are allowed only if the entire codon consists of gaps.
%
%   Example:
%
%           nt2aa('ATGGCTACGCTAGCTCCT')
%
%   See also AA2NT, ALIGNDEMO, AMINOLOOKUP, BASELOOKUP, CODONBIAS, DNDS,
%   DNDSML, GENETICCODE, REVGENETICCODE, SEQVIEWER.

%   Copyright 2002-2012 The MathWorks, Inc.


% check inputs
bioinfochecknargin(nargin,1,mfilename);
[frame, numFrames, outcell, code, alternativeStart,acgtOnly] = parse_inputs(varargin{:});
% If the input is a structure then extract the Sequence data.
if isstruct(nt)
    nt = bioinfoprivate.seqfromstruct(nt);
end

% send individual rows for char array
if size(nt,1) > 1
    for i = size(nt,1):-1:1 %preallocates the resulting array
        aa(i,:) = nt2aa(nt(i,:),varargin{:});
    end
else
    
    % convert to chars for the lookup
    if ~ischar(nt)
        intFlag = true;
        nt = int2nt(nt);
    else
        nt = upper(nt);
        intFlag = false;
        % convert rna2dna
        nt = strrep(nt,'U','T');
    end
    
    % get our genetic code
    if ~isstruct(code)
        code = geneticcode(code);
    end
    %
    seqLen = length(nt);
    try
        for frameNum = frame
            
            numNts = seqLen +1 -frameNum;
            numCodons = floor(numNts/3);
            aa = blanks(numCodons);
            codon = 0;
            % loop through the codons looking up the AAs as we go
            for count =frameNum:3:(seqLen-2)
                codon = codon + 1;
                triplet = nt(count:count+2);
                if any(triplet == '---')
                    if isequal(triplet,'---')
                        aa(codon) = '-';
                    else
                        error(message('bioinfo:nt2aa:GapInCodon', triplet));
                    end
                else
                    try
                        aa(codon) = code.(triplet);
                    catch badTripletException
                        if acgtOnly
                            rethrow(badTripletException)
                        else
                            aa(codon) = checkAmbiguous(triplet,code);
                        end
                    end
                end
            end
            % deal with alternative start codons
            if alternativeStart && (seqLen >= frameNum+2)
                if ismember(nt(frameNum:frameNum+2),code.Starts)
                    aa(1) = 'M';
                end
            end
            
            % convert back to ints if necessary
            if intFlag == true
                aa = aa2int(aa);
            end
            
            if numFrames > 1
                outcell{frameNum} = aa;
            end
            
        end
    catch err
        % rethrow the any expected errors
        if isequal(err.identifier,'bioinfo:nt2aa:GapInCodon')
            rethrow(err);
        end
        % Most likely problem is that we have characters other than ACGTU.
        if acgtOnly && bioinfoprivate.isnt(nt,'ACGTUOnly',true) == false
            if isnumeric(nt)
                try
                    nt = int2nt(nt);
                catch allExceptions %#ok<NASGU>
                    % if we fail then set something so that the error makes some
                    % sense.
                    nt = '?';
                end
            end
            badchars = setxor(unique([upper(nt) 'AGCTU']),'AGCTU'); % add ACGT just to be safe
            if ~isequal(badchars,'-')
                error(message('bioinfo:nt2aa:ACGTOnly', badchars));
            end
        end
        % If we get here then just rethrow the error
        rethrow(err);
    end
    
    % if we have multi frame output then use a cell array
    if numFrames > 1
        aa = outcell;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  codon = checkAmbiguous(triplet,code)
%   When acgtOnly is false we need to handle these cases:
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
%   And map any that are still ambiguous to
%
%   B --> [DN] 	aspartic acid or asparagine
%   Z --> [EQ]	glutamic acid or glutamine
%   X --> [ARNDCQEGHILKMFPSTWYV]

% use seq2regexp to get all possible sequences
for count = 3:-1:1
    nucs{count} = seq2regexp(triplet(count),'AMBIGUOUS',false);
    nucs{count} = strrep(nucs{count},'[','');
    nucs{count} = strrep(nucs{count},']','');
end

% get the first
codon = code.([nucs{1}(1),nucs{2}(1),nucs{3}(1)]);
% now loop through all possible cases until we find
for i = 1:numel(nucs{3})
    for j = 1:numel(nucs{2})
        for k = 1:numel(nucs{1})
            try
                theCodon = code.([nucs{1}(k),nucs{2}(j),nucs{3}(i)]);
                if ~isequal(theCodon,codon)
                    if all(ismember({theCodon,codon},{'B','D','N'}))
                        codon = 'B';
                    elseif all(ismember({theCodon,codon},{'Z','E','Q'}))
                        codon = 'Z';
                    else
                        codon = 'X';
                        return
                    end
                end
            catch allExceptions %#ok<NASGU>
                codon = 'X';
                return
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [frame, numFrames, outcell, code, alternativeStart,acgtOnly] = parse_inputs(varargin)
% Handle inputs

% Check that we have the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:nt2aa:IncorrectNumberOfArguments', mfilename));
end

% Valid inputs
okargs = {'frame','geneticcode','alternativestartcodons','acgtonly'};

% Set defaults
code = 1;
frame = 1;
numFrames = 1;
alternativeStart = true;
acgtOnly = true;
outcell = {};

% Loop over the values
for j=1:2:nargin
    % Lookup the pair
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % frame
            if ~isnumeric(pval) || max(pval) > 3 || min(pval) < 1
                if ischar(pval) && strcmpi(pval,'all')
                    pval = [1, 2, 3];
                else
                    error(message('bioinfo:nt2aa:BadFrameNumber'))
                end
            end
            frame = pval;
            numFrames = length(frame);
            if length(frame) > 1
                outcell = cell(3,1);
            end
        case 2  % genetic code
            code = pval;
            if isstruct(code)
                numFields = numel(fieldnames(code));
                % allow revgeneticcode created struct
                if numFields == 23 && isfield(code, 'Name')
                    code = geneticcode(code.Name);
                    numFields = numel(fieldnames(code));
                end
                if numFields ~= 66 || ~isfield(code, 'Name')
                    error(message('bioinfo:nt2aa:BadCodeStruct'));
                end
            end
        case 3  % alternative start
            alternativeStart = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 4  % acgtOnly
            acgtOnly = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end
end
