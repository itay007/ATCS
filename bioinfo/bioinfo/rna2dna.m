function dna = rna2dna(rna)
% RNA2DNA converts an RNA sequence into a DNA sequence.
%
%   DNA = RNA2DNA(RNA) converts any uracil nucleotides in an RNA sequence
%   into thymine (U-->T).
%
%   DNA is returned in the same format as RNA, so if RNA is an integer
%   sequence then so is DNA.
%
%   Example:
%
%       rna2dna('ACGAUGAGUCAUGCUU')
%
%   See also DNA2RNA, REGEXP, STRREP.

%   Copyright 2002-2012 The MathWorks, Inc.


% If the input is a structure then extract the Sequence data.
if isstruct(rna)
    rna = bioinfoprivate.seqfromstruct(rna);
end
if ~ischar(rna)
    dna = rna;
    return
end

if ~isempty(regexpi(rna,'t','once'))
    warning(message('bioinfo:rna2dna:RNAContainsT'));
end

dna = strrep(rna,'U','T');
dna = strrep(dna,'u','t');

if ~isempty(regexpi(dna,'[^ACGTRYKMSWBDHVN-]','once'))
    warning(message('bioinfo:rna2dna:UnknownSymbols'));
end



