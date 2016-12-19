function nt = aa2nt(aa,varargin)
%AA2NT converts an amino acide sequence to a sequence of nucleotides.
%
%   AA2NT(SEQ) for a cell array calls aa2nt with individual elements of a
%   cell to convert amino acid sequences NT to nucleotide sequences.
%
%   See also AMINOLOOKUP, BASELOOKUP, GENETICCODE, NT2AA, REVGENETICCODE.

% Copyright 2005 The MathWorks, Inc.


nt = cell(size(aa));
for i = 1:numel(aa)
    nt{i} = aa2nt(aa{i},varargin{:});
end
