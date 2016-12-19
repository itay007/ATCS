function rna = dna2rna(dna,varargin)
% DNA2RNA converts a DNA sequence into an RNA sequence.
%
%   RNA = DNA2RNA(DNA) converts any thymine nucleotides in a DNA sequence
%   into uracil (T-->U).
%
%   RNA is returned in the same format as DNA, so if DNA is an integer
%   sequence then so is RNA.
%
%   See also REGEXP, RNA2DNA, STRREP.

% Copyright 2005 The MathWorks, Inc.


rna = cell(size(dna));
for i = 1:numel(dna)
    rna{i} = dna2rna(dna{i},varargin{:});
end
