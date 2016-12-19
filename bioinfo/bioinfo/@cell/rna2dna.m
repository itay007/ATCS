function dna = rna2dna(rna,varargin)
% RNA2DNA converts an RNA sequence into a DNA sequence.
%
%   DNA = RNA2DNA(RNA) converts any uracil nucleotides in an RNA sequence
%   into thymine (U-->T).
%
%   DNA is returned in the same format as RNA, so if RNA is an integer
%   sequence then so is DNA.
%
%   See also DNA2RNA, REGEXP, STRREP.

% Copyright 2005 The MathWorks, Inc.


dna = cell(size(rna));
for i = 1:numel(rna)
    dna{i} = rna2dna(rna{i},varargin{:});
end
