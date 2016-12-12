function aa = nt2aa(nt,varargin)
% NT2AA converts a nucleotide sequence to a sequence of amino acids.
%
%   NT2AA(SEQ) for a cell array calls nt2aa with individual elements of a
%   cell to convert nucleotide sequences NT to amino acid sequences.
%
%   See also AA2NT, AMINOLOOKUP, BASELOOKUP, GENETICCODE, REVGENETICCODE.

% Copyright 2005 The MathWorks, Inc.


aa = cell(size(nt));
for i = 1:numel(nt)
    aa{i} = nt2aa(nt{i},varargin{:});
end
