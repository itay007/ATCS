function  prop = aarelativemutability(aa) 
%AARELATIVEMUTABILITY calculates relative mutability of amino acids (Ala=100).
%
% Author(s) :  Dayhoff M.O., Schwartz R.M., Orcutt B.C.
% Reference :  In "Atlas of Protein Sequence and Structure", Vol.5, Suppl.3 (1978).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Relative mutability';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [100.000; ...% A  alanine
        NaN;     ...% B  aspartic acid or asparagine
        20.000;  ...% C  cysteine
        106.000; ...% D  aspartic acid
        102.000; ...% E  glutamic acid
        41.000;  ...% F  phenylalanine
        49.000;  ...% G  glycine
        66.000;  ...% H  histidine
        96.000;  ...% I  isoleucine
        NaN;     ...% J  not used 
        56.000;  ...% K  lysine
        40.000;  ...% L  leucine
        94.000;  ...% M  methionine
        134.000; ...% N  asparagine
        NaN;     ...% O  not used 
        56.000;  ...% P	 proline
        93.000;  ...% Q  glutamine
        65.000;  ...% R  arginine
        120.000; ...% S  serine
        97.000;  ...% T  threonine
        NaN;     ...% U  not used 
        74.000;  ...% V  valine
        18.000;  ...% W  tryptophan
        NaN;     ...% X  Any amino acid
        41.000;  ...% Y  tyrosine
        NaN;     ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
