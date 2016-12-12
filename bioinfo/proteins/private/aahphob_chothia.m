function  prop = aahphob_chothia(aa) 
%AAHPHOB_CHOTHIA calculates proportion of residues 95% buried (in 12 proteins).
%
% Author(s) :  Chothia C.
% Reference :  J. Mol. Biol. 105:1-14(1976).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Chothia)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.380; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        0.500; ...% C  cysteine
        0.150; ...% D  aspartic acid
        0.180; ...% E  glutamic acid
        0.500; ...% F  phenylalanine
        0.360; ...% G  glycine
        0.170; ...% H  histidine
        0.600; ...% I  isoleucine
        NaN;   ...% J  not used 
        0.030; ...% K  lysine
        0.450; ...% L  leucine
        0.400; ...% M  methionine
        0.120; ...% N  asparagine
        NaN;   ...% O  not used 
        0.180; ...% P  proline
        0.070; ...% Q  glutamine
        0.010; ...% R  arginine
        0.220; ...% S  serine
        0.230; ...% T  threonine
        NaN;   ...% U  not used 
        0.540; ...% V  valine
        0.270; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        0.150; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
