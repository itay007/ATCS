function  prop = aahphob_sweet(aa) 
%AAHPHOB_SWEET calculates optimized matching hydrophobicity (OMH).
%
% Author(s) :  Sweet R.M., Eisenberg D.
% Reference :  J. Mol. Biol. 171:479-488(1983).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Sweet et al.)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [-0.400; ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        0.170;  ...% C  cysteine
        -1.310; ...% D  aspartic acid
        -1.220; ...% E  glutamic acid
        1.920;  ...% F  phenylalanine
        -0.670; ...% G  glycine
        -0.640; ...% H  histidine
        1.250;  ...% I  isoleucine
        NaN;    ...% J  not used 
        -0.670; ...% K  lysine
        1.220;  ...% L  leucine
        1.020;  ...% M  methionine
        -0.920; ...% N  asparagine
        NaN;    ...% O  not used 
        -0.490; ...% P  proline
        -0.910; ...% Q  glutamine
        -0.590; ...% R  arginine
        -0.550; ...% S  serine
        -0.280; ...% T  threonine
        NaN;    ...% U  not used 
        0.910;  ...% V  valine
        0.500;  ...% W  tryptophan
        NaN;    ...% X  any amino acid
        1.670;  ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
