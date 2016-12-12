function  prop = aahphob_ph7_5(aa) 
%AAHPHOB_PH7_5 calculates hydrophobicity indices at ph 7.5 determined by HPLC.
%
% Author(s) :  Cowan R., Whittaker R.G.
% Reference :  Peptide Research 3:75-80(1990).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hphob. HPLC pH7.5';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.350;  ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        0.760;  ...% C  cysteine
        -2.150; ...% D  aspartic acid
        -1.950; ...% E  glutamic acid
        1.690;  ...% F  phenylalanine
        0.000;  ...% G  glycine
        -0.650; ...% H  histidine
        1.830;  ...% I  isoleucine
        NaN;    ...% J  not used 
        -1.540; ...% K  lysine
        1.800;  ...% L  leucine
        1.100;  ...% M  methionine
        -0.990; ...% N  asparagine
        NaN;    ...% O  not used 
        0.840;  ...% P  proline
        -0.930; ...% Q  glutamine
        -1.500; ...% R  arginine
        -0.630; ...% S  serine
        -0.270; ...% T  threonine
        NaN;    ...% U  not used 
        1.320;  ...% V  valine
        1.350;  ...% W  tryptophan
        NaN;    ...% X  any amino acid
        0.390;  ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
