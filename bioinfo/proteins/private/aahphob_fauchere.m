function  prop = aahphob_fauchere(aa) 
%AAHPHOB_FAUCHERE calculates hydrophobicity scale (pi-r).
%
% Author(s) :  Fauchere J.-L., Pliska V.E.
% Reference :  Eur. J. Med. Chem. 18:369-375(1983).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Fauchere & Pliska)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.310;  ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        1.540;  ...% C  cysteine
        -0.770; ...% D  aspartic acid
        -0.640; ...% E  glutamic acid
        1.790;  ...% F  phenylalanine
        0.000;  ...% G  glycine
        0.130;  ...% H  histidine
        1.800;  ...% I  isoleucine
        NaN;    ...% J  not used
        -0.990; ...% K  lysine
        1.700;  ...% L  leucine
        1.230;  ...% M  methionine
        -0.600; ...% N  asparagine
        NaN;    ...% O  not used 
        0.720;  ...% P  proline
        -0.220; ...% Q  glutamine
        -1.010; ...% R  arginine
        -0.040; ...% S  serine
        0.260;  ...% T  threonine
        NaN;    ...% U  not used 
        1.220;  ...% V  valine
        2.250;  ...% W  tryptophan
        NaN;    ...% X  any amino acid
        0.960;  ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
