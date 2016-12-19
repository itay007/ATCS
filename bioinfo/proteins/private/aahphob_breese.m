function  prop = aahphob_breese(aa) 
%AAHPHOB_BREESE calculates hydrophobicity (free energy of transfer to surface in kcal/mole).
%
% Author(s) :  Bull H.B., Breese K.
% Reference :  Arch. Biochem. Biophys. 161:665-670(1974).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Bull & Breese)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.610;  ...% A	alanine
        NaN;    ...% B  aspartic acid or asparagine
        0.360;  ...% C  cysteine
        0.610;  ...% D  aspartic acid
        0.510;  ...% E  glutamic acid
        -1.520; ...% F  phenylalanine
        0.810;  ...% G  glycine
        0.690;  ...% H  histidine
        -1.450; ...% I  isoleucine
        NaN;    ...% J  not used 
        0.460;  ...% K  lysine
        -1.650; ...% L  leucine
        -0.660; ...% M  methionine
        0.890;  ...% N  asparagine
        NaN;    ...% O  not used 
        -0.170; ...% P  proline
        0.970;  ...% Q  glutamine
        0.690;  ...% R  arginine
        0.420;  ...% S  serine
        0.290;  ...% T  threonine
        NaN;    ...% U  not used 
        -0.750; ...% V  valine
        -1.200; ...% W  tryptophan
        NaN;    ...% X  any amino acid
        -1.430; ...% Y  tyrosine
        NaN;    ...% Z	glutamic acid or glutamine
];

prop = data(ndx);
