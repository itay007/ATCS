function  prop = aahphob_guy(aa) 
%AAHPHOB_GUY calculates hydrophobicity scale based on free energy of transfer (kcal/mole). 
%
% Author(s) :  Guy H.R.
% Reference :  Biophys J. 47:61-70(1985).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Guy)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.100;  ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        -1.420; ...% C  cysteine
        0.780;  ...% D  aspartic acid
        0.830;  ...% E  glutamic acid
        -2.120; ...% F  phenylalanine
        0.330;  ...% G  glycine
        -0.500; ...% H  histidine
        -1.130; ...% I  isoleucine
        NaN;    ...% J  not used 
        1.400;  ...% K  lysine
        -1.180; ...% L  leucine
        -1.590; ...% M  methionine
        0.480;  ...% N  asparagine
        NaN;    ...% O  not used 
        0.730;  ...% P  proline
        0.950;  ...% Q  glutamine
        1.910;  ...% R  arginine
        0.520;  ...% S  serine
        0.070;  ...% T  threonine
        NaN;    ...% U  not used 
        -1.270; ...% V  valine
        -0.510; ...% W  tryptophan
        NaN;    ...% X  any amino acid
        -0.210; ...% Y  tyrosine
        NaN;    ...% Z	glutamic acid or glutamine
];

    prop = data(ndx);
