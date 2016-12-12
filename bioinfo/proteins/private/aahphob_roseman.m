function  prop = aahphob_roseman(aa) 
%AAHPHOB_ROSEMAN calculates hydrophobicity scale (pi-r).
%
% Author(s) :  Roseman M.A.
% Reference :  J. Mol. Biol. 200:513-522(1988).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Roseman)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.390;  ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        0.250;  ...% C  cysteine
        -3.810; ...% D  aspartic acid
        -2.910; ...% E  glutamic acid
        2.270;  ...% F  phenylalanine
        0.000;  ...% G  glycine
        -0.640; ...% H  histidine
        1.820;  ...% I  isoleucine
        NaN;    ...% J  not used 
        -2.770; ...% K  lysine
        1.820;  ...% L  leucine
        0.960;  ...% M  methionine
        -1.910; ...% N  asparagine
        NaN;    ...% O  not used 
        0.990;  ...% P  proline
        -1.300; ...% Q  glutamine
        -3.950; ...% R  arginine
        -1.240; ...% S  serine
        -1.000; ...% T  threonine
        NaN;    ...% U  not used 
        1.300;  ...% V  valine
        2.130;  ...% W  tryptophan
        NaN;    ...% X  any amino acid
        1.470;  ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
