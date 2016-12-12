function  prop = aahphob_mobility(aa) 
%AAHPHOB_MOBILITY calculates mobilities of amino acids on chromatography paper (RF).
%
% Author(s) :  Aboderin A.A.
% Reference :  Int. J. Biochem. 2:537-544(1971).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Aboderin)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [5.100;  ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        0.000;  ...% C  cysteine
        0.700;  ...% D  aspartic acid
        1.800;  ...% E  glutamic acid
        9.600;  ...% F  phenylalanine
        4.100;  ...% G  glycine
        1.600;  ...% H  histidine
        9.300;  ...% I  isoleucine
        NaN;    ...% J  not used 
        1.300;  ...% K  lysine
        10.000; ...% L  leucine
        8.700;  ...% M  methionine
        0.600;  ...% N  asparagine
        NaN;    ...% O  not used 
        4.900;  ...% P  proline
        1.400;  ...% Q  glutamine
        2.000;  ...% R  arginine
        3.100;  ...% S  serine
        3.500;  ...% T  threonine
        NaN;    ...% U  not used 
        8.500;  ...% V  valine
        9.200;  ...% W  tryptophan
        NaN;    ...% X  any amino acid
        8.000;  ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];
prop = data(ndx);
