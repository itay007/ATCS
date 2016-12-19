function  prop = aahphob_janin(aa) 
%AAHPHOB_JANIN calculates free energy of transfer from inside to outside of a globular protein.
%
% Author(s) :  Janin J.
% Reference :  Nature 277:491-492(1979).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Janin)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.300;  ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        0.900;  ...% C  cysteine
        -0.600; ...% D  aspartic acid
        -0.700; ...% E  glutamic acid
        0.500;  ...% F  phenylalanine
        0.300;  ...% G  glycine
        -0.100; ...% H  histidine
        0.700;  ...% I  isoleucine
        NaN;    ...% J  not used 
        -1.800; ...% K  lysine
        0.500;  ...% L  leucine
        0.400;  ...% M  methionine
        -0.500; ...% N  asparagine
        NaN;    ...% O  not used 
        -0.300; ...% P  proline
        -0.700; ...% Q  glutamine
        -1.400; ...% R  arginine
        -0.100; ...% S  serine
        -0.200; ...% T  threonine
        NaN;    ...% U  not used 
        0.600;  ...% V  valine
        0.300;  ...% W  tryptophan
        NaN;    ...% X  any amino acid
        -0.400; ...% Y  tyrosine
        NaN;    ...% Z	glutamic acid or glutamine
];

prop = data(ndx);
