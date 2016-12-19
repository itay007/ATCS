function  prop = aahphob_doolittle(aa) 
%AAHPHOB_DOOLITTLE calculates hydropathicity.
%
% Author(s) :  Kyte J., Doolittle R.F.
% Reference :  J. Mol. Biol. 157:105-132(1982).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Kyte & Doolittle)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [1.800;  ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        2.500;  ...% C  cysteine
        -3.500; ...% D  aspartic acid
        -3.500; ...% E  glutamic acid
        2.800;  ...% F  phenylalanine
        -0.400; ...% G  glycine
        -3.200; ...% H  histidine
        4.500;  ...% I  isoleucine
        NaN;    ...% J  not used 
        -3.900; ...% K  lysine
        3.800;  ...% L  leucine
        1.900;  ...% M  methionine
        -3.500; ...% N  asparagine
        NaN;    ...% O  not used 
        -1.600; ...% P  proline
        -3.500; ...% Q  glutamine
        -4.500; ...% R  arginine
        -0.800; ...% S  serine
        -0.700; ...% T  threonine
        NaN;    ...% U  not used 
        4.200;  ...% V  valine
        -0.900; ...% W  tryptophan
        NaN;    ...% X  any amino acid
        -1.300; ...% Y  tyrosine
        NaN;    ...% Z	glutamic acid or glutamine
];

prop = data(ndx);
