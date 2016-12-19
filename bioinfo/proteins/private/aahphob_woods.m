function  prop = aahphob_woods(aa) 
%AAHPHOB_WOODS calculates the hydrophilicity (Hopp & Woods).
%
% Author(s) :  Hopp T.P., Woods K.R.
% Reference :  Proc. Natl. Acad. Sci. U.S.A. 78:3824-3828(1981).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Hopp & Woods) ';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [-0.500; ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        -1.000; ...% C  cysteine
        3.000;  ...% D  aspartic acid
        3.000;  ...% E  glutamic acid
        -2.500; ...% F  phenylalanine
        0.000;  ...% G  glycine
        -0.500; ...% H  histidine
        -1.800; ...% I  isoleucine
        NaN;    ...% J  not used 
        3.000;  ...% K  lysine
        -1.800; ...% L  leucine
        -1.300; ...% M  methionine
        0.200;  ...% N  asparagine
        NaN;    ...% O  not used 
        0.000;  ...% P  proline
        0.200;  ...% Q  glutamine
        3.000;  ...% R  arginine
        0.300;  ...% S  serine
        -0.400; ...% T  threonine
        NaN;    ...% U  not used 
        -1.500; ...% V  valine
        -3.400; ...% W  tryptophan
        NaN;    ...% X  any amino acid
        -2.300; ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
