function  prop = aapolarityzimmerman(aa) 
%AAPOLARITYZIMMERMAN -  Polarity.
%
% Author(s) :  Zimmerman J.M., Eliezer N., Simha R.
% Reference :  J. Theor. Biol. 21:170-201(1968).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Polarity (Zimmerman)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.000;  ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        1.480;  ...% C  cysteine
        49.700; ...% D  aspartic acid
        49.900; ...% E  glutamic acid
        0.350;  ...% F  phenylalanine
        0.000;  ...% G  glycine
        51.600; ...% H  histidine
        0.130;  ...% I  isoleucine
        NaN;    ...% J  not used 
        49.500; ...% K  lysine
        0.130;  ...% L  leucine
        1.430;  ...% M  methionine
        3.380;  ...% N  asparagine
        NaN;    ...% O  not used 
        1.580;  ...% P  proline
        3.530;  ...% Q  glutamine
        52.000; ...% R  arginine
        1.670;  ...% S  serine
        1.660;  ...% T  threonine
        NaN;    ...% U  not used 
        0.130;  ...% V  valine
        2.100;  ...% W  tryptophan
        NaN;    ...% X  any amino acid
        1.610;  ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
