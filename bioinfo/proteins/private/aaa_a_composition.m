function  prop = aaa_a_composition(aa) 
%AAA_A_COMPOSITION calculates the overall amino acid composition (%).
%
% Author(s) :  McCaldon P., Argos P.
% Reference :  Proteins: Structure, Function and Genetics 4:99-122(1988).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Amino Acid Composition (%)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [8.300; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        1.700; ...% C  cysteine
        5.300; ...% D  aspartic acid
        6.200; ...% E  glutamic acid
        3.900; ...% F  phenylalanine
        7.200; ...% G  glycine
        2.200; ...% H  histidine
        5.200; ...% I  isoleucine
        NaN;   ...% J  not used 
        5.700; ...% K  lysine
        9.000; ...% L  leucine
        2.400; ...% M  methionine
        4.400; ...% N  asparagine
        NaN;   ...% O  not used 
        5.100; ...% P  proline
        4.000; ...% Q  glutamine
        5.700; ...% R  arginine
        6.900; ...% S  serine
        5.800; ...% T  threonine
        NaN;   ...% U  not used 
        6.600; ...% V  valine
        1.300; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        3.200; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
