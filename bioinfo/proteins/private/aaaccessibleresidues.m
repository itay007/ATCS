function  prop = aaaccessibleresidues(aa) 
%AAACCESSIBLERESIDUES calculates the molar fraction (%) of 3220 accessible residues.
%
% Author(s) :  Janin J.
% Reference :  Nature 277:491-492(1979).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ % accessible residues';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [6.600; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        0.900; ...% C  cysteine
        7.700; ...% D  aspartic acid
        5.700; ...% E  glutamic acid
        2.400; ...% F  phenylalanine
        6.700; ...% G  glycine
        2.500; ...% H  histidine
        2.800; ...% I  isoleucine
        NaN;   ...% J  not used 
        10.300;...% K  lysine
        4.800; ...% L  leucine
        1.000; ...% M  methionine
        6.700; ...% N  asparagine
        NaN;   ...% O  not used 
        4.800; ...% P  proline
        5.200; ...% Q  glutamine
        4.500; ...% R  arginine
        9.400; ...% S  serine
        7.000; ...% T  threonine
        NaN;   ...% U  not used 
        4.500; ...% V  valine
        1.400; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        5.100; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
