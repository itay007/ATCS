function  prop = aaburiedresidues(aa) 
%AABURIEDRESIDUES calculates molar fraction (%) of 2001 buried residues.
%
% Author(s) :  Janin J.
% Reference :  Nature 277:491-492(1979).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ % buried residues';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [11.200; ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        4.100;  ...% C  cysteine
        2.900;  ...% D  aspartic acid
        1.800;  ...% E  glutamic acid
        5.100;  ...% F  phenylalanine
        11.800; ...% G  glycine
        2.000;  ...% H  histidine
        8.600;  ...% I  isoleucine
        NaN;    ...% J  not used 
        0.500;  ...% K  lysine
        11.700; ...% L  leucine
        1.900;  ...% M  methionine
        2.900;  ...% N  asparagine
        NaN;    ...% O  not used 
        2.700;  ...% P  proline
        1.600;  ...% Q  glutamine
        0.500;  ...% R  arginine
        8.000;  ...% S  serine
        4.900;  ...% T  threonine
        NaN;    ...% U  not used 
        12.900; ...% V  valine
        2.200;  ...% W  tryptophan
        NaN;    ...% X  any amino acid
        2.600;  ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
