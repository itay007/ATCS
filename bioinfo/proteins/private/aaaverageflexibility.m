function  prop = aaaverageflexibility(aa) 
%AAAVERAGEFLEXIBILITY calculates the average flexibility index.
%
% Author(s) :  Bhaskaran R., Ponnuswamy P.K.
% Reference :  Int. J. Pept. Protein. Res. 32:242-255(1988).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Average flexibility';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.360; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        0.350; ...% C  cysteine
        0.510; ...% D  aspartic acid
        0.500; ...% E  glutamic acid
        0.310; ...% F  phenylalanine
        0.540; ...% G  glycine
        0.320; ...% H  histidine
        0.460; ...% I  isoleucine
        NaN;   ...% J  not used 
        0.470; ...% K  lysine
        0.370; ...% L  leucine
        0.300; ...% M  methionine
        0.460; ...% N  asparagine
        NaN;   ...% O  not used 
        0.510; ...% P  proline
        0.490; ...% Q  glutamine
        0.530; ...% R  arginine
        0.510; ...% S  serine
        0.440; ...% T  threonine
        NaN;   ...% U  not used 
        0.390; ...% V  valine
        0.310; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        0.420; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
