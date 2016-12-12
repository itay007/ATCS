function  prop = aahplc2_1(aa) 
%AAHPLC2_1 calculates the retention coefficient in HPLC, pH 2.1.
%
% Author(s) :  Meek J.L.
% Reference :  Proc. Natl. Acad. Sci. USA 77:1632-1636(1980).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ HPLC retention, pH 2.1 (Meek)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [-0.100; ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        -2.200; ...% C  cysteine
        -2.800; ...% D  aspartic acid
        -7.500; ...% E  glutamic acid
        13.900; ...% F  phenylalanine
        -0.500; ...% G  glycine
        0.800;  ...% H  histidine
        11.800; ...% I  isoleucine
        NaN;    ...% J  not used 
        -3.200; ...% K  lysine
        10.000; ...% L  leucine
        7.100;  ...% M  methionine
        -1.600; ...% N  asparagine
        NaN;    ...% O  not used 
        8.000;  ...% P  proline
        -2.500; ...% Q  glutamine
        -4.500; ...% R  arginine
        -3.700; ...% S  serine
        1.500;  ...% T  threonine
        NaN;    ...% U  not used 
        3.300;  ...% V  valine
        18.100; ...% W  tryptophan
        NaN;    ...% X  any amino acid
        8.200;  ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
