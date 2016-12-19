function  prop = aahplc7_4(aa) 
%AAHPLC7_4 calculates the retention coefficient in HPLC, pH 7.4.
%
% Author(s) :  Meek J.L.
% Reference :  Proc. Natl. Acad. Sci. USA 77:1632-1636(1980).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ HPLC retention, pH 7.4 (Meek)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.500;   ...% A  alanine
        NaN;     ...% B  aspartic acid or asparagine
        -6.800;  ...% C  cysteine
        -8.200;  ...% D  aspartic acid
        -16.900; ...% E  glutamic acid
        13.200;  ...% F  phenylalanine
        0.000;   ...% G  glycine
        -3.500;  ...% H  histidine
        13.900;  ...% I  isoleucine
        NaN;     ...% J  not used 
        0.100;   ...% K  lysine
        8.800;   ...% L  leucine
        4.800;   ...% M  methionine
        0.800;   ...% N  asparagine
        NaN;     ...% O  not used 
        6.100;   ...% P  proline
        -4.800;  ...% Q  glutamine
        0.800;   ...% R  arginine
        1.200;   ...% S  serine
        2.700;   ...% T  threonine
        NaN;     ...% U  not used 
        2.700;   ...% V  valine
        14.900;  ...% W  tryptophan
        NaN;     ...% X  any amino acid
        6.100;   ...% Y  tyrosine
        NaN;     ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
