function  prop = aahplchfba(aa) 
%AAHPLCHFBA calculates the retention coefficient in HFBA.
%
% Author(s) :  Browne C.A., Bennett H.P.J., Solomon S.
% Reference :  Anal. Biochem. 124:201-208(1982).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ HFBA retention';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [3.900;   ...% A alanine
        NaN;     ...% B  aspartic acid or asparagine
        -14.300; ...% C  cysteine
        -2.800;  ...% D  aspartic acid
        -7.500;  ...% E  glutamic acid
        14.700;  ...% F  phenylalanine
        -2.300;  ...% G  glycine
        2.000;   ...% H  histidine
        11.000;  ...% I  isoleucine
        NaN;     ...% J  not used 
        -2.500;  ...% K  lysine
        15.000;  ...% L  leucine
        4.100;   ...% M  methionine
        -2.800;  ...% N  asparagine
        NaN;     ...% O  not used 
        5.600;   ...% P  proline
        1.800;   ...% Q  glutamine
        3.200;   ...% R  arginine
        -3.500;  ...% S  serine
        1.100;   ...% T  threonine
        NaN;     ...% U  not used 
        2.100;   ...% V  valine
        17.800;  ...% W  tryptophan
        NaN;     ...% X  any amino acid
        3.800;   ...% Y  tyrosine
        NaN;     ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
