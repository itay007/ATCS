function  prop = aahplctfa(aa) 
%AAHPLCTFA calculates the retention coefficient in TFA.
%
% Author(s) :  Browne C.A., Bennett H.P.J., Solomon S.
% Reference :  Anal. Biochem. 124:201-208(1982).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ TFA retention';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [7.300;  ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        -9.200; ...% C  cysteine
        -2.900; ...% D  aspartic acid
        -7.100; ...% E  glutamic acid
        19.200; ...% F  phenylalanine
        -1.200; ...% G  glycine
        -2.100; ...% H  histidine
        6.600;  ...% I  isoleucine
        NaN;    ...% J   not used 
        -3.700; ...% K  lysine
        20.000; ...% L  leucine
        5.600;  ...% M  methionine
        -5.700; ...% N  asparagine
        NaN;    ...% O  not used 
        5.100;  ...% P  proline
        -0.300; ...% Q  glutamine
        -3.600; ...% R  arginine
        -4.100; ...% S  serine
        0.800;  ...% T  threonine
        NaN;    ...% U  not used 
        3.500;  ...% V  valine
        16.300; ...% W  tryptophan
        NaN;    ...% X  any amino acid
        5.900;  ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
