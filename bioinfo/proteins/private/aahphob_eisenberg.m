function  prop = aahphob_eisenberg(aa) 
%AAHPHOB_EISENBERG calculates normalized consensus hydrophobicity scale.
%
% Author(s) :  Eisenberg D., Schwarz E., Komarony M., Wall R.
% Reference :  J. Mol. Biol. 179:125-142(1984).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Eisenberg et al.) ';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.620;  ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        0.290;  ...% C  cysteine
        -0.900; ...% D  aspartic acid
        -0.740; ...% E  glutamic acid
        1.190;  ...% F  phenylalanine
        0.480;  ...% G  glycine
        -0.400; ...% H  histidine
        1.380;  ...% I  isoleucine
        NaN;    ...% J  not used 
        -1.500; ...% K  lysine
        1.060;  ...% L  leucine
        0.640;  ...% M  methionine
        -0.780; ...% N  asparagine
        NaN;    ...% O  not used 
        0.120;  ...% P  proline
        -0.850; ...% Q  glutamine
        -2.530; ...% R  arginine
        -0.180; ...% S  serine
        -0.050; ...% T  threonine
        NaN;    ...% U  not used 
        1.080;  ...% V  valine
        0.810;  ...% W  tryptophan
        NaN;    ...% X  any amino acid
        0.260;  ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
