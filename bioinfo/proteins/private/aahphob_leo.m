function  prop = aahphob_leo(aa) 
%AAHPHOB_LEO calculates hydrophobicity (delta G1/2 cal)
%
% Author(s) :  Abraham D.J., Leo A.J. 
% Reference :  Proteins: Structure, Function and Genetics 2:130-152(1987).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Abraham & Leo)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.440;  ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        0.580;  ...% C  cysteine
        -0.310; ...% D  aspartic acid
        -0.340; ...% E  glutamic acid
        2.540;  ...% F  phenylalanine
        0.000;  ...% G  glycine
        -0.010; ...% H  histidine
        2.460;  ...% I  isoleucine
        NaN;    ...% J  not used 
        -2.450; ...% K  lysine
        2.460;  ...% L  leucine
        1.100;  ...% M  methionine
        -1.320; ...% N  asparagine
        NaN;    ...% O  not used 
        1.290;  ...% P  proline
        -0.710; ...% Q  glutamine
        -2.420; ...% R  arginine
        -0.840; ...% S  serine
        -0.410; ...% T  threonine
        NaN;    ...% U  not used 
        1.730;  ...% V  valine
        2.560;  ...% W  tryptophan
        NaN;    ...% X  any amino acid
        1.630;  ...% Y  tyrosine
        NaN;    ...% Z	glutamic acid or glutamine
];

prop = data(ndx);
