function  prop = aahphob_parker(aa) 
%AAHPHOB_PARKER calculates hydrophilicity scale derived from HPLC peptide retention times.
%
% Author(s) :  Parker J.M.R., Guo D., Hodges R.S.
% Reference :  Biochemistry 25:5425-5431(1986).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity HPLC (Parker et al.)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [2.100;   ...% A  alanine
        NaN;     ...% B  aspartic acid or asparagine
        1.400;   ...% C  cysteine
        10.000;  ...% D  aspartic acid
        7.800;   ...% E  glutamic acid
        -9.200;  ...% F  phenylalanine
        5.700;   ...% G  glycine
        2.100;   ...% H  histidine
        -8.000;  ...% I  isoleucine
        NaN;     ...% J  not used 
        5.700;   ...% K  lysine
        -9.200;  ...% L  leucine
        -4.200;  ...% M  methionine
        7.000;   ...% N  asparagine
        NaN;     ...% O  not used 
        2.100;   ...% P  proline
        6.000;   ...% Q  glutamine
        4.200;   ...% R  arginine
        6.500;   ...% S  serine
        5.200;   ...% T  threonine
        NaN;     ...% U  not used 
        -3.700;  ...% V  valine
        -10.000; ...% W  tryptophan
        NaN;     ...% X  any amino acid
        -1.900;  ...% Y  tyrosine
        NaN;     ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
