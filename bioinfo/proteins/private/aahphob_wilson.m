function  prop = aahphob_wilson(aa) 
%AAHPHOB_WILSON calculates the hydrophobic constants derived from HPLC peptide retention times.
%
% Author(s) :  Wilson K.J., Honegger A., Stotzel R.P., Hughes G.J.
% Reference :  Biochem. J. 199:31-41(1981).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity HPLC (Wilson et al.)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [-0.300; ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        6.300;  ...% C  cysteine
        -1.400; ...% D  aspartic acid
        0.000;  ...% E  glutamic acid
        7.500;  ...% F  phenylalanine
        1.200;  ...% G  glycine
        -1.300; ...% H  histidine
        4.300;  ...% I  isoleucine
        NaN;    ...% J  not used 
        -3.600; ...% K  lysine
        6.600;  ...% L  leucine
        2.500;  ...% M  methionine
        -0.200; ...% N  asparagine
        NaN;    ...% O  not used 
        2.200;  ...% P  proline
        -0.200; ...% Q  glutamine
        -1.100; ...% R  arginine
        -0.600; ...% S  serine
        -2.200; ...% T  threonine
        NaN;    ...% U  not used 
        5.900;  ...% V  valine
        7.900;  ...% W  tryptophan
        NaN;    ...% X  any amino acid
        7.100;  ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
