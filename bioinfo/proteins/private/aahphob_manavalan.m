function  prop = aahphob_manavalan(aa) 
%AAHPHOB_MANAVALAN calculates average surrounding hydrophobicity.
%
% Author(s) :  Manavalan P., Ponnuswamy P.K.
% Reference :  Nature 275:673-674(1978).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Manavalan et al.)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [12.970; ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        14.630; ...% C  cysteine
        10.850; ...% D  aspartic acid
        11.890; ...% E  glutamic acid
        14.000; ...% F  phenylalanine
        12.430; ...% G  glycine
        12.160; ...% H  histidine
        15.670; ...% I  isoleucine
        NaN;    ...% J  not used 
        11.360; ...% K  lysine
        14.900; ...% L  leucine
        14.390; ...% M  methionine
        11.420; ...% N  asparagine
        NaN;    ...% O  not used 
        11.370; ...% P  proline
        11.760; ...% Q  glutamine
        11.720; ...% R  arginine
        11.230; ...% S  serine
        11.690; ...% T  threonine
        NaN;    ...% U  not used 
        15.710; ...% V  valine
        13.930; ...% W  tryptophan
        NaN;    ...% X  any amino acid
        13.420; ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
