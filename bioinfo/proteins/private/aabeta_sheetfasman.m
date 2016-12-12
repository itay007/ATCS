function  prop = aabeta_sheetfasman(aa) 
%AABETA_SHEETFASMAN calculates the conformational parameter for beta-sheet.
%
%   AABETA_SHEETFASMAN calculates the conformational parameter for
%   beta-sheet (computed from 29 proteins).   
%
% Author(s) :  Chou P.Y., Fasman G.D.
% Reference :  Adv. Enzym. 47:45-148(1978).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Beta-sheet (Chou & Fasman)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.830; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        1.190; ...% C  cysteine
        0.540; ...% D  aspartic acid
        0.370; ...% E  glutamic acid
        1.380; ...% F  phenylalanine
        0.750; ...% G  glycine
        0.870; ...% H  histidine
        1.600; ...% I  isoleucine
        NaN;   ...% J  not used 
        0.740; ...% K  lysine
        1.300; ...% L  leucine
        1.050; ...% M  methionine
        0.890; ...% N  asparagine
        NaN;   ...% O  not used 
        0.550; ...% P  proline
        1.100; ...% Q  glutamine
        0.930; ...% R  arginine
        0.750; ...% S  serine
        1.190; ...% T  threonine
        NaN;   ...% U  not used 
        1.700; ...% V  valine
        1.370; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        1.470; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
