function  prop = aabeta_turnroux(aa) 
%AABETA_TURNROUX calculates the conformational parameter for beta-turn.
%
% Author(s) :  Deleage G., Roux B.
% Reference :  Protein Engineering 1:289-294(1987).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Beta-turn (Roux)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.788; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        0.965; ...% C  cysteine
        1.197; ...% D  aspartic acid
        1.149; ...% E  glutamic acid
        0.624; ...% F  phenylalanine
        1.860; ...% G  glycine
        0.970; ...% H  histidine
        0.240; ...% I  isoleucine
        NaN;   ...% J  not used 
        1.302; ...% K  lysine
        0.670; ...% L  leucine
        0.436; ...% M  methionine
        1.572; ...% N  asparagine
        NaN;   ...% O  not used 
        1.415; ...% P  proline
        0.997; ...% Q  glutamine
        0.912; ...% R  arginine
        1.316; ...% S  serine
        0.739; ...% T  threonine
        NaN;   ...% U  not used 
        0.387; ...% V  valine
        0.546; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        0.795; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
