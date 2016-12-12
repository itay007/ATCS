function  prop = aabeta_turnfasman(aa) 
%AABETA_TURNFASMAN calculates the conformational parameter for beta-turn.
%
%   AABETA_TURNFASMAN(AA) calculates the conformational parameter for
%   beta-turn. 
%
% Author(s) :  Chou P.Y., Fasman G.D.
% Reference :  Adv. Enzym. 47:45-148(1978).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Beta-turn (Chou & Fasman)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.660; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        1.190; ...% C  cysteine
        1.460; ...% D  aspartic acid
        0.740; ...% E  glutamic acid
        0.600; ...% F  phenylalanine
        1.560; ...% G  glycine
        0.950; ...% H  histidine
        0.470; ...% I  isoleucine
        NaN;   ...% J  not used 
        1.010; ...% K  lysine
        0.590; ...% L  leucine
        0.600; ...% M  methionine
        1.560; ...% N  asparagine
        NaN;   ...% O  not used 
        1.520; ...% P  proline
        0.980; ...% Q  glutamine
        0.950; ...% R  arginine
        1.430; ...% S  serine
        0.960; ...% T  threonine
        NaN;   ...% U  not used 
        0.500; ...% V  valine
        0.960; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        1.140; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
