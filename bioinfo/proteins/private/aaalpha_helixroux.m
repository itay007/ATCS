function  prop = aaalpha_helixroux(aa) 
%AAALPHA_HELIXROUX calculates the conformational parameter for alpha helix.
%
% Author(s) :  Deleage G., Roux B.
% Reference :  Protein Engineering 1:289-294(1987).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Alpha helix (Deleage & Roux)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [1.489; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        0.966; ...% C  cysteine
        0.924; ...% D  aspartic acid
        1.504; ...% E  glutamic acid
        1.195; ...% F  phenylalanine
        0.510; ...% G  glycine
        1.003; ...% H  histidine
        1.003; ...% I  isoleucine
        NaN;   ...% J  not used 
        1.172; ...% K  lysine
        1.236; ...% L  leucine
        1.363; ...% M  methionine
        0.772; ...% N  asparagine
        NaN;   ...% O  not used 
        0.492; ...% P  proline
        1.164; ...% Q  glutamine
        1.224; ...% R  arginine
        0.739; ...% S  serine
        0.785; ...% T  threonine
        NaN;   ...% U  not used 
        0.990; ...% V  valine
        1.090; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        0.787; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
