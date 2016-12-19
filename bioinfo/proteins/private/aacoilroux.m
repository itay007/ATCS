function  prop = aacoilroux(aa) 
%AACOILROUX calculates the conformational parameter for coil.
%
% Author(s) :  Deleage G., Roux B.
% Reference :  Protein Engineering 1:289-294(1987).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Coil (Deleage & Roux)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.824; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        0.953; ...% C  cysteine
        1.197; ...% D  aspartic acid
        0.761; ...% E  glutamic acid
        0.797; ...% F  phenylalanine
        1.251; ...% G  glycine
        1.068; ...% H  histidine
        0.886; ...% I  isoleucine
        NaN;   ...% J  not used 
        0.897; ...% K  lysine
        0.810; ...% L  leucine
        0.810; ...% M  methionine
        1.167; ...% N  asparagine
        NaN;   ...% O  not used 
        1.540; ...% P  proline
        0.947; ...% Q  glutamine
        0.893; ...% R  arginine
        1.130; ...% S  serine
        1.148; ...% T  threonine
        NaN;   ...% U  not used 
        0.772; ...% V  valine
        0.941; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        1.109; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
