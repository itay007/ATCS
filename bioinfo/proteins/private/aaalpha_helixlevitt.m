function  prop = aaalpha_helixlevitt(aa) 
%AAALPHA_HELIXLEVITT calculates the normalized frequency for alpha helix.
%
% Author(s) :  Levitt M.
% Reference :  Biochemistry 17:4277-4285(1978).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Alpha helix (Levitt)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [1.290; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        1.110; ...% C  cysteine
        1.040; ...% D  aspartic acid
        1.440; ...% E  glutamic acid
        1.070; ...% F  phenylalanine
        0.560; ...% G  glycine
        1.220; ...% H  histidine
        0.970; ...% I  isoleucine
        NaN;   ...% J  not used 
        1.230; ...% K  lysine
        1.300; ...% L  leucine
        1.470; ...% M  methionine
        0.900; ...% N  asparagine
        NaN;   ...% O  not used 
        0.520; ...% P  proline
        1.270; ...% Q  glutamine
        0.960; ...% R  arginine
        0.820; ...% S  serine
        0.820; ...% T  threonine
        NaN;   ...% U  not used 
        0.910; ...% V  valine
        0.990; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        0.720; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
