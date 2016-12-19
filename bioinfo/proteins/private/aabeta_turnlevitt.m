function  prop = aabeta_turnlevitt(aa) 
%AABETA_TURNLEVITT calculates normalized frequency for beta-turn.
%
% Author(s) :  Levitt M.
% Reference :  Biochemistry 17:4277-4285(1978).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Beta-turn (Levitt)';    
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.770; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        0.810; ...% C  cysteine
        1.410; ...% D  aspartic acid
        0.990; ...% E  glutamic acid
        0.590; ...% F  phenylalanine
        1.640; ...% G  glycine
        0.680; ...% H  histidine
        0.510; ...% I  isoleucine
        NaN;   ...% J  not used 
        0.960; ...% K  lysine
        0.580; ...% L  leucine
        0.410; ...% M  methionine
        1.280; ...% N  asparagine
        NaN;   ...% O  not used 
        1.910; ...% P  proline
        0.980; ...% Q  glutamine
        0.880; ...% R  arginine
        1.320; ...% S  serine
        1.040; ...% T  threonine
        NaN;   ...% U  not used 
        0.470; ...% V  valine
        0.760; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        1.050; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
