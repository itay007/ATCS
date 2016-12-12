function  prop = aahphob_argos(aa) 
%AAHPHOB_ARGOS calculates membrane buried helix parameter.
%
% Author(s) :  Rao M.J.K., Argos P.
% Reference :  Biochim. Biophys. Acta 869:197-214(1986).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Rao & Argos)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [1.360; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        1.270; ...% C  cysteine
        0.110; ...% D  aspartic acid
        0.250; ...% E  glutamic acid
        1.570; ...% F  phenylalanine
        1.090; ...% G  glycine
        0.680; ...% H  histidine
        1.440; ...% I  isoleucine
        NaN;   ...% J  not used 
        0.090; ...% K  lysine
        1.470; ...% L  leucine
        1.420; ...% M  methionine
        0.330; ...% N  asparagine
        NaN;   ...% O  not used 
        0.540; ...% P  proline
        0.330; ...% Q  glutamine
        0.150; ...% R  arginine
        0.970; ...% S  serine
        1.080; ...% T  threonine
        NaN;   ...% U  not used 
        1.370; ...% V  valine
        1.000; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        0.830; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
