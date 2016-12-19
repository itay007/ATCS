function  prop = aahphob_black(aa) 
%AAHPHOB_BLACK calculates hydrophobicity of physiological L-alpha amino acids
%
% Author(s) :  Black S.D., Mould D.R.
% Reference :  Anal. Biochem. 193:72-82(1991). <A HREF="http://psyche.uthct.edu/shaun/SBlack/aagrease.html">http://psyche.uthct.edu/shaun/SBlack/aagrease.html</A>.


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Black & Mould)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.616; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        0.680; ...% C  cysteine
        0.028; ...% D  aspartic acid
        0.043; ...% E  glutamic acid
        1.000; ...% F  phenylalanine
        0.501; ...% G  glycine
        0.165; ...% H  histidine
        0.943; ...% I  isoleucine
        NaN;   ...% J  not used 
        0.283; ...% K  lysine
        0.943; ...% L  leucine
        0.738; ...% M  methionine
        0.236; ...% N  asparagine
        NaN;   ...% O  not used 
        0.711; ...% P  proline
        0.251; ...% Q  glutamine
        0.000; ...% R  arginine
        0.359; ...% S  serine
        0.450; ...% T  threonine
        NaN;   ...% U  not used 
        0.825; ...% V  valine
        0.878; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        0.880; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
