function  prop = aahphob_rose(aa) 
%AAHPHOB_ROSE calculates mean fractional area loss.
%
%   AAHPHOB_ROSE(AA) calculates the mean fractional area loss (f) [average
%   area buried/standard state area]. 
%
% Author(s) :  Rose G.D., Geselowitz A.R., Lesser G.J., Lee R.H., Zehfus M.H.
% Reference :  Science 229:834-838(1985).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Rose & al)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.740; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        0.910; ...% C  cysteine
        0.620; ...% D  aspartic acid
        0.620; ...% E  glutamic acid
        0.880; ...% F  phenylalanine
        0.720; ...% G  glycine
        0.780; ...% H  histidine
        0.880; ...% I  isoleucine
        NaN;   ...% J  not used 
        0.520; ...% K  lysine
        0.850; ...% L  leucine
        0.850; ...% M  methionine
        0.630; ...% N  asparagine
        NaN;   ...% O  not used 
        0.640; ...% P  proline
        0.620; ...% Q  glutamine
        0.640; ...% R  arginine
        0.660; ...% S  serine
        0.700; ...% T  threonine
        NaN;   ...% U  not used 
        0.860; ...% V  valine
        0.850; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        0.760; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
