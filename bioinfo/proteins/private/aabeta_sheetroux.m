function  prop = aabeta_sheetroux(aa) 
%AABETA_SHEETROUX calculates the conformational parameter for beta-sheet.
%
% Author(s) :  Deleage G., Roux B.
% Reference :  Protein Engineering 1:289-294(1987).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Beta-sheet (Roux)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.709; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        1.191; ...% C  cysteine
        0.541; ...% D  aspartic acid
        0.567; ...% E  glutamic acid
        1.393; ...% F  phenylalanine
        0.657; ...% G  glycine
        0.863; ...% H  histidine
        1.799; ...% I  isoleucine
        NaN;   ...% J  not used 
        0.721; ...% K  lysine
        1.261; ...% L  leucine
        1.210; ...% M  methionine
        0.604; ...% N  asparagine
        NaN;   ...% O  not used 
        0.354; ...% P  proline
        0.840; ...% Q  glutamine
        0.920; ...% R  arginine
        0.928; ...% S  serine
        1.221; ...% T  threonine
        NaN;   ...% U  not used 
        1.965; ...% V  valine
        1.306; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        1.266; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
