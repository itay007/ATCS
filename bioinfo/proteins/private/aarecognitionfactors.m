function  prop = aarecognitionfactors(aa) 
%AARECOGNITIONFACTORS calculates recognition factors.
%
% Author(s) :  Fraga S.
% Reference :  Can. J. Chem. 60:2606-2610(1982).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Recognition factors';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [78.000;  ...% A  alanine
        NaN;     ...% B  aspartic acid or asparagine
        89.000;  ...% C  cysteine
        81.000;  ...% D  aspartic acid
        78.000;  ...% E  glutamic acid
        81.000;  ...% F  phenylalanine
        84.000;  ...% G  glycine
        84.000;  ...% H  histidine
        88.000;  ...% I  isoleucine
        NaN;     ...% J  not used 
        87.000;  ...% K  lysine
        85.000;  ...% L  leucine
        80.000;  ...% M  methionine
        94.000;  ...% N  asparagine
        NaN;     ...% O  not used 
        91.000;  ...% P  proline
        87.000;  ...% Q  glutamine
        95.000;  ...% R  arginine
        107.000; ...% S  serine
        93.000;  ...% T  threonine
        NaN;     ...% U  not used 
        89.000;  ...% V  valine
        104.000; ...% W  tryptophan
        NaN;     ...% X  any amino acid
        84.000;  ...% Y  tyrosine
        NaN;     ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
