function  prop = aarefractivity(aa) 
%AAREFRACTIVITY calculates refractivity.
%
% Author(s) :  Jones. D.D.
% Reference :  J. Theor. Biol. 50:167-184(1975).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Refractivity';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [4.340;  ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        35.770; ...% C  cysteine
        13.280; ...% D  aspartic acid
        17.560; ...% E  glutamic acid
        29.400; ...% F  phenylalanine
        0.000;  ...% G  glycine
        21.810; ...% H  histidine
        18.780; ...% I  isoleucine
        NaN;    ...% J  not used 
        21.290; ...% K  lysine
        19.060; ...% L  leucine
        21.640; ...% M  methionine
        12.000; ...% N  asparagine
        NaN;    ...% O  not used 
        10.930; ...% P  proline
        17.260; ...% Q  glutamine
        26.660; ...% R  arginine
        6.350;  ...% S  serine
        11.010; ...% T  threonine
        NaN;    ...% U  not used 
        13.920; ...% V  valine
        42.530; ...% W  tryptophan
        NaN;    ...% X  any amino acid
        31.530; ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
