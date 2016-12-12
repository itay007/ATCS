function  prop = aabulkiness(aa) 
%AABULKINESS calculates bulkiness.
%
% Author(s) :  Zimmerman J.M., Eliezer N., Simha R.
% Reference :  J. Theor. Biol. 21:170-201(1968).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Bulkiness';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [11.500; ...% A  alanine
        NaN;    ...% B  aspartic acid or asparagine
        13.460; ...% C  cysteine
        11.680; ...% D  aspartic acid
        13.570; ...% E  glutamic acid
        19.800; ...% F  phenylalanine
        3.400;  ...% G  glycine
        13.690; ...% H  histidine
        21.400; ...% I  isoleucine
        NaN;    ...% J  not used 
        15.710; ...% K  lysine
        21.400; ...% L  leucine
        16.250; ...% M  methionine
        12.820; ...% N  asparagine
        NaN;    ...% O  not used 
        17.430; ...% P  proline
        14.450; ...% Q  glutamine
        14.280; ...% R  arginine
        9.470;  ...% S  serine
        15.770; ...% T  threonine
        NaN;    ...% U  not used 
        21.570; ...% V  valine
        21.670; ...% W  tryptophan
        NaN;    ...% X  any amino acid
        18.030; ...% Y  tyrosine
        NaN;    ...% Z  glutamic acid or glutamine
];

    prop = data(ndx);
