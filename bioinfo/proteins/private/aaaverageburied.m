function  prop = aaaverageburied(aa) 
%AAAVERAGEBURIED calculates average area buried on transfer from standard state to folded protein.
%
% Author(s) :  Rose G.D., Geselowitz A.R., Lesser G.J., Lee R.H., Zehfus M.H.
% Reference :  Science 229:834-838(1985).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_  Average area buried';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [86.600;  ...% A	 alanine
        NaN;     ...% B  aspartic acid or asparagine
        132.300; ...% C	 cysteine
        97.800;  ...% D	 aspartic acid
        113.900; ...% E	 glutamic acid
        194.100; ...% F	 phenylalanine
        62.900;  ...% G	 glycine
        155.800; ...% H	 histidine
        158.000; ...% I	 isoleucine
        NaN;     ...% J  not used 
        115.500; ...% K	 lysine
        164.100; ...% L	 leucine
        172.900; ...% M	 methionine
        103.300; ...% N	 asparagine
        NaN;     ...% O  not used 
        92.900;  ...% P	 proline
        119.200; ...% Q	 glutamine
        162.200; ...% R	 arginine
        85.600;  ...% S	 serine
        106.500; ...% T	 threonine
        NaN;     ...% U  not used 
        141.000; ...% V	 valine
        224.600; ...% W	 tryptophan
        NaN;     ...% X	 any amino acid
        177.700; ...% Y	 tyrosine
        NaN;     ...% Z	 glutamic acid or glutamine
];

prop = data(ndx);
