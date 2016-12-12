function  prop = aamolecularweight(aa) 
%AAMOLECULARWEIGHT calculates the molecular weight of each amino acid.


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Molecular weight';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [89.000;  ...% A	alanine
        NaN;     ...% B  aspartic acid or asparagine
        121.000; ...% C	cysteine
        133.000; ...% D	aspartic acid
        147.000; ...% E	glutamic acid
        165.000; ...% F	phenylalanine
        75.000;  ...% G	glycine
        155.000; ...% H	histidine
        131.000; ...% I	isoleucine
        NaN;     ...% J  not used 
        146.000; ...% K	lysine
        131.000; ...% L	leucine
        149.000; ...% M	methionine
        132.000; ...% N	asparagine
        NaN;     ...% O  not used 
        115.000; ...% P	proline
        146.000; ...% Q	glutamine
        174.000; ...% R	arginine
        105.000; ...% S	serine
        119.000; ...% T	threonine
        NaN;     ...% U not used 
        117.000; ...% V	valine
        204.000; ...% W	tryptophan
        NaN;     ...% X  any amino acid
        181.000; ...% Y	tyrosine
        NaN;     ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
