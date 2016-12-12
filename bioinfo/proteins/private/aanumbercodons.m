function  prop = aanumbercodons(aa) 
%AANUMBERCODONS calculates the number of codon(s) coding for each residue.
%
% AANUMBERCODONS calculates number of codon(s) coding for each amino acid
% using the universal genetic code. 


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Number of codon(s)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [ 	4.000; ...% A	alanine
NaN;   ...% B  aspartic acid or asparagine
                1.000; ...% C	cysteine
                2.000; ...% D	aspartic acid
                2.000; ...% E	glutamic acid
                2.000; ...% F	phenylalanine
                4.000; ...% G	glycine
                2.000; ...% H	histidine
                3.000; ...% I	isoleucine
                NaN;   ...% J  not used 
                2.000; ...% K	lysine
                6.000; ...% L	leucine
                1.000; ...% M	methionine
                2.000; ...% N	asparagine
                NaN;   ...% O  not used 
                4.000; ...% P	proline
                2.000; ...% Q	glutamine
                6.000; ...% R	arginine
                6.000; ...% S	serine
                4.000; ...% T	threonine
                NaN;   ...% U not used 
                4.000; ...% V	valine
                1.000; ...% W	tryptophan
                NaN;   ...% ANY	Any amino acid
                2.000; ...% Y	tyrosine
                NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
