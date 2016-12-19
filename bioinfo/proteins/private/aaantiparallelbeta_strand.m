function  prop = aaantiparallelbeta_strand(aa) 
%AAANTIPARALLELBETA_STRAND calculates the conformational preference for anti-parallel beta strands.
%
%  AAANTIPARALLELBETA_STRAND(AA) returns a vector of values of the 
%  conformational preference for antiparallel beta strands for the amino
%  acids in sequence AA.
%
% Author(s) :  Lifson S., Sander C.
% Reference :  Nature 282:109-111(1979).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Antiparallel beta strand';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.900; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        1.240; ...% C  cysteine
        0.470; ...% D  aspartic acid
        0.620; ...% E  glutamic acid
        1.230; ...% F  phenylalanine
        0.560; ...% G  glycine
        1.120; ...% H  histidine
        1.540; ...% I  isoleucine
        NaN;   ...% J  not used 
        0.740; ...% K  lysine
        1.260; ...% L  leucine
        1.090; ...% M  methionine
        0.620; ...% N  asparagine
        NaN;   ...% O  not used 
        0.420; ...% P  proline
        1.180; ...% Q  glutamine
        1.020; ...% R  arginine
        0.870; ...% S  serine
        1.300; ...% T  threonine
        NaN;   ...% U  not used 
        1.530; ...% V  valine
        1.750; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        1.680; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
