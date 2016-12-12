function  prop = aatotalbeta_strand(aa) 
%AATOTALBETA_STRAND calculates conformational preference for total beta strand.
%
% AATOTALBETA_STRAND(AA) calculates conformational preference for total
% beta strand (antiparallel+parallel). 
%
% Author(s) :  Lifson S., Sander C.
% Reference :  Nature 282:109-111(1979).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Total beta strand';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.920; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        1.160; ...% C  cysteine
        0.480; ...% D  aspartic acid
        0.610; ...% E  glutamic acid
        1.250; ...% F  phenylalanine
        0.610; ...% G  glycine
        0.930; ...% H  histidine
        1.810; ...% I  isoleucine
        NaN;   ...% J  not used 
        0.700; ...% K  lysine
        1.300; ...% L  leucine
        1.190; ...% M  methionine
        0.600; ...% N  asparagine
        NaN;   ...% O  not used 
        0.400; ...% P  proline
        0.950; ...% Q  glutamine
        0.930; ...% R  arginine
        0.820; ...% S  serine
        1.120; ...% T  threonine
        NaN;   ...% U  not used 
        1.810; ...% V  valine
        1.540; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        1.530; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
