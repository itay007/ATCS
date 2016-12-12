function  prop = aabeta_sheetlevitt(aa) 
%AABETA_SHEETLEVITT calculates the normalized frequency for beta-sheet.
%
% Author(s) :  Levitt M.
% Reference :  Biochemistry 17:4277-4285(1978).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Beta-sheet (Levitt)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [0.900; ...% A  alanine
        NaN;   ...% B  aspartic acid or asparagine
        0.740; ...% C  cysteine
        0.720; ...% D  aspartic acid
        0.750; ...% E  glutamic acid
        1.320; ...% F  phenylalanine
        0.920; ...% G  glycine
        1.080; ...% H  histidine
        1.450; ...% I  isoleucine
        NaN;   ...% J  not used 
        0.770; ...% K  lysine
        1.020; ...% L  leucine
        0.970; ...% M  methionine
        0.760; ...% N  asparagine
        NaN;   ...% O  not used 
        0.640; ...% P  proline
        0.800; ...% Q  glutamine
        0.990; ...% R  arginine
        0.950; ...% S  serine
        1.210; ...% T  threonine
        NaN;   ...% U  not used 
        1.490; ...% V  valine
        1.140; ...% W  tryptophan
        NaN;   ...% X  any amino acid
        1.250; ...% Y  tyrosine
        NaN;   ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
