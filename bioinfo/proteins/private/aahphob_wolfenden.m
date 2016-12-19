function  prop = aahphob_wolfenden(aa) 
%AAHPHOB_WOLFENDEN calculates the hydration potential (kcal/mole) at 25øC.
%
% Author(s) :  Wolfenden R.V., Andersson L., Cullis P.M., Southgate C.C.F.
% Reference :  Biochemistry 20:849-855(1981).


%   Copyright 2003-2004 The MathWorks, Inc.

% no input, return description string 
if nargin == 0, 
  prop = 'AAProp_ Hydrophobicity (Wolfenden et al.)';
  return
end

% create index, and index data 
ndx = double(lower(aa)) - 96;
data = [1.940;   ...% A  alanine
        NaN;     ...% B  aspartic acid or asparagine
        -1.240;  ...% C  cysteine
        -10.950; ...% D  aspartic acid
        -10.200; ...% E  glutamic acid
        -0.760;  ...% F  phenylalanine
        2.390;   ...% G  glycine
        -10.270; ...% H  histidine
        2.150;   ...% I  isoleucine
        NaN;     ...% J  not used 
        -9.520;  ...% K  lysine
        2.280;   ...% L  leucine
        -1.480;  ...% M  methionine
        -9.680;  ...% N  asparagine
        NaN;     ...% O  not used 
        0.000;   ...% P  proline
        -9.380;  ...% Q  glutamine
        -19.920; ...% R  arginine
        -5.060;  ...% S  serine
        -4.880;  ...% T  threonine
        NaN;     ...% U  not used 
        1.990;   ...% V  valine
        -5.880;  ...% W  tryptophan
        NaN;     ...% X  any amino acid
        -6.110;  ...% Y  tyrosine
        NaN;     ...% Z  glutamic acid or glutamine
];

prop = data(ndx);
