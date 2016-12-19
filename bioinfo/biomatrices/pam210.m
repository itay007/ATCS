function [matrix, matrixInfo] = pam210
%PAM210 substitution matrix in 1/3 bit units,
%   Expected score = -1.12, Entropy = 0.470 bits
%   Lowest score = -9, Highest score = 18
%
%   Order:
%   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
%
%   [MATRIX,MATRIXINFO] = PAM210 returns a structure of information about
%   the matrix with fields Name, Scale, Entropy, ExpectedScore, LowestScore,
%   HighestScore and Order.

% Source:  ftp://ftp.ncbi.nih.gov/blast/matrices/


matrix = [...
   2 -2  0  0 -2 -1  0  1 -2 -1 -2 -2 -1 -4  1  1  1 -7 -4  0  0  0  0 -9;...
  -2  7  0 -2 -4  1 -2 -3  2 -2 -4  4 -1 -5  0  0 -1  2 -5 -3 -1  0 -1 -9;...
   0  0  3  2 -4  1  2  0  2 -2 -3  1 -2 -4 -1  1  0 -5 -2 -2  3  1  0 -9;...
   0 -2  2  5 -6  2  4  0  1 -3 -5  0 -3 -7 -1  0  0 -8 -5 -3  4  3 -1 -9;...
  -2 -4 -4 -6 12 -6 -6 -4 -4 -3 -7 -6 -6 -5 -3  0 -3 -9  0 -2 -5 -6 -4 -9;...
  -1  1  1  2 -6  5  3 -2  3 -3 -2  1 -1 -6  0 -1 -1 -6 -5 -2  1  4 -1 -9;...
   0 -2  2  4 -6  3  5  0  0 -2 -4  0 -3 -6 -1  0 -1 -8 -5 -2  3  4 -1 -9;...
   1 -3  0  0 -4 -2  0  5 -3 -3 -5 -2 -3 -5 -1  1  0 -8 -6 -2  0 -1 -1 -9;...
  -2  2  2  1 -4  3  0 -3  8 -3 -2  0 -3 -2  0 -1 -2 -3  0 -3  1  2 -1 -9;...
  -1 -2 -2 -3 -3 -3 -2 -3 -3  5  2 -2  2  1 -3 -2  0 -6 -1  4 -3 -3 -1 -9;...
  -2 -4 -3 -5 -7 -2 -4 -5 -2  2  7 -3  4  2 -3 -3 -2 -2 -1  2 -4 -3 -2 -9;...
  -2  4  1  0 -6  1  0 -2  0 -2 -3  5  1 -6 -2  0  0 -4 -5 -3  0  0 -1 -9;...
  -1 -1 -2 -3 -6 -1 -3 -3 -3  2  4  1  8  0 -3 -2 -1 -5 -3  2 -3 -2 -1 -9;...
  -4 -5 -4 -7 -5 -6 -6 -5 -2  1  2 -6  0 10 -5 -4 -4  0  7 -2 -5 -6 -3 -9;...
   1  0 -1 -1 -3  0 -1 -1  0 -3 -3 -2 -3 -5  7  1  0 -7 -6 -2 -1  0 -1 -9;...
   1  0  1  0  0 -1  0  1 -1 -2 -3  0 -2 -4  1  2  2 -3 -3 -1  1  0  0 -9;...
   1 -1  0  0 -3 -1 -1  0 -2  0 -2  0 -1 -4  0  2  3 -6 -3  0  0 -1  0 -9;...
  -7  2 -5 -8 -9 -6 -8 -8 -3 -6 -2 -4 -5  0 -7 -3 -6 18 -1 -7 -6 -7 -5 -9;...
  -4 -5 -2 -5  0 -5 -5 -6  0 -1 -1 -5 -3  7 -6 -3 -3 -1 11 -3 -4 -5 -3 -9;...
   0 -3 -2 -3 -2 -2 -2 -2 -3  4  2 -3  2 -2 -2 -1  0 -7 -3  5 -2 -2 -1 -9;...
   0 -1  3  4 -5  1  3  0  1 -3 -4  0 -3 -5 -1  1  0 -6 -4 -2  3  2 -1 -9;...
   0  0  1  3 -6  4  4 -1  2 -3 -3  0 -2 -6  0  0 -1 -7 -5 -2  2  4 -1 -9;...
   0 -1  0 -1 -4 -1 -1 -1 -1 -1 -2 -1 -1 -3 -1  0  0 -5 -3 -1 -1 -1 -1 -9;...
  -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9  1;...
    ];

if nargout >1
    matrixInfo.Name = 'PAM210';
    matrixInfo.Scale = 1/3;
    matrixInfo.Entropy = 0.470 ;
    matrixInfo.ExpectedScore = -1.12;
    matrixInfo.LowestScore = -9;
    matrixInfo.HighestScore = 18;
    matrixInfo.Order = 'ARNDCQEGHILKMFPSTWYVBZX*';
end

