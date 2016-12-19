function [matrix, matrixInfo] = pam450
%PAM450 substitution matrix in 1/6 bit units,
%   Expected score = -0.476, Entropy = 0.105 bits
%   Lowest score = -9, Highest score = 30
%
%   Order:
%   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
%
%   [MATRIX,MATRIXINFO] = PAM450 returns a structure of information about
%   the matrix with fields Name, Scale, Entropy, ExpectedScore, LowestScore,
%   HighestScore and Order.

% Source:  ftp://ftp.ncbi.nih.gov/blast/matrices/


matrix = [...
   1 -1  0  1 -2  0  1  1 -1  0 -2  0 -1 -3  1  1  1 -6 -4  0  1  0  0 -9;...
  -1  5  1  0 -4  2  0 -1  2 -2 -3  4  0 -4  0  0  0  4 -4 -2  0  1  0 -9;...
   0  1  1  2 -4  1  2  1  1 -1 -2  1 -1 -4  0  1  0 -4 -3 -1  1  1  0 -9;...
   1  0  2  3 -5  2  3  1  1 -2 -3  1 -2 -5  0  1  0 -7 -5 -1  2  2  0 -9;...
  -2 -4 -4 -5 20 -5 -5 -3 -4 -2 -6 -5 -5 -4 -2  0 -2 -9  1 -2 -4 -5 -2 -9;...
   0  2  1  2 -5  3  2  0  3 -1 -2  1 -1 -4  1  0  0 -5 -4 -1  2  2  0 -9;...
   1  0  2  3 -5  2  3  1  1 -2 -3  1 -2 -5  0  1  0 -7 -5 -1  2  3  0 -9;...
   1 -1  1  1 -3  0  1  4 -1 -2 -3 -1 -2 -5  0  1  1 -8 -5 -1  1  0  0 -9;...
  -1  2  1  1 -4  3  1 -1  5 -2 -2  1 -1 -2  0  0 -1 -2  0 -2  1  2  0 -9;...
   0 -2 -1 -2 -2 -1 -2 -2 -2  4  3 -2  3  2 -1 -1  0 -5  0  3 -2 -2  0 -9;...
  -2 -3 -2 -3 -6 -2 -3 -3 -2  3  7 -2  4  4 -2 -2 -1 -1  1  3 -3 -2 -1 -9;...
   0  4  1  1 -5  1  1 -1  1 -2 -2  4  0 -5  0  0  0 -3 -5 -2  1  1  0 -9;...
  -1  0 -1 -2 -5 -1 -2 -2 -1  3  4  0  4  1 -1 -1  0 -4 -1  2 -2 -1  0 -9;...
  -3 -4 -4 -5 -4 -4 -5 -5 -2  2  4 -5  1 13 -5 -3 -3  2 12  0 -4 -5 -2 -9;...
   1  0  0  0 -2  1  0  0  0 -1 -2  0 -1 -5  5  1  1 -6 -5 -1  0  0  0 -9;...
   1  0  1  1  0  0  1  1  0 -1 -2  0 -1 -3  1  1  1 -3 -3 -1  1  0  0 -9;...
   1  0  0  0 -2  0  0  1 -1  0 -1  0  0 -3  1  1  1 -5 -3  0  0  0  0 -9;...
  -6  4 -4 -7 -9 -5 -7 -8 -2 -5 -1 -3 -4  2 -6 -3 -5 30  2 -6 -6 -6 -4 -9;...
  -4 -4 -3 -5  1 -4 -5 -5  0  0  1 -5 -1 12 -5 -3 -3  2 14 -2 -4 -4 -2 -9;...
   0 -2 -1 -1 -2 -1 -1 -1 -2  3  3 -2  2  0 -1 -1  0 -6 -2  3 -1 -1  0 -9;...
   1  0  1  2 -4  2  2  1  1 -2 -3  1 -2 -4  0  1  0 -6 -4 -1  2  2  0 -9;...
   0  1  1  2 -5  2  3  0  2 -2 -2  1 -1 -5  0  0  0 -6 -4 -1  2  2  0 -9;...
   0  0  0  0 -2  0  0  0  0  0 -1  0  0 -2  0  0  0 -4 -2  0  0  0  0 -9;...
  -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9  1;...
    ];

if nargout >1
    matrixInfo.Name = 'PAM450';
    matrixInfo.Scale = 1/6;
    matrixInfo.Entropy = 0.105 ;
    matrixInfo.ExpectedScore = -0.476;
    matrixInfo.LowestScore = -9;
    matrixInfo.HighestScore = 30;
    matrixInfo.Order = 'ARNDCQEGHILKMFPSTWYVBZX*';
end

