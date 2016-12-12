function [matrix, matrixInfo] = pam180
%PAM180 substitution matrix in 1/3 bit units,
%   Expected score = -1.51, Entropy = 0.591 bits
%   Lowest score = -10, Highest score = 18
%
%   Order:
%   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
%
%   [MATRIX,MATRIXINFO] = PAM180 returns a structure of information about
%   the matrix with fields Name, Scale, Entropy, ExpectedScore, LowestScore,
%   HighestScore and Order.

% Source:  ftp://ftp.ncbi.nih.gov/blast/matrices/


matrix = [...
    3  -3   0   0  -3  -1   0   1  -2  -1  -3  -2  -2  -5   1   1   2  -8  -5   0   0   0  -1 -10;...
   -3   8  -1  -3  -5   1  -2  -4   2  -3  -4   4  -1  -6  -1  -1  -2   2  -6  -4  -2   0  -2 -10;...
    0  -1   4   3  -5   0   2   0   2  -3  -4   1  -3  -5  -1   1   0  -5  -2  -3   3   1  -1 -10;...
    0  -3   3   5  -7   2   4   0   0  -3  -6   0  -4  -8  -2   0  -1  -9  -6  -3   4   3  -1 -10;...
   -3  -5  -5  -7  13  -7  -7  -5  -4  -3  -8  -7  -7  -6  -4   0  -3 -10   0  -3  -6  -7  -4 -10;...
   -1   1   0   2  -7   6   3  -2   4  -3  -2   0  -1  -6   0  -1  -2  -6  -6  -3   1   5  -1 -10;...
    0  -2   2   4  -7   3   5   0   0  -3  -5  -1  -3  -7  -1  -1  -1  -9  -6  -3   3   5  -1 -10;...
    1  -4   0   0  -5  -2   0   6  -3  -4  -6  -3  -4  -6  -1   1  -1  -9  -7  -2   0  -1  -2 -10;...
   -2   2   2   0  -4   4   0  -3   8  -4  -3  -1  -3  -3  -1  -2  -2  -4   0  -3   1   2  -1 -10;...
   -1  -3  -3  -3  -3  -3  -3  -4  -4   6   2  -3   2   1  -3  -2   0  -7  -2   5  -3  -3  -1 -10;...
   -3  -4  -4  -6  -8  -2  -5  -6  -3   2   7  -4   4   1  -4  -4  -3  -3  -2   2  -5  -3  -2 -10;...
   -2   4   1   0  -7   0  -1  -3  -1  -3  -4   6   1  -7  -2  -1   0  -5  -6  -4   0   0  -1 -10;...
   -2  -1  -3  -4  -7  -1  -3  -4  -3   2   4   1   9   0  -3  -2  -1  -6  -4   2  -3  -2  -1 -10;...
   -5  -6  -5  -8  -6  -6  -7  -6  -3   1   1  -7   0  10  -6  -4  -4   0   7  -2  -6  -7  -3 -10;...
    1  -1  -1  -2  -4   0  -1  -1  -1  -3  -4  -2  -3  -6   8   1   0  -7  -7  -2  -2  -1  -1 -10;...
    1  -1   1   0   0  -1  -1   1  -2  -2  -4  -1  -2  -4   1   3   2  -3  -4  -2   1  -1   0 -10;...
    2  -2   0  -1  -3  -2  -1  -1  -2   0  -3   0  -1  -4   0   2   4  -7  -4   0   0  -1  -1 -10;...
   -8   2  -5  -9 -10  -6  -9  -9  -4  -7  -3  -5  -6   0  -7  -3  -7  18  -1  -8  -7  -8  -6 -10;...
   -5  -6  -2  -6   0  -6  -6  -7   0  -2  -2  -6  -4   7  -7  -4  -4  -1  11  -4  -4  -6  -3 -10;...
    0  -4  -3  -3  -3  -3  -3  -2  -3   5   2  -4   2  -2  -2  -2   0  -8  -4   6  -3  -3  -1 -10;...
    0  -2   3   4  -6   1   3   0   1  -3  -5   0  -3  -6  -2   1   0  -7  -4  -3   4   3  -1 -10;...
    0   0   1   3  -7   5   5  -1   2  -3  -3   0  -2  -7  -1  -1  -1  -8  -6  -3   3   5  -1 -10;...
   -1  -2  -1  -1  -4  -1  -1  -2  -1  -1  -2  -1  -1  -3  -1   0  -1  -6  -3  -1  -1  -1  -1 -10;...
  -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10   1;...
    ];

if nargout >1
    matrixInfo.Name = 'PAM180';
    matrixInfo.Scale = 1/3;
    matrixInfo.Entropy = 0.591 ;
    matrixInfo.ExpectedScore = -1.51;
    matrixInfo.LowestScore = -10;
    matrixInfo.HighestScore = 18;
    matrixInfo.Order = 'ARNDCQEGHILKMFPSTWYVBZX*';
end

