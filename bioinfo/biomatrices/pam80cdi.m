function [matrix, matrixInfo] = pam80cdi
%PAM80CDI substitution matrix in 1/10 Bit Units, 
%   Expected score = -12.6, Entropy = 1.44 bits
%   Lowest score = -53, Highest score = 63
%
%   Order:
%   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
%
%   [MATRIX,MATRIXINFO] = PAM80CDI returns a structure of information about
%   the matrix with fields Name, Scale, Entropy, ExpectedScore, LowestScore,
%   HighestScore and Order.

% Source:  ftp://ftp.ncbi.nih.gov/blast/matrices/


matrix = [...
   21 -19  -6  -6 -19 -10  -3   1 -20 -10 -18 -18 -13 -27   1   5   5 -42 -25  -2  -6  -5  -6 -53;...
  -19  37 -12 -25 -26   1 -22 -28   2 -16 -26  10  -9 -31  -8  -6 -16   1 -33 -23 -17  -7 -14 -53;...
   -6 -12  27  13 -31  -5   1  -5   8 -15 -23   2 -22 -27 -14   5  -1 -27 -13 -21  21  -1  -7 -53;...
   -6 -25  13  30 -43   0  18  -5  -6 -22 -37  -9 -30 -45 -20  -6 -10 -49 -34 -23  24  12 -13 -53;...
  -19 -26 -31 -43  46 -43 -43 -29 -24 -19 -47 -44 -42 -39 -24  -5 -22 -52  -8 -17 -36 -43 -27 -53;...
  -10   1  -5   0 -43  34  12 -19  12 -22 -14  -4  -9 -39  -4 -13 -14 -39 -35 -19  -2  25 -10 -53;...
   -3 -22   1  18 -43  12  30  -9  -9 -17 -29 -10 -21 -44 -14  -9 -14 -53 -29 -18  12  24 -11 -53;...
    1 -28  -5  -5 -29 -19  -9  28 -25 -29 -34 -21 -26 -31 -15   0 -12 -49 -42 -15  -5 -12 -13 -53;...
  -20   2   8  -6 -24  12  -9 -25  39 -26 -18 -14 -27 -17  -9 -14 -19 -22  -6 -19   2   4 -12 -53;...
  -10 -16 -15 -22 -19 -22 -17 -29 -26  34   4 -18   5  -1 -24 -18  -3 -42 -16  16 -18 -19 -11 -53;...
  -18 -26 -23 -37 -47 -14 -29 -34 -18   4  30 -25  12  -1 -22 -25 -19 -17 -19   0 -29 -20 -16 -53;...
  -18  10   2  -9 -44  -4 -10 -21 -14 -18 -25  28   0 -43 -18  -8  -6 -33 -31 -26  -3  -7 -13 -53;...
  -13  -9 -22 -30 -42  -9 -21 -26 -27   5  12   0  47  -8 -23 -15  -9 -37 -30   3 -26 -15 -12 -53;...
  -27 -31 -27 -45 -39 -39 -44 -31 -17  -1  -1 -43  -8  41 -33 -21 -26 -10  19 -20 -34 -42 -23 -53;...
    1  -8 -14 -20 -24  -4 -14 -15  -9 -24 -22 -18 -23 -33  35   1  -8 -43 -41 -15 -17  -9 -12 -53;...
    5  -6   5  -6  -5 -13  -9   0 -14 -18 -25  -8 -15 -21   1  22   8 -15 -21 -15   0 -10  -6 -53;...
    5 -16  -1 -10 -22 -14 -14 -12 -19  -3 -19  -6  -9 -26  -8   8  27 -39 -20  -4  -5 -14  -7 -53;...
  -42   1 -27 -49 -52 -39 -53 -49 -22 -42 -17 -33 -37 -10 -43 -15 -39  63 -12 -49 -35 -45 -33 -53;...
  -25 -33 -13 -34  -8 -35 -29 -42  -6 -16 -19 -31 -30  19 -41 -21 -20 -12  45 -23 -20 -31 -22 -53;...
   -2 -23 -21 -23 -17 -19 -18 -15 -19  16   0 -26   3 -20 -15 -15  -4 -49 -23  29 -22 -19 -11 -53;...
   -6 -17  21  24 -36  -2  12  -5   2 -18 -29  -3 -26 -34 -17   0  -5 -35 -20 -22  23   8 -10 -53;...
   -5  -7  -1  12 -43  25  24 -12   4 -19 -20  -7 -15 -42  -9 -10 -14 -45 -31 -19   8  25 -11 -53;...
   -6 -14  -7 -13 -27 -10 -11 -13 -12 -11 -16 -13 -12 -23 -12  -6  -7 -33 -22 -11 -10 -11 -13 -53;...
  -53 -53 -53 -53 -53 -53 -53 -53 -53 -53 -53 -53 -53 -53 -53 -53 -53 -53 -53 -53 -53 -53 -53   1;...
    ];

if nargout >1
    matrixInfo.Name = 'PAM80CDI';
    matrixInfo.Scale = 0.1;
    matrixInfo.Entropy = 1.44 ;
    matrixInfo.ExpectedScore = -12.6;
    matrixInfo.LowestScore = -53;
    matrixInfo.HighestScore = 63;
    matrixInfo.Order = 'ARNDCQEGHILKMFPSTWYVBZX*';
end

