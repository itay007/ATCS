function out = mtimes(a,b)
%MTIMES  Overload matrix multiply for DataMatrix objects
%
%   A = MTIMES(DM1, DM2) returns the matrix product A of DataMatrix objects
%   DM1 and DM2. The number of columns of DM1 must equal the number of rows
%   of DM2.
% 
%   A = MTIMES(DM, B) or A = MTIMES(B, DM) returns the matrix product A of
%   a DataMatrix DM and a numeric or a logical array B. The number of
%   columns of DM or B must equal the number of rows of B or DM.
% 
%   D = MTIMES(A,B) is called for the syntax 'A * B' when A or B is a
%   DataMatrix object.
%
%   See also  MTIMES, DATAMATRIX/TIMES.

%   Copyright 2008-2012 The MathWorks, Inc. 


if isa(a, 'bioma.data.DataMatrix') && isa(b, 'bioma.data.DataMatrix')
    out = doMTimes(a.Matrix, b.Matrix);
    warning(message('bioinfo:DataMatrix:mtimes:OutputNotADataMatrix'));
elseif isa(a, 'bioma.data.DataMatrix') && (isnumeric(b) || islogical(b))
    if isscalar(b)
        out = a;
        out.Matrix = doMTimes(a.Matrix, b);
    else
        out = doMTimes(a.Matrix, b);
        warning(message('bioinfo:DataMatrix:mtimes:OutputNotADataMatrix'));
    end
elseif (isnumeric(a) || islogical(a)) && isa(b, 'bioma.data.DataMatrix') 
    if isscalar(a)
        out = b;
        out.Matrix = doMTimes(a, b.Matrix);
    else
        out = doMTimes(a, b.Matrix);
        warning(message('bioinfo:DataMatrix:mtimes:OutputNotADataMatrix'));
    end
else 
     error(message('bioinfo:DataMatrix:mtimes:InvalidType'))
end
end %DataMatrix/mtimes

function om = doMTimes(p1, p2)
try
    om = p1 * p2;
catch ME
   bioinfoprivate.bioclsrethrow('DataMatrix','mtimes', ME)  
end
end % doMTimes
