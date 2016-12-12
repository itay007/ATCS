function out = mrdivide(a,b)
%MRDIVIDE  Overload slash or right matrix divide for DataMatrix objects
%
%    A = MRDIVIDE(DM1, DM2) performs the matrix division of DataMatrix
%    object DM2 into DataMatrix object DM1.  The output is the numeric
%    results of the right matrix division of DM1 and DM2, and not a
%    DataMatrix object. 
% 
%    A = MRDIVIDE(DM, B) where B is a numeric or a logical array. If B is a
%    scalar, A is a DataMatrix object with the same row names and column
%    names as DM.
%
%    A = MRDIVIDE(B, DM) where B is a numeric or a logical array.
% 
%    D = MRDIVIDE(A,B) is called for the syntax 'A / B' when A or B is a
%    DataMatrix object.
%
%   See also MRDIVIDE, DATAMATRIX/MLDIVIDE, DATAMATRIX/RDIVIDE.

%   Copyright 2008-2012 The MathWorks, Inc. 


if isa(a, 'bioma.data.DataMatrix') && isa(b, 'bioma.data.DataMatrix')
    out = doMRDivide(a.Matrix, b.Matrix);
    warning(message('bioinfo:DataMatrix:mrdivide:OutputNotADataMatrix'));
elseif isa(a, 'bioma.data.DataMatrix') && (isnumeric(b) || islogical(b))
    if isscalar(b)
        out = a;
        out.Matrix = doMRDivide(a.Matrix,b);
    else
        out = doMRDivide(a.Matrix, b);
        warning(message('bioinfo:DataMatrix:mrdivide:OutputNotADataMatrix'));
    end
elseif isa(b, 'bioma.data.DataMatrix') && (isnumeric(a) || islogical(a))
    out = doMRDivide(a, b.Matrix);
    warning(message('bioinfo:DataMatrix:mrdivide:OutputNotADataMatrix'));
else
     error(message('bioinfo:DataMatrix:mrdivide:InvalidType'))
end
end %DataMatrix/mrdivide

function om = doMRDivide(p1, p2)
try
    om = p1 / p2;
catch ME
   bioinfoprivate.bioclsrethrow('DataMatrix','mrdivide', ME)  
end
end % doMRDivide
