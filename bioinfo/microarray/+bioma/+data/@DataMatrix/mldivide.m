function out = mldivide(a,b)
%MLDIVIDE  Overload backslash or left matrix divide for DataMatrix objects
%
%    A = MLDIVIDE(DM1, DM2) performs the matrix division of DataMatrix
%    object DM1 into DataMatrix object DM2.  The output is the numeric
%    results of the left matrix division of DM1 and DM2, and not a
%    DataMatrix object. 
% 
%    A = MLDIVIDE(DM, B) where B is a numeric or a logical array.
%
%    A = MLDIVIDE(B, DM) where B is a numeric or a logical array. If B is a
%    scalar, A is a DataMatrix object with the same row names and column
%    names as DM.
% 
%    D = MLDIVIDE(A,B) is called for the syntax 'A \ B' when A or B is a
%    DataMatrix object.
%
%   See also MLDIVIDE, DATAMATRIX/LDIVIDE, DATAMATRIX/MRDIVIDE.

%   Copyright 2008-2012 The MathWorks, Inc. 


if isa(a, 'bioma.data.DataMatrix') && isa(b, 'bioma.data.DataMatrix')
    out = doMLDivide(a.Matrix, b.Matrix);
    warning(message('bioinfo:DataMatrix:mldivide:OutputNotADataMatrix'));

elseif isa(a, 'bioma.data.DataMatrix') && (isnumeric(b) || islogical(b))
    out = doMLDivide(a.Matrix, b);
    warning(message('bioinfo:DataMatrix:mldivide:OutputNotADataMatrix'));
elseif isa(b, 'bioma.data.DataMatrix') && (isnumeric(a) || islogical(a))
    if isscalar(a)
        out = b;
        out.Matrix = doMLDivide(a, b.Matrix);
    else
        out = doMLDivide(a, b.Matrix);
        warning(message('bioinfo:DataMatrix:mldivide:OutputNotADataMatrix'));
    end
    
else 
     error(message('bioinfo:DataMatrix:mldivide:InvalidType'))
end
end %DataMatrix/mldivide

function om = doMLDivide(p1, p2)
try
    om = p1 \ p2;
catch ME
   bioinfoprivate.bioclsrethrow('DataMatrix','mldivide', ME)  
end
end % doMLDivide
