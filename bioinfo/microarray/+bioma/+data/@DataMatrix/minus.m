function out = minus(a,b)
%MINUS Overloaded array subtraction for DataMatrix objects.
%
%   D = MINUS(DM1, DM2) subtracts DataMatrix object DM2 from DM1. DM1 and
%   DM2 must have the same dimensions, unless one is a scalar (a 1-by-1
%   DataMatrix). Output D is a DataMatrix object. The size, row names and
%   column names of D are the same as DM1, unless DM1 is a scalar; then
%   they are the same as DM2.
% 
%   D = MINUS(DM, A) or MINUS(A, DM) subtracts a numeric or a logical array
%   A from DataMatrix object DM or subtracts DM from A. The size of DM and
%   the size of A must be the same, unless A is a scalar. Output D is a
%   DataMatrix object. The size, row names and column names of D are the
%   same as DM.
% 
%   Note: Arithmetic operations between a scalar DataMatrix object and a
%   non-scalar array is not supported.
% 
%   D = MINUS(A,B) is called for the syntax 'A - B' when A or B is a
%   DataMatrix object.
%
%   See also @DATAMATRIX/PLUS.

%   Copyright 2008-2012 The MathWorks, Inc. 


if isa(a, 'bioma.data.DataMatrix') && isa(b, 'bioma.data.DataMatrix')
    if isscalar(a) && ~isscalar(b)
        out = b;
    else
        out = a;
    end
    out.Matrix = doMinus(a.Matrix, b.Matrix);
elseif isa(a, 'bioma.data.DataMatrix') && (isnumeric(b) || islogical(b))
    if isscalar(a) && ~isscalar(b)
        error(message('bioinfo:DataMatrix:minus:MixScalarDMWithNonScalarArray'))
    end
    out = a;
    out.Matrix = doMinus(a.Matrix, b);
elseif isa(b, 'bioma.data.DataMatrix') && (isnumeric(a) || islogical(a))
    if isscalar(b) && ~isscalar(a)
        error(message('bioinfo:DataMatrix:minus:MixScalarDMWithNonScalarArray'))
    end
    out = b;
    out.Matrix = doMinus(a, b.Matrix);
else 
     error(message('bioinfo:DataMatrix:minus:InvalidType'))
end
end %DataMatrix/minus

function om = doMinus(p1, p2)
try
    om = p1-p2;
catch ME
   bioinfoprivate.bioclsrethrow('DataMatrix','minus', ME)  
end
end % doplus
    