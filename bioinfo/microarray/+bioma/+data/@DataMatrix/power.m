function out = power(a,b)
%POWER  Overloaded array power for DataMatrix objects
%
%   D = POWER(DM1, DM2) powers DataMatrix objects DM1 and DM2
%   element-by-element. DM1 and DM2 must have the same dimensions unless
%   one is a scalar (a 1-by-1 DataMatrix). Output D is a DataMatrix object.
%   The size, row names and column names of D are the same as DM1 unless
%   DM1 is a scalar; then they are the same as DM2.
% 
%   D = POWER(DM, A) or D = POWER(A, DM) powers DataMatrix object DM and a
%   numeric or a logical array A. The size of DM and the size of A must be
%   the same unless A is a scalar. Output D is a DataMatrix object. The
%   size, row names and column names of D are the same as DM. 
% 
%   Note: Arithmetic operations between a scalar DataMatrix object and a
%   non-scalar array is not supported.
% 
%   D = POWER(A,B) is called for the syntax 'A .^ B' when A or B is a
%   DataMatrix object.
%
%   See also DATAMATRIX/TIMES.

%   Copyright 2008-2012 The MathWorks, Inc. 


if isa(a, 'bioma.data.DataMatrix') && isa(b, 'bioma.data.DataMatrix')
    if isscalar(a) && ~isscalar(b)
        out = b;
    else
        out = a;
    end
    out.Matrix = doPower(a.Matrix, b.Matrix);
elseif isa(a, 'bioma.data.DataMatrix') && (isnumeric(b) || islogical(b))
    if isscalar(a) && ~isscalar(b)
        error(message('bioinfo:DataMatrix:power:MixScalarDMWithNonScalarArray'))
    end
    out = a;
    out.Matrix = doPower(a.Matrix, b);
elseif isa(b, 'bioma.data.DataMatrix') && (isnumeric(a) || islogical(a))
    if isscalar(b) && ~isscalar(a)
        error(message('bioinfo:DataMatrix:power:MixScalarDMWithNonScalarArray'))
    end
    out = b;
    out.Matrix = doPower(a, b.Matrix);
else 
     error(message('bioinfo:DataMatrix:power:InvalidType'))
end
end %DataMatrix/power

function om = doPower(p1, p2)
try
    om = p1 .^ p2;
catch ME
   bioinfoprivate.bioclsrethrow('DataMatrix','power', ME)  
end
end % doPower
