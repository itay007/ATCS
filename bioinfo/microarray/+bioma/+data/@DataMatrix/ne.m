function t = ne(a,b)
%NE  Overloaded not equal for DataMatrix objects.
% 
%   T = NE(DM1,DM2) performs element-by-element numerical inequality
%   comparisons between DataMatrix objects DM1 and DM2 and returns a
%   logical matrix T of the same size with elements set to logical 1 where
%   the relationship is true and elements set to logical 0 where it is not.
%   DM1 and DM2 must have the same dimensions, unless one is a scalar (a
%   1-by-1 DataMatrix).
% 
%   T = NE(DM, B) or T = NE(B, DM) performs element-by-element numerical
%   inequality comparisons between DataMatrix object DM and a numeric or a
%   logical array B. DM and B must have the same size, unless one is a
%   scalar.
% 
%   T = NE(A,B) is called for the syntax 'A ~= B' when A or B is a
%   DataMatrix object.
%
%   See also DATAMATRIX/EQ

%   Copyright 2008-2012 The MathWorks, Inc. 


if isa(a, 'bioma.data.DataMatrix') && isa(b, 'bioma.data.DataMatrix')
    t = doNE(a.Matrix,b.Matrix);
elseif isa(a, 'bioma.data.DataMatrix') && (isnumeric(b) || islogical(b))
    t = doNE(a.Matrix, b);
elseif isa(b, 'bioma.data.DataMatrix') && (isnumeric(a) || islogical(a))
    t = doNE(a, b.Matrix);
else 
     error(message('bioinfo:DataMatrix:ne:InvalidType'))
end
end %DataMatrix/ne

function om = doNE(p1, p2)
try
    om = ne(p1,p2);
catch ME
   bioinfoprivate.bioclsrethrow('DataMatrix','ne', ME)  
end
end % doNE
