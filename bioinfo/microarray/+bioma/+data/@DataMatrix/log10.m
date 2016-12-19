function b = log10(a) 
%LOG10  Overload common (base 10) logarithm for DataMatrix object.
% 
%   B = LOG10(A) returns the base 10 logarithm of the elements of the
%   DataMatrix object A. Output B is a DataMatrix object. The size, row
%   names and column names of B are the same as A.
%
%   See also LOG10 DATAMATRIX/LOG, DATAMATRIX/LOG2

%   Copyright 2008-2012 The MathWorks, Inc. 

try
    b = a;
    b.Matrix = log10(a.Matrix);
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','log10', ME);
end

end %DataMatrix/log10
