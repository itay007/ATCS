function b = log2(a) 
%LOG2  Overload base 2 logarithm for DataMatrix object.
% 
%   B = LOG2(A) returns the base 2 logarithm of the elements of the
%   DataMatrix object A. Output B is a DataMatrix object. The size, row
%   names and column names of B are the same as A.
%
%   See also DATAMATRIX/LOG, DATAMATRIX/LOG10, DATAMATRIX/POW2.

%   Copyright 2008 The MathWorks, Inc. 

try
    b = a;
    b.Matrix = log2(a.Matrix);
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix', 'log2', ME);
end

end %DataMatrix/log2
