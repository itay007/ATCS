function b = sqrt(a) 
%SQRT  Overload element-wise square root for DataMatrix object.
% 
%   B = SQRT(A) returns a DataMatrix object B containg the square root of
%   the elements of the DataMatrix object A. The size, row names and column
%   names of B are the same as A.
% 
%   See also DATAMATRIX/POWER.

%   Copyright 2008-2012 The MathWorks, Inc. 

try
    b = a;
    b.Matrix = sqrt(a.Matrix);
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix', 'sqrt', ME);
end

end %DataMatrix/sqrt
