function b = pow2(a) 
%POW2  Overload base 2 power for DataMatrix object.
% 
%   B = POW2(A) for each elements of B is 2 raised to the power of
%   DataMatrix object A. Output B is a DataMatrix object. The size, row
%   names and column names of B are the same as A.
%
%   See also DATAMATRIX/LOG2

%   Copyright 2008 The MathWorks, Inc. 

try
    b = a;
    b.Matrix = pow2(a.Matrix);
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','pow2', ME);
end

end %DataMatrix/pow2
