function b = exp(a) 
%EXP  Overload exponential for DataMatrix object.
% 
%   B = EXP(A) returns the exponential of the elements of the DataMatrix
%   object A. Output B is a DataMatrix object. The size, row names and
%   column names of B are the same as A.
%
%   See also LOG10 DATAMATRIX/LOG, DATAMATRIX/LOG2

%   Copyright 2008-2012 The MathWorks, Inc. 

try
    b = a;
    b.Matrix = exp(a.Matrix);
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','exp', ME);
end

end %DataMatrix/exp
