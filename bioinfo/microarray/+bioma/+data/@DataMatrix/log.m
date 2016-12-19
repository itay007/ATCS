function b = log(a) 
%LOG  Overload natural logarithm for DataMatrix object.
% 
%   B = LOG(A) returns the natural logarithm of the elements of the
%   DataMatrix object A. Output B is a DataMatrix object. The size, row
%   names and column names of B are the same as A.
%
%   See also LOG DATAMATRIX/LOG2, DATAMATRIX/LOG10

%   Copyright 2008-2012 The MathWorks, Inc. 

try
    b = a;
    b.Matrix = log(a.Matrix);
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix', 'log', ME);
end

end %DataMatrix/log
