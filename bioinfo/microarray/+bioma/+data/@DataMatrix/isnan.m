function t = isnan(a)
%ISNAN Overload true for Not-a-Number for DataMatrix objects.
% 
%   T = ISNAN(A) returns a logical array T of the same size of DataMatrix
%   objects A, containing logical 1 where the elements of A are NaN's and 0
%   where they are not.
%
%   See also DATAMATRIX/ISFINITE, DATAMATRIX/ISINF

%   Copyright 2008-2012 The MathWorks, Inc. 


try
    t = isnan(a.Matrix);
catch ME
   bioinfoprivate.bioclsrethrow('DataMatrix','isnan', ME); 
end
end % DataMatrix/isnan
