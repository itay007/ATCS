function t = isinf(a)
%ISINF Overload true for infinite elements for DataMatrix objects.
% 
%   T = ISINF(A) returns a logical array T of the same size of DataMatrix
%   objects A, containing logical 1 where the elements of A are +Inf or
%   -Inf and 0 where they are not.
%
%   See also DATAMATRIX/ISFINITE, DATAMATRIX/ISNAN.

%   Copyright 2008-2012 The MathWorks, Inc. 


try
    t = isinf(a.Matrix);
catch ME
   bioinfoprivate.bioclsrethrow('DataMatrix','isinf', ME); 
end
end % DataMatrix/isinf
