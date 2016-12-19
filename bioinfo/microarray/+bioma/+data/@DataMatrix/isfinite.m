function t = isfinite(a)
%ISFINITE Overload true for finite elements for DataMatrix objects.
% 
%   T = ISFINITE(A) returns a logical array T of the same size of
%   DataMatrix objects A, containing logical 1 where the elements of A are
%   finite and 0 where they are not.
%
%   See also DATAMATRIX/ISINF, DATAMATRIX/ISNAN.

%   Copyright 2008-2012 The MathWorks, Inc. 


try
    t = isfinite(a.Matrix);
catch ME
   bioinfoprivate.bioclsrethrow('DataMatrix','isfinite', ME); 
end
end % DataMatrix/isfinite
