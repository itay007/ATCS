function b = floor(a) 
%FLOOR  Overload rounded towards minus infinity for DataMatrix object.
% 
%   B = FLOOR(A) returns a DataMatrix object B containing the elements of
%   the DataMatrix A rounded to the nearest integers towards minus
%   infinity.
%
%   See also FLOOR, DATAMATRIX/CEIL, DATAMATRIX/ROUND, DATAMATRIX/FIX.

%   Copyright 2008 The MathWorks, Inc. 

try
    b = a;
    b.Matrix = floor(a.Matrix);
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','floor', ME);
end

end %DataMatrix/floor
