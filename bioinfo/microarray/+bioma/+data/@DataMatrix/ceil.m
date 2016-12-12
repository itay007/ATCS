function b = ceil(a) 
%CEIL  Overloaded round towards plus infinity for DataMatrix object.
% 
%   B = CEIL(A) returns a DataMatrix object B containing the elements of the
%   DataMatrix A rounded to the nearest integers towards infinity.
%
%   See also CEIL, DATAMATRIX/FLOOR, DATAMATRIX/ROUND, DATAMATRIX/FIX.

%   Copyright 2008-2012 The MathWorks, Inc. 

try
    b = a;
    b.Matrix = ceil(a.Matrix);
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','ceil', ME);
end

end %DataMatrix/ceil
