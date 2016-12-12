function b = round(a) 
%ROUND  Overloaded round towards nearest integer for DataMatrix object.
% 
%   B = ROUND(A) returns a DataMatrix object B containing the elements of
%   the DataMatrix A rounded to the nearest integers.
%
%   See also ROUND, DATAMATRIX/CEIL, DATAMATRIX/FIX, DATAMATRIX/FLOOR.

%   Copyright 2008-2012 The MathWorks, Inc. 

try
    b = a;
    b.Matrix = round(a.Matrix);
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','round', ME);
end

end %DataMatrix/round
