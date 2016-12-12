function b = fix(a) 
%FIX  Overloaded round towards zero for DataMatrix object.
% 
%   B = FIX(A) returns a DataMatrix object B containing the elements of the
%   DataMatrix A rounded to the nearest integers towards zero.
%
%   See also FIX, DATAMATRIX/CEIL, DATAMATRIX/FLOOR, DATAMATRIX/ROUND.

%   Copyright 2008 The MathWorks, Inc. 

try
    b = a;
    b.Matrix = fix(a.Matrix);
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','fix', ME);
end

end %DataMatrix/fix
