function validateMatrix(data, varName, clsName)
% VALIDATEMATRIX Check the validity of numeric array.
% 
%   VALIDATEMATRIX(DATA, VARNAME, CLSNAME) Check the input DATA is numeric,
%   two diemesional and real array or a DataMatrix object. Throw the error
%   message with VARNAME and CLSNAME.

% Copyright 2009 The MathWorks, Inc.

validateattributes(data, {'numeric', 'bioma.data.DataMatrix'},...
    {'2d', 'real'}, clsName, varName)
end
