function a = cat(dim,varargin)
%CAT Concatenate DataMatrix objects.
%   DS = CAT(DIM, DM1, DM2, ...) concatenates the DataMatrix objects DM1,
%   DM2, ... along dimension DIM by calling the @DATAMATRIX/HORZCAT or
%   @DATAMATRIX/VERTCAT method. DIM must be 1 or 2.
%
%   See also DATAMATRIX/HORZCAT, DATAMATRIX/VERTCAT.

%   Copyright 2008-2012 The MathWorks, Inc. 


if dim == 1
    a = vertcat(varargin{:});
elseif dim == 2
    a = horzcat(varargin{:});
else
    error(message('bioinfo:DataMatrix:cat:InvalidDim'));
end

end % DataMatrix/cat
