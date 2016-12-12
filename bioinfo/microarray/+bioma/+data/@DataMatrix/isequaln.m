function t = isequaln(a,varargin)
%ISEQUALN True if DataMatrix objects are equal.
% 
%   TF = ISEQUALN(A,B) returns true (1) if the DataMatrix objects A and B
%   are numerically equal, have the same size, and have the same row names
%   and column names, and false (0) otherwise. NaNs are considered equal to
%   each other. The Name property of A and B can be different. 
% 
%   TF = ISEQUALN(A, B, C,...) is true (1) if all the input DataMatrix are
%   numerically equal, have the same size, and have the same row names and
%   column names, and false (0) otherwise. NaNs are considered equal to
%   each other. The input DataMatrix objects do not have to have the same
%   Name property.
%
%   See also ISEQUALN, DATAMATRIX/ISEQUAL

%   Copyright 2012 The MathWorks, Inc. 


%== Input check
bioinfochecknargin(nargin,2,mfilename)

t = isa(a, 'bioma.data.DataMatrix');
if ~t
    return;
end

for i=1:nargin-1
    b = varargin{i};
    if isa(b, 'bioma.data.DataMatrix')
        t = isequaln(a.Matrix, b.Matrix) && ...
            isequal(a.NRows, b.NRows) && ...
            isequal(a.NCols, b.NCols) && ...
            isequal(a.RowNames, b.RowNames) && ...
            isequal(a.ColNames, b.ColNames);
        
    else
        t = false;
    end
    
    if ~t 
        return;
    end
end
end %DataMatrix/isequaln
