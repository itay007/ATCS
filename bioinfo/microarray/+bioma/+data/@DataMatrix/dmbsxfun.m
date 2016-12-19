function c = dmbsxfun(fun, a, b)
%DMBSXFUN Binary singleton expansion function for DataMatrix objects.
% 
%   C = DMBSXFUN(FUNC, A, B) applies an element-by-element binary operation
%   to DataMatrix objects A and B, with singleton expansion enabled. FUNC
%   is a function handle. FUNC can either be a function handle for a Matlab
%   function, or one of the built-in functions specified in BSXFUN help.
%   See the help for BSXFUN for more details about FUNC.
%
%   A and B can be DataMatrix objects or MATLAB arrays, but one of them
%   must be DataMatrix class. Each dimension of A and B must be equal to
%   each other, or equal to one. The output C is a DataMatrix object with
%   the size equal to max(size(A), size(B)). * (size(A)>0 & size(B)>0). For
%   example, if size(A)==[2 5] and size(B)==[2 1], then size(C)==[2 5]. If
%   one of A or B is a DataMatrix object then size(C) is equal to the size
%   of that DataMatrix; If both A and B are DataMatrix objects, the
%   dimensional properties of C will be the same as max(size(A), size(B)).
%   If A and B are the same size, C contains the same dimensional
%   properties as in A.
% 
%   Note: Singleton expansion of a DataMatrix object to a non-DataMatrix
%   array returns a non-DataMatrix array. 
% 
%   Example:
% 
%     % Subtract the column means from a DataMatrix object
%     A = bioma.data.DataMatrix(magic(3), 'RowNames', true, 'ColNames',true)
%     A = dmbsxfun(@minus, A, mean(A))
%
%   See also BSXFUN, DATAMATRIX/DMARRAYFUN

%   Copyright 2008-2012 The MathWorks, Inc.


bioinfochecknargin(nargin,3,'DataMatrix:dmbsxfun')

%== Initialize
aIsDM = false;
bIsDM = false;
if isa(a, 'bioma.data.DataMatrix')
    aIsDM = true;
end
if isa(b, 'bioma.data.DataMatrix')
    bIsDM = true;
end

sizeA = size(a);
sizeB = size(b);

try
    if aIsDM && bIsDM
        if sizeA == sizeB
            rowNames = a.RowNames;
            colNames = a.ColNames;
        else
            rowNames = getLargerDimNames(a.RowNames, b.RowNames);
            colNames = getLargerDimNames(a.ColNames, b.ColNames);
        end
        cMatrix = bsxfun(fun, a.Matrix, b.Matrix);
    elseif aIsDM && ~bIsDM
        if any(sizeA == 1)&& ~isequal(sizeA,sizeB)
            c = bsxfun(fun, a.Matrix, b);
            throwOutputWarning()
            return;
        end
        
        rowNames = a.RowNames;
        colNames = a.ColNames;
        cMatrix = bsxfun(fun, a.Matrix, b);
    elseif ~aIsDM && bIsDM 
        if any(sizeB == 1)&& ~isequal(sizeA,sizeB)
            c = bsxfun(fun, a, b.Matrix);
            throwOutputWarning()
            return;
        end
        
        rowNames = b.RowNames;
        colNames = b.ColNames;
        cMatrix = bsxfun(fun, a, b.Matrix);
    end
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','dmbsxfun', ME)        
end

c = bioma.data.DataMatrix(cMatrix, rowNames, colNames);
end %DataMatrix/dmbsxfun

%--------------------
function names = getLargerDimNames(nameA, nameB)
if numel(nameA) > numel(nameB)
    names = nameA;
else
    names = nameB;
end
end %getLargerDimNames

function throwOutputWarning()
%== If one of the input is not a DataMatrix and dimension not equal to one
% the singleton expansion of the DataMatrix throw warning
warning(message('bioinfo:DataMatrix:dmbsxfun:DMDimensionIsOne'));
end %throwOutputWarning


