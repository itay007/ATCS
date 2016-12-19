function varargout = size(a,dim)
%SIZE Size of a DataMatrix object.
% 
%   D = SIZE(A) returns the two-element row vector D = [NROWS,NCOLS]
%   containing the number of rows and columns in the DataMatrix A.
%
%   [NROWS,NCOLS] = SIZE(A) returns the number of rows and columns in the
%   DataMatrix A as separate output variables. 
%
%   [M1,M2,M3,...,MN] = SIZE(A), for N>1, returns the sizes of the first N 
%   dimensions of the DataMatrix A.  If the number of output arguments N
%   does not equal NDIMS(A), then for:
%
%   N > NDIMS(A), SIZE returns ones in the "extra" variables, i.e., outputs
%                 NDIMS(A)+1 through N.
%   N < NDIMS(A), MN contains the product of the sizes of dimensions N
%                 through NDIMS(A).
%
%   M = SIZE(A,DIM) returns the length of the dimension specified by the
%   scalar DIM.  For example, SIZE(A,1) returns the number of rows. If
%   DIM > NDIMS(A), M will be 1.
%
%   See also DATAMATRIX/LENGTH, DATAMATRIX/NDIMS.

%   Copyright 2008-2012 The MathWorks, Inc. 

if nargin == 1
    if nargout < 2
        varargout{1:nargout} = size(a.Matrix);
    elseif nargout == 2        
        [varargout{1}, varargout{2}] = size(a.Matrix);
    else
        [varargout{1}, varargout{2}] = size(a.Matrix);
        varargout(3:nargout) = {1};
    end
        
else
    if nargout > 1
        error(message('bioinfo:DataMatrix:size:TooManyOutputs'));
    end
    
    try
        varargout{1} = size(a.Matrix, dim);
    catch ME
        bioinfoprivate.bioclsrethrow('DataMatrix','size', ME);
    end
end

end % size method
