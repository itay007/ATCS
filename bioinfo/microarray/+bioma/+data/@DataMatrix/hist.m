function varargout = hist(varargin)
%HIST  Overload histogram for DataMatrix objects.
% 
%   N = HIST(Y) bins the elements of the DataMatrix object Y into 10
%   equally spaced containers and returns the number of elements in each
%   container.  If Y is a DataMatrix matrix, HIST works down the columns.
%
%   N = HIST(Y,M), where M is a scalar, uses M bins.
%
%   N = HIST(Y,X), where X is a vector, returns the distribution of Y
%   among bins with centers specified by X. The first bin includes
%   data between -inf and the first center and the last bin
%   includes data between the last bin and inf. Note: Use HISTC if
%   it is more natural to specify bin edges instead. 
%
%   [N,X] = HIST(...) also returns the position of the bin centers in X.
%
%   HIST(...) without output arguments produces a histogram bar plot of
%   the results. The bar edges on the first and last bins may extend to
%   cover the min and max of the data unless a matrix of data is supplied.
%
%   HIST(AX,...) plots into AX instead of GCA.
%
%   See also DATAMATRIX.

%   Copyright 2008-2012 The MathWorks, Inc. 


try 
    if isscalar(varargin{1}) && ishandle(varargin{1})
        if nargout > 0
            varargout{1:nargout} = hist(varargin{1}, varargin{2}.Matrix, varargin{3:end});
        else
             hist(varargin{1}, varargin{2}.Matrix, varargin{3:end});
        end
    else
        if nargout > 0
            varargout{1:nargout} = hist(varargin{1}.Matrix, varargin{2:end});
        else
            hist(varargin{1}.Matrix, varargin{2:end});
        end
    end
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','hist', ME)
end
