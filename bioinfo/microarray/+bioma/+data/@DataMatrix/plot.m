function varargout = plot(a,varargin)
%PLOT  Overload linear plot for DataMatrix objects.
%
%   PLOT(Y) plots the columns of DataMatrix object Y versus their index.
%
%   PLOT(X,Y) plots the DataMatrix vector X versus DataMatrix vector Y.  If
%   X or Y is a DataMatrix matrix, then the DataMatrix vector is plotted
%   versus the rows or columns of the matrix, whichever line up. If X is a
%   scalar and Y is a vector, disconnected line objects are created and
%   plotted as discrete points vertically at X.
% 
%   PLOT(X,Y,S) plots various line types, plot symbols and colors specified
%   by a character string S. See the help for MATLAB PLOT function for more
%   details on line types, plot symbols and colors. 
%
%   See also PLOT, DATAMATRIX.

%   Copyright 2008-2012 The MathWorks, Inc.


%== Find all the DataMatrix objects
if isa(a, 'bioma.data.DataMatrix')
    a = a.Matrix;
end 
for i = 1:nargin -1
    if isa(varargin{i}, 'bioma.data.DataMatrix')
        varargin{i} = varargin{i}.Matrix;
    end
end

try
    varargout{1:nargout} = plot(a, varargin{:});
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','plot', ME)
end
end %DataMatrix/plot
    
    
    
    
