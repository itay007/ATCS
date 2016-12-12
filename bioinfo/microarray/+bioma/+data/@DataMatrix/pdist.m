function Y = pdist(X,dist,varargin)
%PDIST Overload pairwise distance between observations for DataMatrix objects.
% 
%   Y = PDIST(X) returns a vector Y containing the Euclidean distances
%   between each pair of observations in the N-by-P DataMatrix X.  Rows of
%   X correspond to observations, columns correspond to variables.  Y is a
%   1-by-(N*(N-1)/2) row vector, corresponding to the N*(N-1)/2 pairs of
%   observations in X.
%
%   Y = PDIST(X, DISTANCE) computes Y using DISTANCE. See the help for
%   PDIST for more details on choices for the distance. 
% 
%   See also PDIST, DATAMATRIX/PCA.

%   Copyright 2008-2012 The MathWorks, Inc.


try
    Y = pdist(X.Matrix, dist, varargin{:});
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','pdist', ME);
end
end %DataMatrix/pdist
