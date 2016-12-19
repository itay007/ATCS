function [idx, varargout] = kmeans(X, k, varargin)
%KMEANS Overload K-mean clustering for DataMatrix objects.
% 
%   IDX = KMEANS(X, K) partitions the points in the N-by-P DataMatrix
%   object X into K clusters.  This partition minimizes the sum, over all
%   clusters, of the within-cluster sums of point-to-cluster-centroid
%   distances.  Rows of X correspond to points, columns correspond to
%   variables.  KMEANS returns an N-by-1 vector IDX containing the cluster
%   indices of each point.  By default, KMEANS uses squared Euclidean
%   distances.
%
%   See the help for KMEANS for more details on output and input option
%   parameter name/value pairs.  
% 
%   See also KMEANS, DATAMATRIX/PDIST, DATAMATRIX/PCA.

%   Copyright 2008-2012 The MathWorks, Inc.


try
    [idx, varargout{1:nargout-1}] = kmeans(X.Matrix, k, varargin{:});
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','kmeans', ME);
end
end %DataMatrix/kmeans
