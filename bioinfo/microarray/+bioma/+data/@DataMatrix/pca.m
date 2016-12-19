function varargout = pca(X, varargin)
%PCA Overload Principal Components Analysis for DataMatrix objects.
%
%   COEFF = PCA(X) performs principal components analysis on the
%   N-by-P DataMatrix object X, and returns the principal component
%   coefficients, also known as loadings.  Rows of X correspond to
%   observations, columns to variables.  COEFF is a P-by-P matrix, each
%   column containing coefficients for one principal component.  The
%   columns are in order of decreasing component variance.
%
%   See the help for PCA for more details on output and additional
%   inputs.
% 
%   See also PCA, DATAMATRIX/KMEANS, DATAMATRIX/PCA, DATAMATRIX/PDIST.

%   Copyright 2008-2012 The MathWorks, Inc.


try
    [varargout{1:nargout}] = pca(X.Matrix, varargin{:});
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','pca', ME);
end
end %DataMatrix/pca
