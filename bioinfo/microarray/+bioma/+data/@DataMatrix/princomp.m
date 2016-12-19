function varargout = princomp(X, varargin)
%PRINCOMP Overload Principal Components Analysis for DataMatrix objects.
% 
%   PRINCOMP is not recommended. Use PCA instead.
%
%   COEFF = PRINCOMP(X) performs principal components analysis on the
%   N-by-P DataMatrix object X, and returns the principal component
%   coefficients, also known as loadings.  Rows of X correspond to
%   observations, columns to variables.  COEFF is a P-by-P matrix, each
%   column containing coefficients for one principal component.  The
%   columns are in order of decreasing component variance.
%
%   See the help for PRINCOMP for more details on output and additional
%   inputs.
% 
%   See also PCA, DATAMATRIX/KMEANS, DATAMATRIX/PCA, DATAMATRIX/PDIST.

%   Copyright 2008-2012 The MathWorks, Inc.


try
    [varargout{1:nargout}] = princomp(X.Matrix, varargin{:});
catch ME
    bioinfoprivate.bioclsrethrow('DataMatrix','princomp', ME);
end
end %DataMatrix/princomp
