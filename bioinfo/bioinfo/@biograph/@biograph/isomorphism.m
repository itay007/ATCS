function varargout =  isomorphism(varargin)
%ISOMORPHISM finds an isomorphism between two graphs.
%
%  ISOMORPHISM method for BIOGRAPH objects extracts the adjacency matrix
%  and calls GRAPHISOMORPHISM function. All other input and output
%  arguments are analogous to GRAPHISOMORPHISM.
% 
%  See also graphisomorphism.

% Copyright 2006 The MathWorks, Inc.


[varargout{1:nargout}] = graphisomorphism(getmatrix(varargin{1}),getmatrix(varargin{2}),varargin{3:end});
