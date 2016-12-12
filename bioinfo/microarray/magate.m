function varargout = magate(varargin)
%MAGATE Gateway routine to call bioinfo/microarray private functions.
%
%   [OUT1, OUT2,...] = MAGATE(FCN, VAR1, VAR2,...) calls FCN in 
%   the bioinfo/microarray private directory with input arguments VAR1, 
%   VAR2,... and returns the output, OUT1, OUT2,....
%
%   MAGATE is an internal method used by bioinfo/microarray. This method should
%   not be called directly by users and will likely change in a future
%   release.
%

%   Copyright 2006 The MathWorks, Inc. 

if nargin == 0
   error(message('bioinfo:magate:NoInput'));
end

nout = nargout;
if nout==0,
   feval(varargin{:});
else
   [varargout{1:nout}] = feval(varargin{:});
end
