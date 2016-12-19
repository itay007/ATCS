function set(tr,varargin) %#ok
%SET  Set object properties of a phylogenetic tree object.
%
%  Properties in a phylogenetic tree object cannot be manually set.
%  A PHYTREE object must be created by its constructor method PHYTREE
%  or by using one of the functions: PHYTREEREAD, SEQLINKAGE, SEQNEIGHJOIN.
%
%  See also: PHYTREE, PHYTREEREAD, SEQLINKAGE, SEQNEIGHJOIN.

% Copyright 2003-2006 The MathWorks, Inc.


error(message('bioinfo:phytree:set:NotAllowedMethod'))
