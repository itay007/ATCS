function nodes = getrelatives(h,gen)
%GETRELATIVES walks through a biograph to find relatives
%
%   NODES = GETRELATIVES(BIOGRAPH_NODE) returns the current node and its
%   direct relatives.   
%
%   NODES = GETRELATIVES(BIOGRAPH_NODE, GEN) returns all the BIOGRAPH_NODE
%   relatives up to GEN generations (including itself). 
%
%   Example:
%
%       % Create a BIOGRAPH object.
%       cm = [0 1 1 0 0;1 0 0 1 1;1 0 0 0 0;0 0 0 0 1;1 0 1 0 0];
%       bg = biograph(cm)
%
%       % Find all nodes interacting with node 1.
%       intNodes = getrelatives(bg.nodes(1));
%       set(intNodes,'Color',[.7 .7 1])
%       view(bg)
%
%   See also BIOGRAPH/BIOGRAPH, BIOGRAPH.BIOGRAPH/DOLAYOUT,
%   BIOGRAPH.BIOGRAPH/GETEDGESBYNODEID, BIOGRAPH.BIOGRAPH/GETNODESBYID,
%   BIOGRAPH.BIOGRAPH/VIEW, BIOGRAPH.NODE/GETANCESTORS,
%   BIOGRAPH.NODE/GETDESCENDANTS, GET, SET.

% Copyright 2003-2006 The MathWorks, Inc.


if nargin == 1
    gen = 1;
else
    gen = abs(round(gen)); % make sure it is a positive integer
end

nodes = unique([getancestors(h,gen);getdescendants(h,gen)]);

