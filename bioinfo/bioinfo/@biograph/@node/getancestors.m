function nodes = getancestors(h,gen)
%GETANCESTORS walks through a biograph to find ancestors
%
%   NODES = GETANCESTORS(BIOGRAPH_NODE) returns the current node and its
%   direct ancestors.   
%
%   NODES = GETANCESTORS(BIOGRAPH_NODE, GEN) returns all the BIOGRAPH_NODE
%   ancestors up to GEN generations (including itself).
%
%   Example: 
%
%       % Create a BIOGRAPH object.
%       cm = [0 1 1 0 0;1 0 0 1 1;1 0 0 0 0;0 0 0 0 1;1 0 1 0 0];
%       bg = biograph(cm)
%
%       % Find one generation of ancestors for node 2.
%       ancNodes = getancestors(bg.nodes(2));
%       set(ancNodes,'Color',[1 .7 .7])
%       view(bg)
%
%       % Find two generations of ancestors for node 2.
%       ancNodes = getancestors(bg.nodes(2),2);
%       set(ancNodes,'Color',[.7 1 .7])
%       view(bg)
%
%   See also BIOGRAPH/BIOGRAPH, BIOGRAPH.BIOGRAPH/DOLAYOUT,
%   BIOGRAPH.BIOGRAPH/GETEDGESBYNODEID, BIOGRAPH.BIOGRAPH/GETNODESBYID,
%   BIOGRAPH.BIOGRAPH/VIEW, BIOGRAPH.NODE/GETDESCENDANTS,
%   BIOGRAPH.NODE/GETRELATIVES, GET, SET.

% Copyright 2003-2006 The MathWorks, Inc.


if nargin == 1
    gen = 1;
else
    gen = abs(round(gen)); % make sure it is a positive integer
end

bg = h(1).up;
nl = ismember(bg.Nodes,h);
findancestors(gen);
nodes = bg.Nodes(nl);

     function findancestors(gen)
        if gen == 0
            return
        end
        if gen>1
            findancestors(gen-1);
        end
        nl(any(bg.to(:,nl),2)) = true;
     end
end


