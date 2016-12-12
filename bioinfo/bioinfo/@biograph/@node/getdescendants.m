function nodes = getdescendants(h,gen)
%GETDESCENDANTS walks through a biograph to find descendants
%
%   NODES = GETDESCENDANTS(BIOGRAPH_NODE) returns the current node and its
%   direct descendants.  
%
%   NODES = GETDESCENDANTS(BIOGRAPH_NODE, GEN)  returns all the
%   BIOGRAPH_NODE descendants up to GEN generations (including itself).
%
%   Example: 
%
%       % Create a BIOGRAPH object.
%       cm = [0 1 1 0 0;1 0 0 1 1;1 0 0 0 0;0 0 0 0 1;1 0 1 0 0];
%       bg = biograph(cm)
%
%       % Find one generation of descendants for node 4.
%       desNodes = getdescendants(bg.nodes(4));
%       set(desNodes,'Color',[1 .7 .7])
%       view(bg)
%
%       % Find two generations of descendants for node 4.
%       desNodes = getdescendants(bg.nodes(4),2);
%       set(desNodes,'Color',[.7 1 .7])
%       view(bg)
%
%   See also BIOGRAPH/BIOGRAPH, BIOGRAPH.BIOGRAPH/DOLAYOUT,
%   BIOGRAPH.BIOGRAPH/GETEDGESBYNODEID, BIOGRAPH.BIOGRAPH/GETNODESBYID,
%   BIOGRAPH.BIOGRAPH/VIEW, BIOGRAPH.NODE/GETANCESTORS,
%   BIOGRAPH.NODE/GETRELATIVES, GET, SET.

% Copyright 2003-2006 The MathWorks, Inc.


if nargin == 1
    gen = 1;
else
    gen = abs(round(gen)); % make sure it is a positive integer
end

bg = h(1).up;
nl = ismember(bg.Nodes,h);
finddescendants(gen);
nodes = bg.Nodes(nl);

    function finddescendants(gen)
        if gen == 0
            return
        end
        if gen>1
            finddescendants(gen-1);
        end
        nl(any(bg.from(:,nl),2)) = true;
    end

end
