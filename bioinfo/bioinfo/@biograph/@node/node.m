function h = node(hPar,id,i)
%NODE class constructor.


% Copyright 2003-2006 The MathWorks, Inc.

h = biograph.node;
h.idx   = i;
h.ID    = id; 

% by connecting up after setting the ID's we remove the uniqueness ID
% checking to save time
connect(h,hPar,'up') 



