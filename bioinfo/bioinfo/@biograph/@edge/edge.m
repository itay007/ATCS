function h = edge(hPar,from,to,id,len,i)
%EDGE class constructor.


% Copyright 2003-2006 The MathWorks, Inc.

h = biograph.edge;
h.idx   = i;
h.ID    = id; 
h.Weight = len; 

% by connecting up after setting the ID's we remove the uniqueness ID
% checking to save time
connect(h,hPar,'up') 

h.FromNode = from;
h.ToNode = to;





