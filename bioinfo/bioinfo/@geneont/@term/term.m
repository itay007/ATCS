function GOterm = term(hpar) 
%TERM class constructor.


% Copyright 2003-2006 The MathWorks, Inc.
GOterm = geneont.term; % instantiate object
GOterm.connect(hpar,'up');

