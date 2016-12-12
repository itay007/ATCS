function schema
%schema

%   Copyright 2003-2006 The MathWorks, Inc.


pkg              = findpackage('seqviewer');
c                = schema.class(pkg, 'MSAUtil');
m                = schema.method(c, 'msaMatlab', 'static'); %#ok
