function schema
%schema

%   Copyright 2005-2006 The MathWorks, Inc.


pkg              = findpackage('seqviewer');
c                = schema.class(pkg, 'SVUtil');

m                = schema.method(c, 'svMatlab', 'static'); %#ok
m                = schema.method(c, 'svSaveFile', 'static'); %#ok
m                = schema.method(c, 'svExportToWS', 'static'); %#ok
m                = schema.method(c, 'svImportFromWS', 'static'); %#ok
m                = schema.method(c, 'bioinfoabout', 'static'); %#ok
