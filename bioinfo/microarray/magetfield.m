function data = magetfield(maStruct,fieldname)
%MAGETFIELD extracts data for a given field from a microarray structure.
%
%   MAGETFIELD(MASTRUCT,FIELDNAME) extracts data column FIELDNAME from a
%   microarray structure MASTRUCT.
%
%   Examples:
%
%       maStruct = gprread('mouse_a1wt.gpr');
%       cy5data = magetfield(maStruct,'F635 Median');
%       cy3data = magetfield(maStruct,'F532 Median');
%       mairplot(cy5data,cy3data,'title','R vs G IR plot');
%
%   See also AGFEREAD, GPRREAD, IMAGENEREAD, MABOXPLOT, MAIRPLOT, MALOGLOG,
%   MALOWESS, MOUSEDEMO, SPTREAD.

% Copyright 2003-2008 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);

if ~isstruct(maStruct) || ~isfield(maStruct,'ColumnNames') || ~isfield(maStruct,'Data')
    error(message('bioinfo:magetfield:NotStruct'))
end

theCol = strmatch(fieldname,maStruct.ColumnNames);
if isempty(theCol)
    error(message('bioinfo:magetfield:UnknownField', fieldname))
end
if numel(theCol) > 1
    theCol = strmatch(fieldname,maStruct.ColumnNames,'exact');
    if numel(theCol) == 0
        error(message('bioinfo:magetfield:AmbiguousField', fieldname))
    end
end
try
    data = maStruct.Data(:,theCol);
catch allExceptions
    error(message('bioinfo:magetfield:BadColumn', fieldname))
end


