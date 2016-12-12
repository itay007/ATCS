function p = get(a,name)
%GET Get a DataMatrix object property.
% 
%   GET(A) displays the list of properties and their current values of
%   DataMatrix object A.
%
%   S = GET(A) returns all properties of DataMatrix object A in a scalar
%   structure, in which each field name is a property of A, and each field
%   contains the value of that property.
%
%   P = GET(A,'PropertyName') returns the property value of 'PropertyName'
%   of DataMatrix object A. If 'PropertyName' is replaced by a 1-by-N or
%   N-by-1 cell array of strings containing property names, GET will return
%   a 1-by-N or N-by-1 cell array of values.
%
%   See also DATAMATRIX/SET.

%   Copyright 2008 The MathWorks, Inc. 


if nargin == 1
    propNames = fieldnames(a);
    s   = cell2struct(cell(size(propNames)), propNames, 1);    
    for i = 1:length(propNames)
        s.(propNames{i}) = getProperty(a,propNames{i});
    end
    
    if nargout == 1
        p = s;
    else
        disp(s);
    end
elseif nargin == 2
    if iscellstr(name)
        p = cell(1,numel(name));
        for i = 1:length(p)
            p{i} = getProperty(a,name{i});
        end
    else
        p = getProperty(a,name);
    end
end
end % get method
