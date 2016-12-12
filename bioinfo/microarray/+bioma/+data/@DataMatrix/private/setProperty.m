function obj = setProperty(obj,name,p)
%SETPROPERTY Set a DataMatrix property.

%   Copyright 2008-2012 The MathWorks, Inc. 


%== Subscript expression
hasSubscript = false;
if isstruct(name)
    s = name;
    if s(1).type == '.'
        name = s(1).subs;
    else
        error(message('bioinfo:DataMatrix:setProperty:InvalidSubscript'));
    end
    hasSubscript = true;
end

if ~(ischar(name) && isvector(name) && (size(name,1)==1))
    error(message('bioinfo:DataMatrix:setProperty:InvalidPropertyName'));
end
propertyNames = [fieldnames(obj); {'Matrix'}];
k = bioinfoprivate.pvpair(name, [], propertyNames, 'DataMatrix:setProperty');
name = propertyNames{k};

if hasSubscript && ~isscalar(s)
    oldp = getProperty(obj, name);
    p = subsasgn(oldp, s(2:end), p);
else
    
end
    
% Assign the new property value into the DataMatrix.
switch (k)
    case 1 % Names
        obj = setName(obj,p);
    case 2 % RowNames
        if any([obj.NRows obj.NCols]==0)
            obj.RowNames = [];
        else
            obj.RowNames = setDimNames(p, obj.NRows, 1);
        end
    case 3 % ColumnNames
        if any([obj.NRows obj.NCols]==0)
            obj.ColNames = [];
        else
            obj.ColNames = setDimNames(p, obj.NCols, 2);
        end
    otherwise
        error(message('bioinfo:DataMatrix:setProperty:SetProhibited', name));
end
end % DataMatrix/setProperty

%------------ Helper Functions-----------------
function obj = setName(obj, x)
if isempty(x)
    return;
end

if ischar(x)
    obj.Name = x;
elseif iscellstr(x)
    obj.Name = x{1};
else
    error(message('bioinfo:DataMatrix:setProperty:BadNameInputFormat'))
end
end % setName

function newNames = setDimNames(names, numdims, dim)
% Check if the names for row or columns are unique and match the number of
% rows or columns in data matrix       

if (iscellstr(names) || isnumeric(names) || ischar(names)) && isvector(names)
    if ischar(names) 
        if size(names, 1) > 1
            names = cellstr(names);
        else
            newNames = cellstr([ repmat(names,numdims,1)  num2str((1:numdims)')]);
            if dim == 2
                newNames = newNames';
            end
            return;
        end 
    elseif isnumeric(names)
        names = cellstr(num2str(names(:)));
    end
    
    numNames = numel(names);
    if dim == 1
        if numNames == numdims
            newNames = names(:);
        else
            error(message('bioinfo:DataMatrix:setProperty:NumberOfRowNamesNotMatch', numNames, numdims))
        end
    elseif dim == 2
        if numNames == numdims
            newNames = names(:)';
        else
            error(message('bioinfo:DataMatrix:setProperty:NumberOfColumnNamesNotMatch', numNames, numdims))
        end
    end
elseif islogical(names)
    %default unique names.
   if dim == 1
       newNames = cellstr([ repmat('Row',numdims,1)  num2str((1:numdims)')]);
   else
       newNames = cellstr([ repmat('Col',numdims,1)  num2str((1:numdims)')]);
   end
elseif isempty(names)
    newNames = names;
else
    switch dim
        case 1
            error(message('bioinfo:DataMatrix:setProperty:InvalidRowNames'))
        case 2
            error(message('bioinfo:DataMatrix:setProperty:InvalidColumnNames'))
    end
end
end %setDimNames

