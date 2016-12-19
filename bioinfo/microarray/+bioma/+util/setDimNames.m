function newNames = setDimNames(names, numdims, dim, defName, allowDup)
% SETDIMNAMES Sets the names for rows and columns.
% 
% NEWNAMES = SETDIMNAMES(NAMES, NUMDIMS, DIM, DEFNAME, ALLOWDUP) check if the NAMES
% for row or columns are unique and match the number of rows or columns,
% NUMDIMS, in object, and return the new names. DEFNAME is used for
% construct an error message.

%   Copyright 2009-2012 The MathWorks, Inc.

% Errors given by this helper do not comply with error ids for classes,
% therefore when used it should be within a try-catch block and a correct
% message id should be given.

% defName can be any of the following: Element, Feature, Sample, DMName
if (iscellstr(names) || isnumeric(names) || ischar(names)) && isvector(names)
    if ischar(names) 
        if size(names, 1) > 1
            names = cellstr(names);
        else
            newNames = bioma.util.appendUniqueNumToNames(names, numdims, 1);
            if dim == 2
                newNames = newNames';
            end
            return;
        end 
    elseif isnumeric(names)
        names = cellstr(num2str(names(:)));
    end
    
    if ~allowDup && (numel(unique(names))~=numel(names))
        error(message('bioinfo:setDimNames:DuplicatedNames', defName))
    end
    
    numNames = numel(names);
    if numNames == numdims
        newNames = names(:);
        if dim == 2
            newNames = newNames';
        end
    else
        switch defName
            case 'Element'
                error(message('bioinfo:setDimNames:NumberOfElementNamesNotMatch', numNames, numdims))
            case 'Feature'
                error(message('bioinfo:setDimNames:NumberOfFeatureNamesNotMatch', numNames, numdims))
            case 'Sample'
                error(message('bioinfo:setDimNames:NumberOfSampleNamesNotMatch', numNames, numdims))
            case 'DMName'
                error(message('bioinfo:setDimNames:NumberOfDMNameNamesNotMatch', numNames, numdims))
        end
    end
elseif islogical(names)
    %= Default unique names.
   newNames = bioma.util.appendUniqueNumToNames(defName, numdims, 1);
elseif isempty(names)
    newNames = names;
else
    switch defName
        case 'Element'
           error(message('bioinfo:setDimNames:InvalidElementNames'))
        case 'Feature'
           error(message('bioinfo:setDimNames:InvalidFeatureNames'))
        case 'Sample'
           error(message('bioinfo:setDimNames:InvalidSampleNames'))
        case 'DMName'
           error(message('bioinfo:setDimNames:InvalidDMNameNames'))
    end
end
end %setDimNames
