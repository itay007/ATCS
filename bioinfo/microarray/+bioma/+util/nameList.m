function out = nameList(obj, propName, idx, names, varargin)
%NAMELIST Set and retrieve a list of names of an object property.
%   NAMES = NAMELIST(OBJ, PROPNAME, FNCNAME, CLSNAME) returns a cell array
%   of names in a property, PROPNAME, in an object OBJ. The property,
%   PROPNAME, must contain a cell array of strings. FNCNAME is the method
%   name that call this function. CLSNAME is the class name of OBJ. FNCNAME
%   and CLSNAME are needed for throw error messages.
%
%   VNAMES = NAMELIST(OBJ, ..., IDX) returns the names in PROPNAME
%   specified by IDX of OBJ. IDX can be a positive integer, a vector of
%   positive integers, a string specifying a name, a cell array containing
%   one or more names, or a logical vector.
%
%   B = NAMELISTS(OBJ, ...,IDX, NAMES) returns an object B with the names
%   in specified property PROPNAME set to NAMES. The number of names in
%   NAMES must be equal to the number of names specified by I. NAMES can be
%   a cell array of strings. a character array, or a numeric vector. NAMES
%   can also be a single string as a prefix for names; numbers will be
%   appended to the prefix. NAMES can also be a logical true or false
%   (default); if true, default unique names will be assigned to the name
%   list in PROPNAME.
% 
%   VNAMES = NAMELIST(OBJ, ..., IDX, NAMES, CANSETEMPTY) specify if the
%   names in PROPNAME can be set to empty. Default is FALSE.
% 
%   VNAMES = NAMELIST(OBJ, ..., IDX, NAMES, CANSETEMPTY, PREFIXSTR) specify
%   the prefix string if NAMES is a logical. The default is TMP.

%   Copyright 2009 The MathWorks, Inc.


%== Name List

try
    nameList = get(obj, propName);
catch ME %#ok
    nameList = obj.(propName);
end

nameListLen = numel(nameList);
canSetEmptyFlag = false;
prefixStr = 'TMP';

if nargin >= 5
    canSetEmptyFlag = varargin{1};
end

if nargin == 6
    prefixStr = varargin{2};
end

if nargin == 2
    out = nameList;
elseif nargin == 3
    if isempty(idx)
        idx = ':';
    end
    
    indices = bioma.util.findLabelIndices(nameList,idx);
    
    if any(indices == 0) && ~islogical(indices)
        bioma.util.errorUnknownNames(idx(indices == 0))
    elseif isempty(indices) && ~isempty(nameList)
        bioma.util.errorUnknownNames(idx)
    end
    
    if ~isempty(nameList)
        out = nameList(indices);
    else
        out = [];
    end
elseif nargin >= 4 %== Set names
    if isempty(idx) || bioma.util.isColon(idx)
        out = set(obj, propName, names);
    else
        [indices, numIndices] = bioma.util.findLabelIndices(nameList, idx);
        if any(indices == 0) && ~islogical(indices)
            bioma.util.errorUnknownNames(idx(indices == 0))
        elseif isempty(indices) && ~isempty(names)
            bioma.util.errorUnknownNames(idx)
        end
        %== NAMES is empty case
        if isempty(names)
            if canSetEmptyFlag           
                if numIndices == nameListLen
                    obj = set(obj, PropName, []);
                else
                    nameList(indices) = cellstr(repmat('', 1, numIndices));
                    obj = set(obj, PropName, nameList);
                end
                out = obj;
                return;
            else
                error(message('bioinfo:nameList:CanNotSetEmptyNames'))
            end
        end
        
        if (iscellstr(names) || isnumeric(names) || ischar(names)) && isvector(names)
            if ischar(names)
                idxStr = bioma.util.indicesToStr(indices, nameListLen);
                names = cellstr([ repmat(names, numIndices, 1)  idxStr]);
            elseif isnumeric(names)
                names = cellstr(num2str(names(:)));
            end
            
            if numel(names) ~= numIndices
                error(message('bioinfo:nameList:NumIndicesNotMatchNumNames'))
            end
        elseif islogical(names)
            %default unique names.
            idxStr = bioma.util.indicesToStr(indices, nameListLen);
            names = cellstr([ repmat(prefixStr, numIndices, 1)  idxStr]);
        else
            error(message('bioinfo:nameList:InvalidNames'))
        end
        
        nameList(indices) = names;
        obj = set(obj, propName, nameList);
        out = obj;
    end
end
end
