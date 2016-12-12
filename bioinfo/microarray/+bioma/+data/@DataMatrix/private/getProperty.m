function [varargout] = getProperty(a,name)
%GETPROPERTY Get a DataMatrix property.

%   Copyright 2008 The MathWorks, Inc. 


% We may be given a name, or a subscript expression that starts with a
% '.name' subscript.  Get the name and validate it in any case.
if isstruct(name)
    s = name;
    if s(1).type == '.'
        name = s(1).subs;
    else
        error(message('bioinfo:DataMatrix:getProperty:InvalidSubscript'));
    end
    haveSubscript = true;
else
    haveSubscript = false;
end

if ~(ischar(name) && isvector(name) && (size(name,1)==1))
    error(message('bioinfo:DataMatrix:getProperty:InvalidPropertyName'));
end
propertyNames = [fieldnames(a); {'Matrix'}];
k = bioinfoprivate.pvpair(name, [], propertyNames, 'DataMatrix:getProperty');
name = propertyNames{k};

% Get the property out of the DataMatrix.  Some properties "should be" 0x0
% cell arrays if they're not present, but can become 1x0 or 0x1 through
% subscripting.  Make those cosmetically nice.
switch name
    case 'Name'
        p = a.Name;
    case 'RowNames'
        p = a.RowNames;
    case 'ColNames'
        p = a.ColNames;
    case 'NRows'
        p = a.NRows;
    case 'NCols'
        p = a.NCols;
    case 'NDims'
        p = a.NDims;
    case 'ElementClass'
        p = class(a.Matrix);
    case 'Matrix'
        p = a.Matrix;
end

if haveSubscript && ~isscalar(s)
    % If there's cascaded subscripting into the property, let the property's
    % subsasgn handle the reference.  This may return a comma-separated list,
    % so ask for and assign to as many outputs as we're given.  If there's no
    % LHS to the original expression (nargout==0), this only assigns one
    % output and drops everything else in the CSL.
    [varargout{1:nargout}] = subsref(p,s(2:end));
else
    % If there's no cascaded subscripting, only ever assign the property
    % itself.
    varargout{1} = p;
end
end %getProperty private function
