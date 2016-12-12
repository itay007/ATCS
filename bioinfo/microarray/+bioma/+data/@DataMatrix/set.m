function b = set(a,varargin)
%SET  Set DataMatrix object property value.
% 
%   B = SET(A,'PropertyName',VALUE) sets the property 'PropertyName' of the
%   DataMatrix object A to the value VALUE. 
%  
%   B = SET(A,'Property1',Value1,'Property2',Value2,...) sets multiple
%   property values of a DataMatrix object A with a single statement.
%
%   SET(A,'Property') displays possible values for the specified property
%   in DataMatrix object A.
%
%   SET(A) displays all properties of DataMatrix A and their possible
%   values. 
%
%   See also DataMatrix/GET.

%   Copyright 2008-2012 The MathWorks, Inc. 


if nargin < 3
    propertyNames = fieldnames(a);
    propertyNames = propertyNames(1:3);
    propVals   = cell2struct(cell(size(propertyNames)), propertyNames, 1);
    propDescrs = cell2struct(cell(size(propertyNames)), propertyNames, 1);
    
    propVals.Name   = '';
    propDescrs.Name = 'A DataMatrix''s ''Name'' property does not have a fixed set of values.';
    propVals.RowNames   = {};
    propDescrs.RowNames = 'Empty, a cell array of strings or a numeric vector.';
    propVals.ColNames   = {};
    propDescrs.ColNames = 'Empty, a cell array of strings or a numeric vector.';
 
    if nargin == 2
        
        name = varargin{1};
        if ~(ischar(name) && isvector(name) && (size(name,1)==1))
            error(message('bioinfo:DataMatrix:set:InvalidPropertyName'));
        end
        allPropertyNames = fieldnames(a);
        k = bioinfoprivate.pvpair(name, [], allPropertyNames, 'DataMatrix:set');
        name = allPropertyNames{k};
                       
        if isempty(find(strcmp(name, propertyNames), 1))
            error(message('bioinfo:DataMatrix:set:SetProhibited', name));
        end
        if nargout == 1
            b = propVals.(name);
        else
            disp(propDescrs.(name));
        end        
    else
        if nargout == 1
            b = propVals;
        else
            disp(propDescrs);
        end
    end
    
elseif mod(nargin,2) == 1 %name/value pairs
        b = a;
        for i = 1:(nargin-1)/2
            name = varargin{2*i-1};
            p = varargin{2*i};
            b = setProperty(b,name,p);
        end
else
    error(message('bioinfo:DataMatrix:set:WrongNumberArgs'));
end
end % DataMatrix/set
