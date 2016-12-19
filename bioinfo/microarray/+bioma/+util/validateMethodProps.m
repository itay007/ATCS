function varargout = validateMethodProps(obj, name)
% VALIDATEMETHODPROPS Validate input name is a method or public property name.
% 
%  VALIDATEMETHODPROPS(OBJ,NAME, CLSNAME, FUNNAME) validates if NAME is a
%  method or public property name of OBJ class. Useful in subsref and
%  subsasgn.

%   Copyright 2009 The MathWorks, Inc.


methodNames = methods(obj);
propNames = properties(obj);
methodsProp = [methodNames; propNames];
if ~any(strcmp(name, methodsProp))
    error(message('bioinfo:validateMethodProps:InvalidMethodOrProps', name, class( obj )));
end

if nargout >= 1
    varargout{1} = methodNames;
end 

if nargout >= 2
    varargout{2} = propNames;
end

end
