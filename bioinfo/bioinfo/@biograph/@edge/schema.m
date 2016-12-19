function schema
%

% Copyright 2003-2010 The MathWorks, Inc.

%% add edge class to package
bgpk = findpackage('biograph');
cls = schema.class(bgpk,'edge');

%% public properties
p = schema.prop(cls,'ID','string');                  p.FactoryValue =  '';
                                                     p.SetFunction = @verifyUniqueID;
p = schema.prop(cls,'Label','string');               p.FactoryValue =  '';
p = schema.prop(cls,'Description','string');         p.FactoryValue =  '';
p = schema.prop(cls,'Weight','double');              p.FactoryValue = 1;
p = schema.prop(cls,'LineWidth','double');           p.FactoryValue = 0.5;
                                                     p.SetFunction = @verifyPositiveValue;
p = schema.prop(cls,'LineColor','color');            p.FactoryValue = [.5 .5 .5];
p = schema.prop(cls,'UserData','MATLAB array');      p.FactoryValue = [];

%% private properties
p = schema.prop(cls,'Visible','bool');               p.FactoryValue = true;
                                                     p.Visible = 'off';
p = schema.prop(cls,'ControlPoints','MATLAB array'); p.FactoryValue = [];
                                                     p.Visible = 'off';
p = schema.prop(cls,'FromNode','handle');            p.FactoryValue = handle([]);
                                                     p.Visible = 'off';
p = schema.prop(cls,'ToNode','handle');              p.FactoryValue = handle([]);
                                                     p.Visible = 'off';
p = schema.prop(cls,'idx','long');                   p.FactoryValue = 0;
                                                     p.Visible = 'off';                                                     

%% handles to HG
p = schema.prop(cls,'hgLine','MATLAB array');        p.FactoryValue = [];
                                                     p.Visible = 'off';
p = schema.prop(cls,'hgArrow','MATLAB array');       p.FactoryValue = [];
                                                     p.Visible = 'off';
p = schema.prop(cls,'isSelected','bool');            p.FactoryValue = false;  
                                                     p.Visible = 'off';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function valueStored = verifyPositiveValue(obj, valueProposed)
if valueProposed>0
    valueStored = valueProposed;
else
    error(message('bioinfo:biographedgeschema:ParameterNotPositive'))
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function valueStored = verifyUniqueID(obj, valueProposed)
if isempty(obj.up) %object being copied
    valueStored = valueProposed;
    return
end
h = valueProposed(:)<' ' | (valueProposed(:)>=127 & valueProposed(:)<=159);
if any(h)
    warning(message('bioinfo:biographedgeschema:ControlCharacter'))
    valueProposed(h) = ' ';
end

if isempty(valueProposed)
   error(message('bioinfo:biographedgeschema:ParameterEmpty')) 
end
current = setdiff(obj.up.edges,obj);
if isempty(current)
    valueStored = valueProposed;
elseif numel(current)==1 
    if isequal(valueProposed,current.ID)
        error(message('bioinfo:biographedgeschema:IDNotUnique'))
    else
        valueStored = valueProposed;
    end
else
    if ismember(valueProposed,get(current,'ID'))
        error(message('bioinfo:biographedgeschema:IDNotUnique'))
    else
        valueStored = valueProposed;
    end
end
