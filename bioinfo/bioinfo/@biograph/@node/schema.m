function schema
%

% Copyright 2003-2010 The MathWorks, Inc.

%% add node class to package
bgpk = findpackage('biograph');
cls = schema.class(bgpk,'node');

%% public properties
p = schema.prop(cls,'ID','string');                  p.FactoryValue =  '';
                                                     p.SetFunction = @verifyUniqueID;
p = schema.prop(cls,'Label','string');               p.FactoryValue =  '';       
p = schema.prop(cls,'Description','string');         p.FactoryValue =  '';
p = schema.prop(cls,'Position','MATLAB array');      p.FactoryValue = [];
                                                     p.SetFunction = @verifyProperPosition;
p = schema.prop(cls,'Shape','BiographShape');        p.FactoryValue = 'box';
p = schema.prop(cls,'Size','MATLAB array');          p.FactoryValue = [10 10];
                                                     p.SetFunction = @verifyProperSize;
p = schema.prop(cls,'Color','color');                p.FactoryValue = [1 1 0.7];
p = schema.prop(cls,'LineWidth','double');           p.FactoryValue = 1;
                                                     p.SetFunction = @verifyPositiveValue;
p = schema.prop(cls,'LineColor','color');            p.FactoryValue = [0.3 0.3 1];
p = schema.prop(cls,'FontSize','double');            p.FactoryValue = 9;
                                                     p.SetFunction = @verifyPositiveValue;
p = schema.prop(cls,'TextColor','color');            p.FactoryValue = [0 0 0];
p = schema.prop(cls,'UserData','MATLAB array');      p.FactoryValue = [];

%% private properties   
p = schema.prop(cls,'Visible','bool');               p.FactoryValue = true;
                                                     p.Visible = 'off';
p = schema.prop(cls,'idx','int');                    p.FactoryValue = 0;
                                                     p.Visible = 'off';                                                    

%% handles to HG
p = schema.prop(cls,'hgPatch','MATLAB array');      p.FactoryValue = [];
                                                    p.Visible = 'off';
p = schema.prop(cls,'hgText','MATLAB array');       p.FactoryValue = [];
                                                    p.Visible = 'off';
p = schema.prop(cls,'isSelected','bool');           p.FactoryValue = false;  
                                                    p.Visible = 'off';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function valueStored = verifyProperPosition(obj, valueProposed)
if ~isnumeric(valueProposed) || ~isreal(valueProposed) || ...
    numel(valueProposed)~=2  
    error(message('bioinfo:biographnodeschema:ParameterNotTwoElemReal'))
else
    valueStored = valueProposed([1 2]);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function valueStored = verifyProperSize(obj, valueProposed)
if ~isnumeric(valueProposed) || ~isreal(valueProposed) || ...
    numel(valueProposed)~=2  || ~all(valueProposed>0)
    error(message('bioinfo:biographnodeschema:ParameterNotTwoElemPositive'))
else
    valueStored = valueProposed([1 2]);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function valueStored = verifyPositiveValue(obj, valueProposed)
if valueProposed>0
    valueStored = valueProposed;
else
    error(message('bioinfo:biographnodeschema:ParameterNotPositive'))
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function valueStored = verifyUniqueID(obj, valueProposed)
if isempty(obj.up) %object being copied
    valueStored = valueProposed;
    return
end
h = valueProposed(:)<' ' | (valueProposed(:)>=127 & valueProposed(:)<=159);
if any(h)
    warning(message('bioinfo:biographnodeschema:ControlCharacter'))
    valueProposed(h) = ' ';
end

if isempty(valueProposed)
   error(message('bioinfo:biographnodeschema:ParameterEmpty')) 
end
current = setdiff(obj.up.nodes,obj);
if isempty(current)
    valueStored = valueProposed;
elseif numel(current)==1 
    if isequal(valueProposed,current.ID)
        error(message('bioinfo:biographnodeschema:IDNotUnique'))
    else
        valueStored = valueProposed;
    end
else
    if ismember(valueProposed,get(current,'ID'))
        error(message('bioinfo:biographnodeschema:IDNotUnique'))
    else
        valueStored = valueProposed;
    end
end
