function schema
%

% Copyright 2003-2010 The MathWorks, Inc.

%% set layout types
if isempty(findtype('BiographLayouts'))
    schema.EnumType('BiographLayouts',{'hierarchical','equilibrium','radial'});
end
%% set layout types
if isempty(findtype('BiographEdgeType'))
    schema.EnumType('BiographEdgeType',{'straight','curved','segmented'});
end 
%% set text in Nodes
if isempty(findtype('BiographTextInNodes'))
    schema.EnumType('BiographTextInNodes',{'none','id','label'});
end 

%% add biograph class to package
bgpk = findpackage('biograph'); 
cls = schema.class(bgpk,'biograph');  
 
%% public properties
p = schema.prop(cls,'ID','string');                  p.FactoryValue =  '';
p = schema.prop(cls,'Label','string');               p.FactoryValue =  '';
p = schema.prop(cls,'Description','string');         p.FactoryValue =  '';

p = schema.prop(cls,'LayoutType','BiographLayouts'); p.FactoryValue = 'hierarchical';
p = schema.prop(cls,'LayoutScale','double');         p.FactoryValue = 1;
                                                     p.SetFunction = @verifyPositiveValue;

p = schema.prop(cls,'Scale','double');               p.FactoryValue = 1;
                                                     p.SetFunction = @verifyPositiveValue;
p = schema.prop(cls,'NodeAutoSize','on/off');        p.FactoryValue = 'on';                                                     
p = schema.prop(cls,'ShowTextInNodes','BiographTextInNodes'); p.FactoryValue = 'label';


p = schema.prop(cls,'EdgeType','BiographEdgeType');  p.FactoryValue = 'curved';
p = schema.prop(cls,'EdgeTextColor','color');        p.FactoryValue = [0 0 0];
p = schema.prop(cls,'ShowArrows','on/off');          p.FactoryValue = 'on';
p = schema.prop(cls,'ArrowSize','double');           p.FactoryValue = 8;
                                                     p.SetFunction = @verifyProperSize;
p = schema.prop(cls,'ShowWeights','on/off');         p.FactoryValue = 'off';
p = schema.prop(cls,'EdgeFontSize','double');        p.FactoryValue = 8;
                                                     p.SetFunction = @verifyPositiveValue;
 
p = schema.prop(cls,'NodeCallbacks','MATLAB array');  p.FactoryValue = @(node) inspect(node);
p = schema.prop(cls,'EdgeCallbacks','MATLAB array');  p.FactoryValue = @(edge) inspect(edge);
p = schema.prop(cls,'CustomNodeDrawFcn','MATLAB array');  p.FactoryValue = [];

p = schema.prop(cls,'Nodes','handle vector');        p.FactoryValue = handle([]);
p = schema.prop(cls,'Edges','handle vector');        p.FactoryValue = handle([]);

%% private properties
p = schema.prop(cls,'IsLaidout','bool');             p.FactoryValue = false;
                                                     p.Visible = 'off';
p = schema.prop(cls,'BoundingBox','MATLAB array');   p.FactoryValue = [];
                                                     p.Visible = 'off';
p = schema.prop(cls,'DestroyListener','handle');     p.FactoryValue = [];
                                                     p.Visible = 'off';
p = schema.prop(cls,'to','MATLAB array');            p.FactoryValue = [];
                                                     p.Visible = 'off';
p = schema.prop(cls,'from','MATLAB array');          p.FactoryValue = [];
                                                     p.Visible = 'off';
%% handles to HG 
p = schema.prop(cls,'hgAxes','MATLAB array');        p.FactoryValue = [];
                                                     p.Visible = 'off';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function valueStored = verifyPositiveValue(obj, valueProposed)
if valueProposed>0
    valueStored = valueProposed;
else
    error(message('bioinfo:biographschema:ParameterNotPositive'))
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function valueStored = verifyProperSize(obj, valueProposed)
if rem(valueProposed,1)==0
    valueStored = valueProposed;
else
    error(message('bioinfo:biographschema:ParameterNotInteger'))
end
