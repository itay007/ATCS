function schema
%

% Copyright 2003-2010 The MathWorks, Inc.

%% add bggui class to package
bgpk = findpackage('biograph');
cls = schema.class(bgpk,'bggui');

%% private properties
p = schema.prop(cls,'biograph','handle');     p.FactoryValue = handle([]);
                                              p.Visible = 'off';
p = schema.prop(cls,'hgFigure','MATLAB array');      p.FactoryValue = [];
                                              p.Visible = 'off';
p = schema.prop(cls,'selecting','bool');      p.FactoryValue = false;
                                              p.Visible = 'off';
p = schema.prop(cls,'dragging','bool');       p.FactoryValue = false;
                                              p.Visible = 'off';
p = schema.prop(cls,'moving','bool');         p.FactoryValue = false;
                                              p.Visible = 'off';
p = schema.prop(cls,'innode','bool');         p.FactoryValue = false;
                                              p.Visible = 'off';
p = schema.prop(cls,'opening','bool');        p.FactoryValue = false;
                                              p.Visible = 'off';
p = schema.prop(cls,'currentNode','handle');  p.FactoryValue = handle([]);
                                              p.Visible = 'off';
p = schema.prop(cls,'dragBox','MATLAB array');       p.FactoryValue = [];
                                                     p.Visible = 'off';
p = schema.prop(cls,'dataTip','MATLAB array');       p.FactoryValue = [];
                                                     p.Visible = 'off';
p = schema.prop(cls,'ghostNodes','MATLAB array');    p.FactoryValue = [];
                                                     p.Visible = 'off';
p = schema.prop(cls,'movNodesPos','MATLAB array');   p.FactoryValue = [0,0];
                                                     p.Visible = 'off';
p = schema.prop(cls,'movPatchXData','MATLAB array'); p.FactoryValue = [];
                                                     p.Visible = 'off';
p = schema.prop(cls,'movPatchYData','MATLAB array'); p.FactoryValue = [];
                                                     p.Visible = 'off';



