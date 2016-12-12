function schema
%

% Copyright 2003-2006 The MathWorks, Inc.

pkg = findpackage('geneont');
cls = schema.class(pkg,'term');

%% public properties
p = schema.prop(cls,'id','double');              p.FactoryValue =  0;
p = schema.prop(cls,'name','string');            p.FactoryValue =  '';       
p = schema.prop(cls,'ontology','OntologyTypes'); p.FactoryValue = 'gene ontology';
p = schema.prop(cls,'definition','string');      p.FactoryValue =  '';
p = schema.prop(cls,'comment','string');         p.FactoryValue =  '';
p = schema.prop(cls,'synonym','MATLAB array');   p.FactoryValue =  [];  
p = schema.prop(cls,'is_a','MATLAB array');      p.FactoryValue = [];
p = schema.prop(cls,'part_of','MATLAB array');   p.FactoryValue = [];
p = schema.prop(cls,'is_a_i','MATLAB array');    p.FactoryValue = [];
                                                 p.Visible = 'off';
p = schema.prop(cls,'part_of_i','MATLAB array'); p.FactoryValue = [];
                                                 p.Visible = 'off';
p = schema.prop(cls,'obsolete','bool');          p.FactoryValue =  false;   
