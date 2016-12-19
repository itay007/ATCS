function schema
%

% Copyright 2003-2012 The MathWorks, Inc.

pkg = findpackage('geneont');
cls = schema.class(pkg,'geneont');

p = schema.prop(cls,'default_namespace','string');    p.FactoryValue =  '';
p = schema.prop(cls,'format_version','string');       p.FactoryValue =  '';
p = schema.prop(cls,'data_version','string');         p.FactoryValue =  '';
p = schema.prop(cls,'version','string');              p.FactoryValue =  '';
p = schema.prop(cls,'date','string');                 p.FactoryValue =  '';
p = schema.prop(cls,'saved_by','string');             p.FactoryValue =  '';
p = schema.prop(cls,'auto_generated_by','string');    p.FactoryValue =  '';
p = schema.prop(cls,'subsetdef','MATLAB array');      p.FactoryValue =  [];
p = schema.prop(cls,'import','string');               p.FactoryValue =  '';
p = schema.prop(cls,'synonymtypedef','string');       p.FactoryValue =  '';
p = schema.prop(cls,'idspace','string');              p.FactoryValue =  '';
p = schema.prop(cls,'default_relationship_id_prefix','string'); 
                                                      p.FactoryValue =  '';
p = schema.prop(cls,'id_mapping','string');           p.FactoryValue =  '';
p = schema.prop(cls,'remark','string');               p.FactoryValue =  '';
p = schema.prop(cls,'typeref','string');              p.FactoryValue =  '';
p = schema.prop(cls,'unrecognized_tag','MATLAB array');      
                                                      p.FactoryValue =  [];
p = schema.prop(cls,'Terms','handle vector');         p.FactoryValue = handle([]);

%% private properties
p = schema.prop(cls,'hashID','MATLAB array');       p.FactoryValue = [];
                                                    p.Visible = 'off';
p = schema.prop(cls,'relationmatrix','MATLAB array'); p.FactoryValue = [];
                                                    p.Visible = 'off';
                                                    p.GetFunction = @getrelationmatrix;
                                                    p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'relationtype','MATLAB array'); p.FactoryValue = [];
                                                    p.Visible = 'off';
p = schema.prop(cls,'gonumlist','MATLAB array');    p.FactoryValue =  [];   
                                                    p.Visible = 'off';                                                    
p = schema.prop(cls,'isarelationmatrix','MATLAB array'); p.FactoryValue = [];
                                                    p.Visible = 'off';
p = schema.prop(cls,'partofrelationmatrix','MATLAB array'); p.FactoryValue = [];
                                                    p.Visible = 'off';
p = schema.prop(cls,'isarelationmatrixi','MATLAB array'); p.FactoryValue = [];
                                                    p.Visible = 'off';
p = schema.prop(cls,'partofrelationmatrixi','MATLAB array'); p.FactoryValue = [];
                                                    p.Visible = 'off';                                                    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = getrelationmatrix(h,w) %#ok
% In Bioinfo 3.0 RELATIONMATRIX is substituted by ISARELATIONMATRIX and
% PARTOFRELATIONMATRIX, we tie this property to this function for backwards
% compatibility
v = (2*get(h,'partofrelationmatrix')+sparse(1,1,0))+get(h,'isarelationmatrix');
