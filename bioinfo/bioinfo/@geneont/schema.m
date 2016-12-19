function schema
%

% Copyright 2003-2005 The MathWorks, Inc.
schema.package('geneont');

%% set ontology types
if isempty(findtype('OntologyTypes'))
    schema.EnumType('OntologyTypes',{'gene ontology','cellular component','molecular function','biological process','unknown'});
end
