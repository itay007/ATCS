%% Connecting to the KEGG API Web Service
% This example shows how to access the KEGG system using the SOAP/WSDL
% based web service from within MATLAB(R). We show some examples using the
% KEGG Pathway API. Refer to
% <http://www.genome.jp/kegg/soap/doc/keggapi_manual.html the KEGG API
% Manual> for more examples and information.
%
% Note: The pathway information and maps shown in this example may differ
% from the results you get due to the frequent update to the KEGG pathway
% database.

%   Copyright 2006-2012 The MathWorks, Inc.


%% Setting up Web Connection
% Create a variable containing the definition URL of the Web Service
% Description Language (WSDL) web service for the KEGG database
wsdlURL = 'http://soap.genome.jp/KEGG.wsdl';

%%
% In MATLAB, use the |createClassFromWSDL| function to call web service
% methods. Create the classes from KEGG WSDL. This will create a directory
% called @KEGG in the current directory.
className = createClassFromWsdl(wsdlURL);

%%
% The @KEGG directory contains automatically generated files that implement
% the KEGG web service methods
dir @KEGG
%%
% You can also use the |methods| command to view the list of methods
% available for the KEGG class. You will notice that there are more methods
% available than files in the @KEGG directory. These are inherited methods
% that are available for all objects in MATLAB. 
methods(KEGG)

%% Using the Web Service
% In order to use the web service, you must first create an instance of the
% KEGG object. 
kegg = KEGG;
%%
% You can confirm that |kegg| is an instance of the KEGG using the |class|
% command
classType = class(kegg)

%% Using KEGG Pathway API
% The KEGG API functions work with KEGG IDs only. The |bconv| method
% converts NCBI GIs, NCBI GeneIDs, GenBank(R) IDs, and UniProt IDs to KEGG
% IDs
kegg_ids_conv = bconv(kegg, 'ncbi-gi:10047086 ncbi-geneid:14751')
%%
% Use the |list_organisms| method to get a structure array containing the
% organisms in the KEGG/GENES database. The organism full names are in the
% definition field.
organisms = list_organisms(kegg)
%%
% Find the entry where definition has the string _Homo sapiens_ .
homo_idx = find(~cellfun(@isempty,regexpi({organisms.definition},'Homo sapiens')))
homo = organisms(homo_idx)

%%
% You can get a list of pathway ids for _Homo sapiens_ in KEGG PATHWAY
% database.
pathway_list = list_pathways(kegg, homo.entry_id);
num_pathways = length(pathway_list)
%%
% This example uses the Fatty acid metabolism pathway. Find the entry where
% definition has the string _Fatty acid metabolism_.
famp_idx =  find(~cellfun(@isempty,regexpi({pathway_list.definition},'Fatty acid metabolism')))
famp_pathway = pathway_list(famp_idx)

%%
% Get lists of genes, compounds, enzymes and reactions involved in the
% Fatty acid metabolism pathway.
famp_genes = get_genes_by_pathway(kegg, famp_pathway.entry_id);
num_genes = length(famp_genes)
%%
famp_compounds= get_compounds_by_pathway(kegg, famp_pathway.entry_id);
num_compounds = length(famp_compounds)
%%
famp_enzymes = get_enzymes_by_pathway(kegg, famp_pathway.entry_id);
num_enzymes = length(famp_enzymes)
%%
famp_reactions = get_reactions_by_pathway(kegg, famp_pathway.entry_id);
num_reactions = length(famp_reactions)

%%
% You can get KEGG pathway IDs from a list of KEGG gene IDs, enzymes or
% reactions. For example, you can get the pathways IDs for the first 10
% genes and reactions.
pathway_by_genes = get_pathways_by_genes(kegg, famp_genes(1:10))

%%
pathway_by_react = get_pathways_by_reactions(kegg, famp_reactions(1:10))

%% Coloring Pathways
% In KEGG pathway maps, a gene or enzyme is represented by a rectangle, and
% a compound is shown as a small circle. In this example, the Fatty acid
% metabolism pathway map returned by KEGG has already colored the gene
% products related to _Homo sapiens_ in green. In addition, you can specify
% colors for specific genes or compounds, for example, the highly expressed
% genes from a gene expression experiment. KEGG returns an URL of the given
% pathway map with the elements corresponding to the specified colored
% genes.
%%
% List the objects to be colored, and specify the colors for each object
obj_list = [famp_genes(30:31); {'hsa:8310'}; {'cpd:C00823'}; famp_compounds(5:6)];
fg_list  = {'#ff6600', '#0000ff', '#0000ff', '#0000ff', '#ff0000', '#ff0000'};
bg_list  = {'#99ccff', 'yellow',  '#ff6633', '#ff0000', '#ccffff', '#ccffff'};

%%
% KEGG colors the given objects on the pathway map with the specified
% colors and return the URL of a static colored pathway map in _gif_ format. 
pathway_map_colored = color_pathway_by_objects(kegg, famp_pathway.entry_id,...
                        obj_list, fg_list, bg_list)
%%
% The |get_html_of_colored_pathway_by_objects| method returns an URL of an
% interactive pathway map.
pathway_map_html = get_html_of_colored_pathway_by_objects(kegg,...
                    famp_pathway.entry_id, obj_list, fg_list, bg_list)
                
%% Displaying Pathway Maps
%%
% You can simply display the pathway map in a browser.
web(pathway_map_html)

%%
% You can display the static pathway map in a figure. Use |imread| function
% to read the image data and the associate color map from the pathway map
% in _gif_ format.
[x,cmap] = imread(pathway_map_colored);

%%
% Display the image in a figure. You need to set the color map of the
% figure to cmap. Note: If you have Image Processing Toolbox(TM), just use
% |imshow(x,map)| to display the pathway map.
hfig = figure('Colormap', cmap);
hax = axes('Parent', hfig);
himg = image(x, 'Parent', hax);
set(hax, 'Visible', 'off')
scaleimagefigure(hfig, hax, himg);

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20CONNECTKEGG%20in%20Bioinformatics%20Toolbox%204.2 Suggest an enhancement for this example.>*

displayEndOfDemoMessage(mfilename)
