%% Visually Representing Interconnected Data
% This example shows how to use the BIOGRAPH object to visually represent
% interconnected data.
%
% The need for representing interconnected data appears in several
% bioinformatics applications. For example, protein-protein interactions,
% network inference, reaction pathways, cluster data, Bayesian networks,
% and phylogenetic trees can be represented with interconnected graphs. The
% BIOGRAPH object allows you to create a comprehensive and graphical layout
% of this type of data. In this example you learn how to populate a
% BIOGRAPH object, render it, and then modify its properties in order to
% customize its display.

%   Copyright 2003-2012 The MathWorks, Inc.


%% Representing a Phylogenetic Tree as a Graph
% Read a phylogenetic tree into a PHYTREE object.

tr = phytreeread('pf00002.tree');

%%
% Reduce the tree to only the human proteins (to make the example smaller,
% you can also use the full tree by omitting the following lines).

sel = getbyname(tr,'human');
tr = prune(tr,~sel(1:33))

%%
% The |plot| method for a PHYTREE object can create a basic layout of the
% phylogenetic tree; however, the graph elements are static. 

plot(tr)

%%
% The PHYTREE object information can be put into a BIOGRAPH object, so you
% can create a dynamic layout. First, pull some information from the PHYTREE
% object. Use |get| to obtain the PHYTREE object properties and the
% |getmatrix| method to obtain the connection matrix.

[names,nn,nb,nl] = get(tr,'NodeNames','NumNodes','NumBranches','NumLeaves');
cm = getmatrix(tr);

%%
% A connection matrix of a low degree graph is best represented by a sparse
% matrix. The average degree of a phylogenetic tree is approximately equal
% to two. You can use the function |spy| to visualize the sparsity pattern,
% every mark represents an edge in the graph.

figure
spy(cm)
colormap(flipud(bone))
title('Top-down direct edges')

%%
% A different approach to visualize this information is to look at pairwise
% distances between all the nodes (branches and leaves) in the tree. Use
% the |pdist| method for PHYTREE objects to find the pairwise distances.

dm = pdist(tr,'criteria','levels','nodes','all','square',true);
figure
imagesc(dm)
colormap(flipud(bone))
axis image
title('Pairwise distances (levels)')

%%
% Phylogenetic tree datasets provide the root of the tree as the last
% element in the set. When building a BIOGRAPH connection matrix, it helps
% to reverse the order of both the data and names, so that the graph is
% built from the root out. This will create a more logical and visually
% appealing presentation.

cm = flipud(fliplr(cm));
names = flipud(names);

%%
% Call the BIOGRAPH object constructor with the connection matrix and the
% node IDs. To explore its properties you can use the |get| function.
bg = biograph(cm,names)
get(bg)

%%
% Once a BIOGRAPH object has been created with the essential information
% (the connection matrix and node IDs), you can modify its properties. For
% example, change the layout type to 'radial', which is best for
% phylogenetic data and the scale.

bg.LayoutType = 'radial'; 
bg.LayoutScale = 3/4; 
get(bg)

%%
% Although nodes and edges have been created, the BIOGRAPH object does not
% have the coordinates in which the graph elements should be drawn such
% that its rendering results into a pretty and uncluttered display.
% Before rendering a BIOGRAPH object, you need to calculate an appropriate
% location for every node. |dolayout| is the method that calls the layout
% engine.

%%
% Some properties of the BIOGRAPH object interact with the layout engine,
% among them 'LayoutType' which selects the layout algorithm. 
dolayout(bg)

%%
% Draw the BIOGRAPH object in a viewer window. The |view| method creates a
% Graphical User Interface (GUI) with the interconnected graph returning
% a handle to a deep copy of the BIOGRAPH object which is contained by
% the figure. With this object handle you can later change some of the
% rendering properties.
bgInViewer = view(bg)

%% Changing the BIOGRAPH Object Properties
% You might want to change the color of all the nodes that represent
% branches. Knowing that the first 'nb' nodes are branches, you can use the 
% vectorized form of |set| to change the 'Color' property of these nodes.
nodeHandlers = bgInViewer.Nodes;
branchHandlers = nodeHandlers(1:nb);
leafHandlers = nodeHandlers(nb+1:end);
set(branchHandlers,'Color',[.7 1 .7])
set(leafHandlers,'Color',[1 .7 .7])

%%
% Changing some geometrical properties requires you to call the |dolayout|
% method again to update the graph to the desired specifications. 
 
% First change the 'Shape' of the branches to circles.
set(branchHandlers,'Shape','circle')

%%
% Notice that the new shape is an ellipse and the edges do not connect
% nicely to the limits of new shapes.

%%
% Now, run the layout engine over the BIOGRAPH object contained by the
% viewer to correct the shapes and the edges.
dolayout(bgInViewer)

%%
% The extent (size) of the nodes is estimated automatically using the
% node 'FontSize' and 'Label' properties. You can force the nodes to have
% any size by turning off the BIOGRAPH 'NodeAutoSize' property and then
% refreshing the layout.
bgInViewer.NodeAutoSize = 'off';
set(branchHandlers,'Size',[20 20])
dolayout(bgInViewer)

%%
% To remove the labels from the branch node, we need to manually copy the
% text strings from the 'ID' property to the 'Label' property. |dolayout|
% automatically sets the 'ShowTextInNodes' property to 'ID' if all nodes
% have their 'Label' property empty. By default when a new biograph object
% is created the 'Label' properties are empty.

for i = 1:numel(leafHandlers)
    leafHandlers(i).Label = leafHandlers(i).ID;
end
bgInViewer.ShowTextInNodes = 'label';

%% Drawing Customized Nodes
% You can draw your own customized nodes in the layout; for example, pie
% charts or histograms may be embedded into the nodes. In this example you
% will use the function |customnodedraw| (an example in the biodemos
% directory) to display the atomic composition of each protein as a pie
% chart. Use this function as a template to create your own customized
% nodes. 

%%
% Get the sequences of the current human proteins you are working with and
% put the sequences into the "UserData" property of their respective nodes.
% Also store the vector with the respective atomic composition.
seqs = fastaread('pf00002.fa','ignoregaps',true)
idxs = seqmatch(get(leafHandlers,'ID'),{seqs.Header});
for i = 1:numel(leafHandlers)
    seq = seqs(idxs(i));
    comp = struct2cell(atomiccomp(seq));
    leafHandlers(i).UserData = seq;
    leafHandlers(i).UserData.Distribution = [comp{:}];
end

%%
% Point the BIOGRAPH object to the customized function that draw nodes. In
% this example |customnodedraw| looks into the 'UserData.Distribution'
% property for the data used in the pie chart.
set(leafHandlers,'Size',[40 40],'shape','circle')
bgInViewer.ShowTextInNodes = 'none';
bgInViewer.CustomNodeDrawFcn = @(node) customnodedraw(node);
bgInViewer.dolayout

%%
% You may attach to nodes additional functionality, such as open links in
% the web browser or perform some calculations, in this case we open the
% aminoacid sequence with |seqviewer|
bgInViewer.NodeCallbacks = {@(x) seqviewer(x.UserData)}

%%
% Place an empty sequence in the branch nodes to avoid an error when the
% callback function looks for the field "Sequence".
set(branchHandlers,'UserData',struct('Sequence','-')) 

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20SEQSTATSDEMO%20in%20Bioinformatics%20Toolbox%204.2 Suggest an enhancement for this example.>*

displayEndOfDemoMessage(mfilename)
