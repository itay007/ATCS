%% Working with Graph Theory Functions
% This example shows how Bioinformatics Toolbox(TM) can be used to work
% with and vizualize graphs.
%
% Graphs, in the sense of graph theory, are a mathematical way of
% representing connections or relationships between objects. There are many
% applications in bioinformatics where understanding relationships between
% objects is very important. Such applications include phylogenetic
% analysis, protein-protein interactions, pathway analysis and many more. 
% Bioinformatics Toolbox provides a set of generic functions for working
% with and visualizing graphs. 

%  Copyright 2007-2012 The MathWorks, Inc.


%% Creating a Graph from a SimBiology(R) Model
% The graph theory functions in Bioinformatics Toolbox work on sparse
% matrices. The only restriction is that the matrix be square. In this
% example, a graph was created from a SimBiology(R) model of a
% Repressilator [1] oscillatory network . In this model, protein A
% represses protein B, protein B represses protein C, which in turn
% represses protein A.  

load oscillatorgraph

%%
% There are two variables: |g|, a sparse matrix, and |names|, a list of the
% names of the nodes of the graph. 
whos g names

%%
% If you have SimBiology you can create the graph using the following
% commands:

% sbioloadproject oscillator
% class(m1) 
% Now get the adjacency matrix
% [g names] = getadjacencymatrix(m1);


%% Visualizing the Graph
% There are many functions in MATLAB(R) for working with sparse matrices.
% The |spy| function displays as * wherever there is a non-zero element of
% the matrix. 

spy(g)

%%
% This gives some indication of the number of edges of the graph and also
% shows that the graph is not symmetric and, hence, is a directed graph.
% However, it is difficult to visualize what is going on. The |biograph|
% object is another way of representing a graph in Bioinformatics Toolbox.

gObj = biograph(g,names)

%%
% The |view| method lays out the graph and displays it in a figure. 
gObj = view(gObj);

%%
% You can interact with the graph using the mouse. You can also
% programmatically modify the way that the graph is displayed.

% find the nodes pA, pB, and pC
pANode = find(strcmp('pA', names));
pBNode = find(strcmp('pB',names));
pCNode = find(strcmp('pC', names));
% Color these red, green, and blue
set(gObj.nodes(pANode),'Color',[1 0 0],'size',[40 30]);
set(gObj.nodes(pBNode),'Color',[0 1 0],'size',[40 30]);
set(gObj.nodes(pCNode),'Color',[0 0 1],'size',[40 30]);
dolayout(gObj);

%% Using the Graph Theory Functions
% There are several functions in Bioinformatics Toolbox for working with
% graphs. These include |graphshortestpath|, which finds the shortest path
% between two nodes, |graphisspantree|, which checks if a graph is a
% spanning tree, and |graphisdag|, which checks if a graph is a directed
% acyclic graph. 

graphisdag(g)

%%
% There are also corresponding methods of the biograph object. These have
% names similar to the functions for working with sparse matrices but
% without the prefix 'graph'.
isdag(gObj)

%% Finding the Shortest Path Between Nodes pA and pC
% A common question to ask about a graph is what is the shortest path
% between two nodes. Note that in this example all the edges have length 1.

[dist,path,pred] = shortestpath(gObj,pANode,pCNode);

%%
% Color the nodes and edges of the shortest path
set(gObj.Nodes(path),'Color',[1 0.4 0.4])
edges = getedgesbynodeid(gObj,get(gObj.Nodes(path),'ID'));
set(edges,'LineColor',[1 0 0])
set(edges,'LineWidth',1.5)

%%
% You can use |allshortestpaths| to calculate the shortest paths from each
% node to all other nodes. 
allShortest = allshortestpaths(gObj);

%%
% A heatmap of these distances shows some interesting patterns.
imagesc(allShortest)
colormap(pink);
colorbar
title('All Shortest Paths for Oscillator Model');

%% Traversing the Graph
% Another common problem with graphs is finding an efficient way to
% traverse a graph by moving between adjacent nodes. The |traverse| method
% uses a depth-first search by default but you can also choose to use a
% breadth-first search.

order = traverse(gObj,pANode);

%%
% The return value |order| shows the order in which the nodes were
% traversed starting at pA. You can use this to find an alternative path
% from pA to pC. 

alternatePath = order(1:find(order == pCNode));
set(gObj.Nodes(alternatePath),'Color',[0.4 0.4 1])
edges = getedgesbynodeid(gObj,get(gObj.Nodes(alternatePath),'ID'));
set(edges,'LineColor',[0 0 1])
set(edges,'LineWidth',1.5)

%% Finding Connected Components in the Graph
% The oscillator model is cyclic with pA, pB, and pC all connected. The
% method |conncomp| finds connected components.  A strongly connected
% component of a graph is a maximal group of nodes that are mutually
% reachable without violating the edge directions. You can use the
% |conncomp| method to determine which nodes are not part of the main
% cycle.

[S,C] = conncomp(gObj);

% Mark the nodes for each component with different color
colors = flipud(jet(S));
for i = 1:numel(gObj.nodes)
    gObj.Nodes(i).Color = colors(C(i),:);
end

%%
% You will notice that the "trash" node is a sink. Several nodes connect to
% this node but there is no path from "trash" to any other node.

%% Simulating Knocking Out a Reaction
% In biological pathways it is common to find that while some reactions are
% essential to the survival of the behavior of the pathway, others are not.
% You can use the sparse graph representation of the pathway to investigate
% whether Reaction1 and Reaction2 in the model are essential to the
% survival of the oscillatory properties.

%%
% Find the nodes in which you are interested.
r1Node= find(strcmp( 'Reaction1', names));
r2Node= find(strcmp( 'Reaction2', names));

%%
% Create copies of the sparse matrix and remove all edges associated with
% the reactions.
gNoR1 = g;
gNoR1(r1Node,:) = 0; gNoR1(:,r1Node)=0;
gNoR2 = g;
gNoR2(r2Node,:) = 0; gNoR2(:,r2Node)=0;

%%
% In the case where we remove Reaction2, there are still paths from pA to
% pC and back and the structure has not changed very much.

distNoR2CA = graphshortestpath(gNoR2,pCNode,pANode)
distNoR2AC = graphshortestpath(gNoR2,pANode,pCNode)
% Display the graph from which Reaction2 was removed.
gNoR2Obj = view(biograph(gNoR2,names));
[S,C] = conncomp(gNoR2Obj);
% Mark the nodes for each component with different color
colors = flipud(jet(S));
for i = 1:numel(gNoR2Obj.nodes)
    gNoR2Obj.Nodes(i).Color = colors(C(i),:);
end
%%
% However, in the case where we remove Reaction1, there is no longer a path
% from pC back to pA.
distNoR1AC = graphshortestpath(gNoR1,pANode,pCNode)
distNoR1CA = graphshortestpath(gNoR1,pCNode,pANode)


%%
% When you visualize the graph from which Reaction1 was removed you will see a
% significant change in the structure of the graph.

% Display the graph from which Reaction1 was removed.
gNoR1Obj = view(biograph(gNoR1,names));
[S,C] = conncomp(gNoR1Obj);
% Mark the nodes for each component with different color
colors = flipud(jet(S));
for i = 1:numel(gNoR1Obj.nodes)
    gNoR1Obj.Nodes(i).Color = colors(C(i),:);
end

%% References
% [1] Elowitz, M.B, and Leibler, S.(2000) A Synthetic Oscillatory Network
%     of Transcriptional Regulators; Nature. 403(6767):335-8. 

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Enhancement%20for%20GRAPHTHEORYDEMO%20in%20Bioinformatics%20Toolbox%204.2 Suggest an enhancement for this example.>*

displayEndOfDemoMessage(mfilename)
