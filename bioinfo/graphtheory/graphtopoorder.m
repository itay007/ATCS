function order = graphtopoorder(G)
%GRAPHTOPOORDER performs topological sort of a directed acyclic graph.
%  
% ORDER = GRAPHTOPOORDER(G) returns an index vector with the order of the
% nodes sorted topologically. In a topological order there can only exist
% edges between node u to node v, if and only if u appears first in ORDER
% than v. G is a directed acyclic graph represented by an n-by-n sparse
% matrix in which all nonzero entries indicate the presence of an edge. 
%
% Example:
%   % Create a DAG with 6 nodes and 8 edges
%   DG = sparse([6 6 6 2 2 3 5 1],[2 5 1 3 4 5 1 4],true,6,6)
%   view(biograph(DG))
%   % Find the topological order
%   order = graphtopoorder(DG)
%   % Permute the nodes, so they appear ordered when showing the graph
%   DG = DG(order,order)
%   view(biograph(DG))
%
% See also: GRAPHALLSHORTESTPATHS, GRAPHCONNCOMP, GRAPHISDAG,
% GRAPHISOMORPHISM, GRAPHISSPANTREE, GRAPHMAXFLOW, GRAPHMINSPANTREE,
% GRAPHPRED2PATH, GRAPHSHORTESTPATH, GRAPHTHEORYDEMO, GRAPHTRAVERSE.

%   Copyright 2006 The MathWorks, Inc.


 order = graphalgs('topoorder',0,true,G);
