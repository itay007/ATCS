function out = graphisspantree(G)
%GRAPHISSPANTREE test a spanning tree.
%  
% GRAPHISSPANTREE(G) returns true if G is spanning tree. A spanning tree
% must touch all the nodes and must be acyclic. G is an undirected graph
% represented by the lower triangle of an n-by-n sparse matrix in which all
% nonzero entries indicate the presence of an edge. 
%
% Example:
%   % Create a connection matrix from a phylogenetic tree
%   tr = phytreeread('pf00002.tree')
%   [CM,labels,dist] = getmatrix(tr);
%   % check if CM is a spanning tree
%   graphisspantree(CM)
%   % Add an edge between the root and the first leaf
%   CM(end,1) = 1;
%   % check if CM is a spanning tree
%   graphisspantree(CM)
%
% See also: GRAPHALLSHORTESTPATHS, GRAPHCONNCOMP, GRAPHISDAG,
% GRAPHISOMORPHISM, GRAPHMAXFLOW, GRAPHMINSPANTREE, GRAPHPRED2PATH,
% GRAPHSHORTESTPATH, GRAPHTHEORYDEMO, GRAPHTOPOORDER, GRAPHTRAVERSE.

%   Copyright 2006-2008 The MathWorks, Inc.


out = graphalgs('istree',0,false,G);
