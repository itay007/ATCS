function out = graphisdag(G)
%GRAPHISDAG tests for cycles in a directed graph.
%  
% GRAPHISDAG(G) returns true if G is a directed acyclic graph. G is a
% directed graph represented by an n-by-n sparse matrix in which all
% nonzero entries indicate the presence of an edge. 
%
% Examples:
%
%   % Create a DAG and check if it has cycles
%   DG = sparse([1 1 1 2 2 3 4 6],[2 4 6 3 5 4 6 5],true,6,6)
%   view(biograph(DG))
%   graphisdag(DG)
%   % Place a cycle and test it again
%   DG(5,1) = true;
%   view(biograph(DG))
%   graphisdag(DG)
%   
%   % Build a random acyclic graph with 15 nodes and 20 edges
%   g = sparse([],[],true,15,15);
%   while nnz(g) < 20
%      edge = randsample(15*15,1); % get a random edge
%      g(edge) = true;
%      g(edge) = graphisdag(g);
%   end
%   graphisdag(g)
%   view(biograph(g))
%
%   % Check if the gene ontology is a directed acyclic graph
%   GO = geneont('live',true); % may take some time
%   CM = getmatrix(GO);
%   graphisdag(CM)
%
% See also: GRAPHALLSHORTESTPATHS, GRAPHCONNCOMP, GRAPHISOMORPHISM,
% GRAPHISSPANTREE, GRAPHMAXFLOW, GRAPHMINSPANTREE, GRAPHPRED2PATH,
% GRAPHSHORTESTPATH, GRAPHTHEORYDEMO, GRAPHTOPOORDER, GRAPHTRAVERSE.

%   Copyright 2006-2008 The MathWorks, Inc.


out = graphalgs('isdag',0,true,G);
