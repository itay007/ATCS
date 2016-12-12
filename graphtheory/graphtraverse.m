function [disc,pred,close,depthCount] = graphtraverse(G,S,varargin)
%GRAPHTRAVERSE performs traversal of the graph by following adjacent nodes.
% 
% [DISC,PRED,CLOSE] = GRAPHTRAVERSE(G,S) traverses graph G starting from
% the node S. G is a directed graph represented by an n-by-n sparse matrix
% in which all nonzero entries indicate the presence of an edge. DISC is a
% list of nodes in the order in which they are discovered. PRED is a vector
% of predecessor node indices of the resulting spanning tree, and, CLOSE is
% a list of node indices in the order in which they are closed. 
% 
% GRAPHTRAVERSE(...,'METHOD',METHOD) selects the algorithm to use, options
% are:
%    'BFS'         - Breadth First Search.
%   ['DFS']        - Depth First Search. 
% 
% The time complexity of either of these algorithms is O(n+e), where n and
% e are number of nodes and edges respectively.
% 
% GRAPHTRAVERSE(...,'DIRECTED',false) indicates that the graph G is
% undirected, upper triangle of the sparse matrix is ignored. Default is
% true.
% 
% GRAPHTRAVERSE(...,'DEPTH',D) limits the depth of the search. Default is
% Inf.
%  
% Examples:
%
%   % Create a directed graph with 10 nodes and 12 edges
%   DG = sparse([1 2 3 4 5 5 5 6 7 8 8 9],[2 4 1 5 3 6 7 9 8 1 10 2],true,10,10)
%
%   % Find the DFS discovery order starting at node 4
%   order = graphtraverse(DG,4)
%   % Label the nodes with the DFS discovery order
%   h = view(biograph(DG))
%   for i = 1:10
%       h.Nodes(order(i)).Label = sprintf('%s:%d',h.Nodes(order(i)).ID,i);
%   end
%   h.ShowTextInNodes = 'label'
%   dolayout(h)
%
%   % Find the BFS discovery order starting at node 4
%   order = graphtraverse(DG,4,'Method','BFS')
%   % Label the nodes with the BFS discovery order
%   h = view(biograph(DG))
%   for i = 1:10
%       h.Nodes(order(i)).Label = sprintf('%s:%d',h.Nodes(order(i)).ID,i);
%   end
%   h.ShowTextInNodes = 'label'
%   dolayout(h)
%
%   % Find and mark nodes that are close (2 edges) from node 4
%   h = view(biograph(DG))
%   node_idxs = graphtraverse(DG,4,'depth',2)
%   set(h.nodes(node_idxs),'Color',[1 0 0])
%
%
% See also: GRAPHALLSHORTESTPATHS, GRAPHCONNCOMP, GRAPHISDAG,
% GRAPHISOMORPHISM, GRAPHISSPANTREE, GRAPHMAXFLOW, GRAPHMINSPANTREE,
% GRAPHPRED2PATH, GRAPHSHORTESTPATH, GRAPHTHEORYDEMO, GRAPHTOPOORDER.
%
% Reference: 
%  [1]	R. Sedgewick "Algorithms in C++, Part 5 Graph Algorithms"
%       Addison-Wesley, 2002.  

%   Copyright 2006-2008 The MathWorks, Inc.


algorithms = {'bfs','dfs'};
debug_level = 0;

% set defaults of optional input arguments
algorithm  = 2;
directed = true;
depth = inf;

% read in optional PV input arguments
nvarargin = numel(varargin);
if nvarargin
    if rem(nvarargin,2) == 1
        error(message('bioinfo:graphtraverse:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'method','directed','depth'};
    for j=1:2:nvarargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:graphtraverse:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:graphtraverse:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % 'method'
                    algorithm = find(strncmpi(pval,algorithms,numel(pval))); 
                    if isempty(algorithm)
                        error(message('bioinfo:graphtraverse:NotValidMethod', pval))
                    elseif numel(algorithm)>1
                         error(message('bioinfo:graphtraverse:AmbiguousMethod', pval))
                    end
                case 2 % 'directed'
                    directed = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 3 % 'depth'
                    depth = pval;
            end
        end
    end
end

% call the mex implementation of the graph algorithms
if nargout>3 % fourth output is undocummneted now and may change in the future
    [disc,pred,close,depthCount] = graphalgs(algorithms{algorithm},debug_level,directed,G,S,depth);
elseif nargout>2
    [disc,pred,close] = graphalgs(algorithms{algorithm},debug_level,directed,G,S,depth);
elseif nargout>1
    [disc,pred] = graphalgs(algorithms{algorithm},debug_level,directed,G,S,depth);
else
    disc = graphalgs(algorithms{algorithm},debug_level,directed,G,S,depth);
end
