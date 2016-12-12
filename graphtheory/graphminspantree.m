function [T,pred] = graphminspantree(G,varargin)
%GRAPHMINSPANTREE finds the minimal spanning tree in graph.
% 
% [T, PRED] = GRAPHMINSPANTREE(G) finds an acyclic subset of edges that
% connects all the nodes in the undirected graph G and for which the total
% weight is minimized. Weights of the edges are all nonzero entries in the
% lower triangle of the n-by-n sparse matrix G. T is a spanning tree
% represented by a sparse matrix. The output PRED contains the predecessor
% nodes of the minimal spanning tree with the root node indicated by a
% zero. The root defaults to the first node in the largest connected
% component, which requires an extra call to the graphconncomp function.
% 
% [T, PRED] = GRAPHMINSPANTREE(G,R) sets the root of the minimal spanning
% tree to node R.
%
% GRAPHMINSPANTREE(...,'METHOD',METHOD) selects the algorithm to use,
% options are: 
%    ['Prim']     - Prim's algorithm grows the MST one edge at a time by
%                   adding a minimal edge that connects a node in the
%                   growing MST with any other node. Time complexity is
%                   O(e*log(n)). 
%    'Kruskal'    - Kruskal's algorithm grows the MST one edge at a time by
%                   finding an edge that connects two trees in a spreading
%                   forest of growing MSTs. Time complexity is
%                   O(e+x*log(n)) where x is the number of edges no longer
%                   than the longest edge in the MST. 
% 
% Note: n and e are number of nodes and edges respectively.
% 
% GRAPHMINSPANTREE(...,'WEIGHTS',W) provides custom weights for the edges,
% useful to indicate zero valued weights. W is a column vector with one
% entry for every edge in G, traversed column-wise.
% 
% Remarks: When the graph is unconnected, Prim's algorithm only returns the
% tree that contains R, while Kruskal's algorithm returns an MST for every
% component. 
% 
% Example:
%   % Create an undirected graph with 6 nodes
%   W = [.41 .29 .51 .32 .50 .45 .38 .32 .36 .29 .21];
%   DG = sparse([1 1 2 2 3 4 4 5 5 6 6],[2 6 3 5 4 1 6 3 4 2 5],W)
%   UG = tril(DG + DG')
%   view(biograph(UG,[],'ShowArrows','off','ShowWeights','on'))
%   % Find the minimum spanning tree of UG
%   [ST,pred] = graphminspantree(UG)
%   view(biograph(ST,[],'ShowArrows','off','ShowWeights','on'))
%
% See also: GRAPHALLSHORTESTPATHS, GRAPHCONNCOMP, GRAPHISDAG,
% GRAPHISOMORPHISM, GRAPHISSPANTREE, GRAPHMAXFLOW, GRAPHPRED2PATH,
% GRAPHSHORTESTPATH, GRAPHTHEORYDEMO, GRAPHTOPOORDER, GRAPHTRAVERSE.
%
% References: 
%  [1]	J. B. Kruskal. "On the shortest spanning subtree of a graph and the
%       traveling salesman problem" In Proceedings of the American
%       Mathematical Society, 7:48-50, 1956. 
%  [2]	R. Prim. "Shortest connection networks and some generalizations"
%       Bell System Technical Journal, 36:1389-1401, 1957.

%   Copyright 2006-2008 The MathWorks, Inc.


algorithms = {'prim','kruskal'};
algorithmkeys = {'pri','kru'};
debug_level = 0;

% set defaults of optional input arguments
W = []; % no custom weights
R = []; % no root given
algorithm  = 1; % defaults to prim

% find out signature of input arguments
if nargin>1 && isnumeric(varargin{1})
    R = varargin{1};
    varargin(1) = [];
end

% read in optional PV input arguments
nvarargin = numel(varargin);
if nvarargin
    if rem(nvarargin,2) == 1
        error(message('bioinfo:graphminspantree:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'method','weights'};
    for j=1:2:nvarargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:graphminspantree:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:graphminspantree:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % 'method'
                    algorithm = find(strncmpi(pval,algorithms,numel(pval))); 
                    if isempty(algorithm) 
                        error(message('bioinfo:graphminspantree:NotValidMethod', pval))
                    elseif numel(algorithm)>1
                         error(message('bioinfo:graphminspantree:AmbiguousMethod', pval))
                    end
                case 2 % 'weights'
                    W = pval(:);
            end
        end
    end
end

% find manually the best root (if it was not given)
if isempty(R)
    [num_comp,classes] = graphconncomp(G,'directed',false);
    if num_comp==1
        R = 1;
    else
        R = find(classes==mode(classes),1,'first');
    end
end

% call the mex implementation of the graph algorithms
if nargout>1
    if isempty(W)
        [T,pred] = graphalgs(algorithmkeys{algorithm},debug_level,false,G,R);
    else
        [T,pred] = graphalgs(algorithmkeys{algorithm},debug_level,false,G,R,W);
    end    
else
    if isempty(W)
        T = graphalgs(algorithmkeys{algorithm},debug_level,false,G,R);
    else
        T = graphalgs(algorithmkeys{algorithm},debug_level,false,G,R,W);
    end 
end
