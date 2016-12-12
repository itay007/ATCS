function [m,flow,cuts] = graphmaxflow(G,S,D,varargin)
%GRAPHMAXFLOW calculates the maximum flow in a directed graph.
% 
% [M,F,K] = GRAPHMAXFLOW(G,S,D) Calculates the maximum flow of the directed
% graph G from node S to node D. G is an n-by-n sparse matrix that
% represents a directed graph with non reciprocal edges. Nonzero entries in
% G determine the capacity of the edges. M is the maximum flow and F is a
% sparse matrix with all the flow values for every edge. F(i,j) is the flow
% from node i to j. K is a logical row vector indicating the nodes
% connected to S after calculating the minimum between S and D. If several
% solutions to the minimum cut problem exist then K is a matrix.
% 
% GRAPHMAXFLOW(...,'METHOD',METHOD) selects the algorithm to use, options
% are: 
%    'Edmonds'    - Edmonds and Karp algorithm. Implementation is based on
%                   a variation called the "labeling algorithm". Time
%                   complexity is O(n*e^2). 
%   ['Goldberg']  - Goldberg algorithm uses the generic method known as
%                   "preflow-push". Time complexity is O(n^2*sqrt(e)).
% 
% Notes: (1) n and e are number of nodes and edges respectively. (2) The
% algorithm that finds all minimum cuts has a time complexity of O(2^n),
% avoid using the third output argument if it is not required.
% 
% GRAPHMAXFLOW(...,'CAPACITY',C) provides custom capacities for the edges.
% C is a column vector with one entry for every edge in G, traversed
% column-wise. 
%
% Example:
%   % Create a directed graph with 6 nodes
%   cm = sparse([1 1 2 2 3 3 4 5],[2 3 4 5 4 5 6 6],[2 3 3 1 1 1 2 3],6,6)
%   % Call the maximum flow algorithm between 1 and 6
%   [M,F,K] = graphmaxflow(cm,1,6)
%   % View graph with original flows
%   h = view(biograph(cm,[],'ShowWeights','on'))
%   % View graph with actual flows
%   view(biograph(F,[],'ShowWeights','on'))
%   % Show in the original graph one solution of the mincut problem
%   set(h.Nodes(K(1,:)),'Color',[1 0 0])
%
% See also: GRAPHALLSHORTESTPATHS, GRAPHCONNCOMP, GRAPHISDAG,
% GRAPHISOMORPHISM, GRAPHISSPANTREE, GRAPHMINSPANTREE, GRAPHPRED2PATH,
% GRAPHSHORTESTPATH, GRAPHTHEORYDEMO, GRAPHTOPOORDER, GRAPHTRAVERSE.
%
% References: 
%  [1]	J. Edmonds and R.M. Karp "Theoretical improvements in the
%       algorithmic efficiency for network flow problems" Journal of the
%       ACM, 19:248-264, 1972. 
%  [2]	A.V. Goldberg "A New Max-Flow Algorithm" MIT Technical report
%       MIT/LCS/TM-291, 1985. 

%   Copyright 2006-2008 The MathWorks, Inc.


algorithms = {'edmonds','goldberg'};
algorithmkeys = {'edm','gol'};
debug_level = 0;

% set defaults of optional input arguments
C = []; % no custom capacities
algorithm  = 2; % defaults to Goldberg

% read in optional PV input arguments
nvarargin = numel(varargin);
if nvarargin
    if rem(nvarargin,2) == 1
        error(message('bioinfo:graphmaxflow:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'method','capacity'};
    for j=1:2:nvarargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:graphmaxflow:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:graphmaxflow:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % 'method'
                    algorithm = find(strncmpi(pval,algorithms,numel(pval))); 
                    if isempty(algorithm) 
                        error(message('bioinfo:graphmaxflow:NotValidMethod', pval))
                    elseif numel(algorithm)>1
                         error(message('bioinfo:graphmaxflow:AmbiguousMethod', pval))
                    end
                case 2 % 'capacities'
                    C = pval(:);
            end
        end
    end
end

% call the mex implementation of the graph algorithms
if nargout>1
    if isempty(C)
        [m,flow] = graphalgs(algorithmkeys{algorithm},debug_level,true,G,S,D);
    else
        [m,flow] = graphalgs(algorithmkeys{algorithm},debug_level,true,G,S,D,C);
    end    
else
    if isempty(C)
        m = graphalgs(algorithmkeys{algorithm},debug_level,true,G,S,D);
    else
        m = graphalgs(algorithmkeys{algorithm},debug_level,true,G,S,D,C);
    end
end

% find the mincut solutions
if nargout>2   
    n = length(G);
    cuts = false(0,n);
    mask = true(1,n);
    mask([S D]) = false;
    cut = false(1,n);
    cut(S) = true;
    idx = true(1,n-2);
    h = 0;
    if ~isempty(C)
        if islogical(G)
            G = double(G);
        end
        G(G>0)=C;
    end
    qq = nonzeros(G);
    n_times_eps = n*eps*max(qq(~isinf(qq)));
    counter = 10000000;
    ws = warning('query','bioinfo:graphmaxflow:mincutsLimit');
    conterDecrease = strcmp(ws.state,'on');
    while ~isempty(h) && counter
       cut(mask)=idx;
       if (sum(nonzeros(G(cut,~cut)))-m < n_times_eps)
           cuts = [cuts;cut]; %#ok<AGROW>
       end
       h = find(idx,1);
       idx(1:h-1) = true;
       idx(h) = false;
       counter = counter - conterDecrease;
    end
    if ~counter
        warning(message('bioinfo:graphmaxflow:mincutsLimit'))
    end
end
