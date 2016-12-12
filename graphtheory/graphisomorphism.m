function [F,M]=graphisomorphism(G1,G2,varargin)
%GRAPHISOMORPHISM finds an isomorphism between two graphs.
% 
% [F,M] = GRAPHISOMORPHISM(G1,G2) returns true in F if G1 and G2 are
% isomorphic graphs. A graph isomorphism is a 1-to-1 mapping of the nodes
% in the graph G1 and the nodes in the graph G2 such that adjacency is
% preserved. G1 and G2 are both directed graphs represented by n-by-n
% sparse matrices. F is Boolean, when F is true M contains the node mapping
% from G1 to G2 for one possible isomorphism, and when F is false M is
% empty. The worst-case time complexity is O(n!) where n is the number of
% nodes.
% 
% GRAPHISOMORPHISM(...,'DIRECTED',false) assumes both G1 and G2 are
% undirected graphs, the upper triangles of the sparse matrices G1 and G2
% are ignored. 
% 
% Examples:
%    % Create a directed graph with 8 nodes
%    m('ABCDEFGH') = [1 2 3 4 5 6 7 8];
%    g1 = sparse(m('ABDCDCGEFFG'),m('BCBDGEEFHGH'),true,8,8);
%    view(biograph(g1,'ABCDEFGH'))
%    % Set a random permutation and create a new permuted graph
%    p = randperm(8)
%    g2 = g1(p,p);
%    view(biograph(g2,'12345678'))
%    % Check if both graphs are isomorphic
%    [F,M] = graphisomorphism(g2,g1)
%    % Reverse the D-G edge and check again
%    g1(m('DG'),m('GD')) = g1(m('GD'),m('DG'));
%    view(biograph(g1,'ABCDEFGH'))
%    [F,M] = graphisomorphism(g2,g1)
%    % If undirected graphs are considered they should be isomorphic
%    [F,M] = graphisomorphism(g2+g2',g1+g1','directed',false)
%
% See also: GRAPHALLSHORTESTPATHS, GRAPHCONNCOMP, GRAPHISDAG,
% GRAPHISSPANTREE, GRAPHMAXFLOW, GRAPHMINSPANTREE, GRAPHPRED2PATH,
% GRAPHSHORTESTPATH, GRAPHTOPOORDER, GRAPHTRAVERSE.  
%
% References: 
%  [1]	S. Fortin "The Graph Isomorphism Problem" TR 96-20, Dept. of
%       Computer Science, University of Alberta, 1996. 
%  [2]	B.D. McKay "Practical Graph Isomorphism" Congressus Numerantium,
%       1981. 

%   Copyright 2006 The MathWorks, Inc.


debug_level = 0;
directed = true;

% read in optional PV input arguments
nvarargin = numel(varargin);
if nvarargin
    if rem(nvarargin,2) == 1
        error(message('bioinfo:graphisomorphism:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'directed'};
    for j=1:2:nvarargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:graphisomorphism:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:graphisomorphism:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % 'directed'
                    directed = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            end
        end
    end
end

% call the mex implementation of the graph algorithms
if nargout>1
   [F,M] = graphalgs('isomorphism',debug_level,directed,G1,G2);
else
   F = graphalgs('isomorphism',debug_level,directed,G1,G2);
end    

