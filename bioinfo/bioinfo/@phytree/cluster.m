function [clus,nclus,steps] = cluster(tr,v,varargin)
%CLUSTER Construct clusters from a phylogenetic tree.
%
% I = CLUSTER(T,V) returns a column vector containing a cluster index for
% each species (leaf) in a phylogenetic tree object. It determines the
% optimal number of clusters as follows:
%    - Starting with two clusters (k = 2), selects the tree cut that 
%      optimizes an specified criterion
%    - Increments k by 1 and again selects a new cut that added to
%      the previous partition optimizes the criterion
%    - Continues incrementing k and selecting the optimal cuts for each k
%      until a criterion value satisfies the scalar threshold value V or k
%      reaches the maximum number of clusters (which by default is the
%      number of leaves)
%    - From all the optimal cuts at every k value, selects the partition
%      that optimizes the criterion
%
% CLUSTER(...,'CRITERION',C) sets the criterion used for determining the
% number of clusters as a function of the species pairwise distances. The
% available options are: 
%
%   'maximum'    - Maximum within cluster pairwise distance (Wmax). Smaller
%   (default)      is better, cluster splitting stops when Wmax<=V.  
%      
%   'median'     - Median within cluster pairwise distance (Wmed). Smaller
%                  is better, cluster splitting stops when Wmed<=V. 
%
%   'average'    - Average within cluster pairwise distance (Wavg). Smaller
%                  is better, cluster splitting stops when Wavg<=V. 
%
%   'ratio'      - Between/within cluster pairwise distance ratio defined
%                  as:   
%                        BWrat = (trace(B)/(k-1)) / (trace(W)/(n-k))
%
%                  where B and W are the between/within-scatter matrices, k
%                  is number of clusters, and n is number of species in the
%                  tree. Larger is better, cluster splitting stops when
%                  BWrat>=V.  
%
%   'gain'       - Within cluster pairwise distance gain defined as:
%
%                        Wgain = (trace(Wold)/trace(W)-1) * (n-k-1)
%
%                  where W and Wold are the within-scatter matrices for k
%                  and k-1, k is number of clusters, and n is number of
%                  species in the tree. Larger is better, cluster splitting
%                  stops when Wgain<=V.
%
%   'silhouette' - Average silhouette width (SWavg). SWavg ranges from -1
%                  to +1, larger is better, cluster splitting stops when
%                  SWavg>=V. 
%
% CLUSTER(...,'MAXCLUST',N) sets the maximum number of possible clusters
% for the tested partitions. N defaults to number of leaves in the tree.
% 
% Notes:
% 1. When using the 'ratio', 'gain', or 'silhouette' criterions, it is hard
% to know a good threshold value V in advance. Set V to [] (empty) to find
% the optimal number of clusters below N. Set also N to a small value to
% avoid expensive computation by testing all possible number of clusters.  
% 2. When using the 'maximum', 'median', or 'average' criterions, the
% optimal number of clusters is equal to N because such metric monotonically
% decrease as k increases. Set V to [] (empty) to force CLUSTER to return N
% clusters.
% 
% CLUSTER(...,'DISTANCES',D) substitutes the patristic (tree) distances
% with a user provided pairwise distance matrix D. For example, D can be 
% the real sample pairwise distances.
%
% [I,J] = CLUSTER(...) returns cluster indices for all the nodes in the
% tree, leaf nodes and branch nodes. 
%
% [I,J,S] = CLUSTER(...) returns the branch being considered for the cut
% and the value of criterion at each step of the algorithm. Set V to []
% (empty) and N to its default to obtain the whole curve of the criterion
% versus the number of clusters. Observe that the computation of some
% criterions may be computationally intensive.
%
% Example:
%    % Find the best partition and number of clusters for a neighbor
%    % joining tree built from a multiple alignment:
%    gagaa = multialignread('aagag.aln');
%    gag_tree = seqneighjoin(seqpdist(gagaa),'equivar',gagaa);
%    [i,j] = cluster(gag_tree,[],'criterion','gain','maxclust',10);
%    h = plot(gag_tree);
%    set(h.BranchLines(j==2),'Color','b')
%    set(h.BranchLines(j==1),'Color','r')
%
%   See also CLUSTER, PHYTREE, PHYTREEREAD, PHYTREE/SUBTREE, SEQLINKAGE,
%   SEQNEIGHJOIN, SEQPDIST, SILHOUETTE.

% Copyright 2009-2010 The MathWorks, Inc.


% References:
%
% [1] Dudoit S, Fridlyan J, A prediction-based resampling method for
%     estimating the number of clusters in a dataset. Genome Biology,
%     3(7):research0036.10036.21, 2002.
%
% [2] Theodoridis, S. and Koutroumbas, K. Pattern Recognition, Academic
%     Press, pp. 434-435, 1999. 
%
% [3] Kaufman L, Rousseeuw PJ, Finding Groups in Data: An Introduction to
%     Cluster Analysis. New York, Wiley, 1990. 
%
% [4] Calinski R, Harabasz J, A dendrite method for cluster analysis.
%     Commun Statistics, 3:1-27, 1974.
%
% [5] Hartigan JA, Statistical theory in clustering, J Classification,
%     2:63-76, 1985.

numBranches = size(tr.tree,1);
numLeaves = numBranches + 1;
numLabels = numBranches + numLeaves;

[criteria,P,n] = parse_inputs(tr,numLeaves,varargin{:});

% Check v input and initialize properly in case v=[]
if ~isempty(v) && (~isnumeric(v) || ~isscalar(v))
    error(message('bioinfo:phytree:cluster:InvalidThreshold'))
end
if isempty(v)
    if criteria=='s'  || criteria=='r'
        v = inf;
    else
        v = -inf;
    end
end
if criteria=='s'  || criteria=='r'
    v = -v;
end

% Find all possible binary clusterizations based on the input tree
bc = false(numLabels,numLabels);
for i = 1:numLeaves
    bc(i,i) = true;
end
for i = 1:numBranches
    bc(:,i+numLeaves) = any(bc(:,tr.tree(i,:)),2);
    bc(i+numLeaves,i+numLeaves) = true;
end
selLeaves = 1:numLabels <= numLeaves;

% Initialize all species and edges to the first cluster
eclus = ones(numLabels,1);

% Initialize output 'steps'
steps = zeros(numLeaves,2);
switch criteria % for k = 1
    case 'm'
        steps(1,:) = [NaN max(squareform(P))]; 
    case 'd'
        steps(1,:) = [NaN median(squareform(P))]; 
    case 'a'
        steps(1,:) = [NaN mean(squareform(P))];
    case 's'
        steps(1,:) = [NaN 1];
    case 'r'
        steps(1,:) = [NaN 0];
        P2 = P.^2;
        TP2 = sum(P2(:));
    case 'g'
        P2 = P.^2;
        Wold = sum(squareform(P2))./numLeaves;
        steps(1,:) = [NaN inf];
end

% Initialize bcsel; used for marking edges that had been removed for
% creating partitions. bcsel(end) is the root and is invalidated from the
% beginning, no partition is valid at the root.
bcsel = false(numLabels,1);
bcsel(end) = true;

cache =  0;
k = 1;
% Main loop: repeat for every k>1 until exit conditions apply
while (k<numLeaves) && (k<n) && (steps(k,2)>v)
   
    k = k+1;
    pospart = find(~bcsel);
    mas = zeros(numel(pospart),1);
    nma = zeros(numel(pospart),2);
    oma = bsxfun(@times,bsxfun(@ne,1:k-1,eclus(pospart)),cache);
    switch criteria
        case 'm'
            for i = 1:numel(pospart)
                thisclu = eclus(selLeaves)==eclus(pospart(i));
                idx = thisclu & bc(selLeaves,pospart(i));
                if sum(idx)>1
                    nma(i,1) = max(squareform(P(idx,idx)));
                end
                idx = thisclu & ~bc(selLeaves,pospart(i));
                if sum(idx)>1
                    nma(i,2) = max(squareform(P(idx,idx)));
                end
                
            end
            mas = max([oma nma],[],2);
        case 'd'
            for i = 1:numel(pospart)
                thisclu = eclus(selLeaves)==eclus(pospart(i));
                idx = thisclu & bc(selLeaves,pospart(i));
                if sum(idx)>1
                    nma(i,1) = median(squareform(P(idx,idx)));
                end
                idx = thisclu & ~bc(selLeaves,pospart(i));
                if sum(idx)>1
                    nma(i,2) = median(squareform(P(idx,idx)));
                end
                
            end
            mas = max([oma nma],[],2);            
        case 'a'
            for i = 1:numel(pospart)
                thisclu = eclus(selLeaves)==eclus(pospart(i));
                idx = thisclu & bc(selLeaves,pospart(i));
                if sum(idx)>1
                    nma(i,1) = mean(squareform(P(idx,idx)));
                end
                idx = thisclu & ~bc(selLeaves,pospart(i));
                if sum(idx)>1
                    nma(i,2) = mean(squareform(P(idx,idx)));
                end
                
            end
            mas = max([oma nma],[],2);             
        case 's' % Reference [3] (we use the negative of the criterion since 
                   % this implementation searches for the cut that
                   % minimizes the criterion while the reference indicates
                   % that larger values of the silhouette  are better)
            for i = 1:numel(pospart)
                tclus = eclus(selLeaves) + k.*bc(selLeaves,pospart(i));
                sili = silhouette((1:numLeaves)',tclus,@(x,y) P(x,y));
                mas(i) = - mean(sili);
            end
        case 'r' % Reference [4] (we use the negative of the criterion since 
                   % this implementation searches for the cut that
                   % minimizes the criterion while the reference indicates
                   % that larger values of the ratio are better)
            for i = 1:numel(pospart)
                thisclu = eclus(selLeaves)==eclus(pospart(i));
                idx = thisclu & bc(selLeaves,pospart(i));
                if sum(idx)>1
                    nma(i,1) = sum(squareform(P2(idx,idx)))./sum(idx);
                end
                idx = thisclu & ~bc(selLeaves,pospart(i));
                if sum(idx)>1
                    nma(i,2) = sum(squareform(P2(idx,idx)))./sum(idx);
                end
            end
            mas = (1-TP2./sum([oma nma],2)) .* (numLeaves-k) ./ (k-1);
            mas(isnan(mas)) = -inf;
        case 'g' % Reference [5] 
            for i = 1:numel(pospart)
                thisclu = eclus(selLeaves)==eclus(pospart(i));
                idx = thisclu & bc(selLeaves,pospart(i));
                if sum(idx)>1
                    nma(i,1) = sum(squareform(P2(idx,idx)))./sum(idx);
                end
                idx = thisclu & ~bc(selLeaves,pospart(i));
                if sum(idx)>1
                    nma(i,2) = sum(squareform(P2(idx,idx)))./sum(idx);
                end
            end
            mas = (numLeaves-k-1).*(Wold./sum([oma nma],2)-1);
    end
    
    % Find the cut that leads to the minimization of the criterion
    if criteria == 'g'
        bcri = max(mas);
    else
        bcri = min(mas);
    end
    mins = pospart(bcri == mas);
    if numel(mins)>1
        % From all the possible cuts that give the minimum criterion select
        % the longest edge for the cut:
        [~,h] = max(tr.dist(mins));
        bidx = mins(h);
    else
        bidx = mins;
    end
    % Find cluster being cut:
    kk = eclus(bidx);
    % Create new cluster:  
    eclus((eclus==kk) & bc(:,bidx)) = k;
    % Update pre-calculated criteria on current clusters
    cache([k kk]) = nma(pospart==bidx,:);
        
    % In case of criteria is 'gain' recalculate the best within-scatter
    % matrix:
    if criteria == 'g'
        Wold = sum([oma(pospart==bidx,:) nma(pospart==bidx,:)]);
    end
    
    % Mark edge being cut
    bcsel(bidx) = true;
    % Remove edges in branches where the other two adjacent edges had
    % already been cut:
    tbcsel = ~bcsel;
    while ~isequal(tbcsel,bcsel)
        tbcsel = bcsel;
        q = sum([bcsel(tr.tree) bcsel((numLeaves+1):end)],2)==2;
        bcsel(tr.tree(q,:)) = true;
        bcsel(numLeaves+find(q)) = true;
    end
    % Update output 'steps'
    steps(k,:) = [bidx,bcri];
end

steps = steps(1:k,:);
if criteria == 'g'
    [~,optk] = max(steps(2:end,2));
    optk = optk +1;
else
    [~,optk] = min(steps(:,2));
end
optk = max(2,optk);

% calculate cluster for each node with propagation towards leaves
nclus = ones(numLabels,1);
cuts = sort(steps(2:optk,1),'descend');
for i = 1:numel(cuts)
    k = i+1;
    nclus(cuts(i))= k;
    for j = cuts(i)-numLeaves : -1 :1
        if nclus(j+numLeaves)==k
            nclus(tr.tree(j,:))= k;
        end
    end
end

% first output argument is only cluster index for each leaf
clus = nclus(1:numLeaves);

if criteria == 's' || criteria == 'r'
   steps(:,2) = -steps(:,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [criteria,P,n] = parse_inputs(tr,numLeaves,varargin)
% Parse the varargin parameter/value inputs

% Check that we have the right number of PVP inputs
if rem(numel(varargin),2)~= 0
    error(message('bioinfo:phytree:cluster:IncorrectNumberOfArguments', mfilename));
end

% The allowed inputs
okargs = {'criterion','distances','maxclust'};

% Set default values
P = squareform(pdist(tr)); % Pairwise distances, stored in square form 
                           % since it will be indexed by row and columns
n = inf; % default maximum number of clusters
criteria = 'm'; % maximum is default criteria

for j=1:2:numel(varargin)-1
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, ['phytree:' mfilename]);
    switch(k)
        case 1 % 'criterion'
            crits = {'maximum','average','median','silhouette','gain','ratio'};
            crit = strmatch(lower(pval),crits);
            if isempty(crit)
                error(message('bioinfo:phytree:cluster:NotValidCriterion'))
            elseif numel(crit)>1
                error(message('bioinfo:phytree:cluster:AmbiguousCriterion'))
            else
                codes = 'madsgr';
                criteria = codes(crit);
            end
        case 2 % 'distances'
            if isnumeric(pval) && isequal(size(pval),[numLeaves,numLeaves]) && all(~diag(pval)) && isequal(pval,pval') 
                P = pval;
            elseif isnumeric(pval) && isvector(pval) && numel(pval)==(numLeaves*(numLeaves-1)/2)
                P = squareform(pval);
            else
                error(message('bioinfo:phytree:cluster:InvalidDistances')) 
            end
        case 3 % 'maxclust'
            if ~isnumeric(pval) || ~isscalar(pval) || rem(pval,1) || pval<1
                 error(message('bioinfo:phytree:cluster:InvalidMaxClust'))
            end
            n = pval;
    end
end




 
