function [tr,order] = reorder(tr,order,varargin)
%REORDER changes the leaf order of a phylogenetic tree.
%
%   T2 = REORDER(T1,ORDER) changes the order of the leaves of the
%   phylogenetic tree T1 without modifying its structure and distances.
%   ORDER is a vector with position indices for each leaf.
%
%   [T2 OPT_ORDER] = REORDER(T1,ORDER,'APPROXIMATE',TRUE) uses optimal
%   leaf-ordering calculation to find the closest order possible to the
%   suggested one without dividing the clades (or having crossing branches).
%   Default is FALSE, and REORDER returns an error if the suggested order
%   is not valid.
%
%   [T3 OPT_ORDER] = REORDER(T1,T2) uses optimal leaf-ordering calculation
%   to find the closest order possible of the leaves in T1 to match the
%   order of the leaves in T2 without dividing the clades (or having
%   crossing branches). 
%
%   Examples:
%
%   % Permute the leaves of a tree with valid order:
%   b = [1 2; 3 4; 5 6; 7 8;9 10];
%   TR = phytree(b);   
%   view(TR)
%   TR_REORDERED = reorder(TR,[5 6 3 4 1 2])
%   view(TR_REORDERED)
%
%   % Reorder a tree the best possible to match an alphabetical order: 
%   TREE = phytreeread('pf00002.tree');
%   [dummy,order] = sort(get(TREE,'LeafNames'));
%   TREE_REORDERED = reorder(TREE,order,'approximate',true)
%   view(TREE)
%   view(TREE_REORDERED)
%
%   % Build a phylogenetic tree with two different methods and then
%   % reorder one of them for an easier visual comparison:
%   seqs = fastaread('pf00002.fa');
%   dist = seqpdist(seqs,'method','jukes-cantor','indels','pair');
%   NJ = seqneighjoin(dist,'equivar',seqs)  % neighbor joining 
%   HC = seqlinkage(dist,'single',seqs)     % single hierarchical clustering
%   % reorder HC to match the best possible the leaves of NJ:
%   HC_reordered = reorder(HC,NJ)
%   view(NJ)
%   view(HC_reordered)
%
%   See also PHYTREE, PHYTREE/GET, PHYTREE/GETBYNAME, PHYTREE/PRUNE,
%   PHYTREE/REROOT, PHYTREE/SELECT, SEQNEIGHJOIN. 

% Copyright 2006-2012 The MathWorks, Inc.


if numel(tr)~=1
     error(message('bioinfo:phytree:reorder:NoMultielementArrays'));
end
numBranches = size(tr.tree,1);
numLeaves = numBranches + 1;
numLabels = numBranches + numLeaves; 

% validate second input (is it a tree or a vector with the siggested
% order?)
guidingTree = false;
if isa(order,'phytree')
    guidingTree = true;
    trg = order;
elseif numel(order)~=numLeaves || any(sort(order(:))~=(1:numLeaves)')
    error(message('bioinfo:phytree:reorder:InvalidOrder'));
end

doOptimalLeafOrdering = false;
% read in optional PV input arguments
nvarargin = numel(varargin);
if nvarargin
    if rem(nvarargin,2) == 1
        error(message('bioinfo:phytree:reorder:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'approximate',''};
    for j=1:2:nvarargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:phytree:reorder:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:phytree:reorder:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % 'approximate'
                    if guidingTree
                        warning(message('bioinfo:phytree:reorder:PVignored'))
                    else
                        doOptimalLeafOrdering = bioinfoprivate.opttf(pval);
                        if isempty(doOptimalLeafOrdering)
                            error(message('bioinfo:phytree:reorder:doOptimalLeafOrderingOptionNotLogical', upper( char( okargs( k ) ) )));
                        end
                    end
            end
        end
    end
end

if guidingTree

    % match leaf names of the tree to the guiding tree
    [lastwmsg,lastwid] = lastwarn;
    numLeavesG = (numel(trg.dist)+1)/2;
    warnState = warning('off','bioinfo:seqmatch:StringNotFound');
    ptrs = seqmatch(tr.names(1:numLeaves),trg.names(1:numLeavesG),'exact',true);
    warning(warnState);
    lastwarn(lastwmsg,lastwid);

    % propagate most desired location of groups to the branches
    dl = [ptrs;zeros(numBranches,1)];
    sz = [ptrs>0;zeros(numBranches,1)];
    for i = 1:numBranches
        sz(numLeaves+i) = sum(sz(tr.tree(i,:)).*(dl(tr.tree(i,:))~=0));
        if sz(numLeaves+i)>0
            dl(numLeaves+i) = sum(dl(tr.tree(i,:)).*sz(tr.tree(i,:)))/sz(numLeaves+i);
        end
    end
    % propagate values to leaves with missing values
    for i = numBranches:-1:1
        if dl(tr.tree(i,1)) == 0
            dl(tr.tree(i,1)) = dl(numLeaves+i);
        end
        if dl(tr.tree(i,2)) == 0
            dl(tr.tree(i,2)) = dl(numLeaves+i);
        end
    end
    % normalize to account for different sized trees
    ptrs = dl(1:numLeaves)*numLeavesG/numLeaves;
    % call optimal leaf order with the suggested order
    order = optimalleaforder(tr.tree,pdist(ptrs),'criteria','group');
    % reverse order to better match suggested order if necessary
    if sum(abs(diff([ptrs,order'],[],2)))>sum(abs(diff([ptrs,fliplr(order)'],[],2)))
        order = fliplr(order);
    end

    % reorder the tree
    permuta = [order(:);(numLeaves+1:numLabels)'];
    ipermuta(permuta) = 1:numLabels;
    tr.tree = ipermuta(tr.tree);
    tr.dist = tr.dist(permuta);
    tr.names = tr.names(permuta);


else % guidingTree==false, user suggests an arbitrarily order

    % reorder the tree
    permuta = [order(:);(numLeaves+1:numLabels)'];
    ipermuta(permuta) = 1:numLabels;
    tr.tree = ipermuta(tr.tree);
    tr.dist = tr.dist(permuta);
    tr.names = tr.names(permuta);

    % check if the reordered tree structure leads to a tree with non-divided
    % clusters (no crossing branches when displayed):
    isInvalidOrder = false;
    mi = [1:numLeaves,zeros(1,numBranches)];
    ma = [1:numLeaves,zeros(1,numBranches)];
    sz = ones(numBranches+numLeaves,1);
    for i = 1:numBranches
        j = i+numLeaves;
        mi(j) = min(mi(tr.tree(i,:)));
        ma(j) = max(ma(tr.tree(i,:)));
        sz(j) = sum(sz(tr.tree(i,:)));
        if ma(j)-mi(j)+1 ~= sz(j)
            isInvalidOrder = true;
            break;
        end
    end

    if isInvalidOrder
        % suggested order is invalid, try to fix it or error out:
        if doOptimalLeafOrdering
            order = optimalleaforder(tr.tree,pdist((1:numLeaves)'),'criteria','group');
            permuta = [order(:);(numLeaves+1:numLabels)'];
            ipermuta(permuta) = 1:numLabels;
            tr.tree = ipermuta(tr.tree);
            tr.dist = tr.dist(permuta);
            tr.names = tr.names(permuta);
        else
            error(message('bioinfo:phytree:reorder:UnorderedClusters'))
        end
    end

end



    
