function [tr order] = phytree(varargin)
%PHYTREE Phylogenetic tree object.
%
%  TREE = PHYTREE(B) creates an ultrametric phylogenetic tree object. B is
%  a numeric array of size [NUMBRANCHES X 2] in which every row represents
%  a branch of the tree and it contains two pointers to the branches 
%  or leaves nodes which are its children. Leaf nodes are numbered from 1
%  to NUMLEAVES and branch nodes are numbered from NUMLEAVES + 1 to
%  NUMLEAVES + NUMBRANCHES. Note that since only binary trees are allowed,
%  then NUMLEAVES = NUMBRANCHES + 1. Branches are defined in chronological
%  order, i.e. B(i,:) > NUMLEAVES + i. As a consequence, the  first row can
%  only have pointers to leaves and the last row must represent the 'root'
%  branch. Parent-child distances are set to the unit or by the ultrametric
%  condition if child is a leaf.
%
%  TREE = PHYTREE(B,D) creates an additive phylogenetic tree object with
%  branch distances defined by D. D is a numeric array of size [NUMNODES X
%  1] with the distances of every child node (leaf or branch) to its parent
%  branch. NUMNODES = NUMLEAVES + NUMBRANCHES. D(end), the distance
%  associated to the root node, is meaningless.
%
%  TREE = PHYTREE(B,C) creates an ultrametric phylogenetic tree object with
%  branch distances defined by C. C is a numeric array of size [NUMBRANCHES
%  X 1] with the coordinates of every branch node. In ultrametric tress all
%  the leaves are at the same location (i.e. same distance to the root).
%
%  TREE = PHYTREE(BC) creates an ultrametric phylogenetic binary tree
%  object with branch pointers in BC(:,[1 2]) and branch coordinates in
%  BC(:,3). Same as PHYTREE(B,C).
%
%  TREE = PHYTREE(...,N) specifies the names for the leaves and/or the
%  branches. N is a cell of strings. If NUMEL(N)==NUMLEAVES then the names
%  are assigned chronologically to the leaves. If NUMEL(N)==NUMBRANCHES the
%  names are assigned to the branch nodes. If NUMEL(N)==NUMLEAVES +
%  NUMBRANCHES all the nodes are named. Unassigned names default to 'Leaf
%  #' and/or 'Branch #' as required.
% 
%  TREE = PHYTREE creates an empty phylogenetic tree object.
%
%  Example: 
%
%      % create an ultrametric tree
%      b = [1 2; 3 4; 5 6; 7 8;9 10];
%      t = phytree(b);   
%      view(t)
% 
%      % create an ultrametric tree with specified branch distances
%      bd = [.1 .2 .3 .3 .4 ]';
%      b = [1 2; 3 4; 5 6; 7 8;9 10];
%      t = phytree(b,bd);  
%      view(t)
%
%  See also PHYTREE/GET, PHYTREE/SELECT, PHYTREEREAD, PHYTREETOOL,
%   PHYTREEWRITE, SEQLINKAGE, SEQNEIGHJOIN, SEQPDIST.


% Copyright 2003-2012 The MathWorks, Inc.


justVerifyValidity = false;

switch nargin
    case 0
       tr.tree = zeros(0,2);
       tr.dist = zeros(0,1);
       tr.names = {};
    case 1
        B = varargin{1};
    case 2
        B = varargin{1};
        if iscell(varargin{2})
            N = varargin{2};
        else
            D = varargin{2};
        end
    case 3
        B = varargin{1};
        D = varargin{2};
        N = varargin{3};
    otherwise
        error(message('bioinfo:phytree:phytree:IncorrectNumberOfArguments', mfilename));
end

if nargin==1 && isstruct(B) && isfield(B,'tree') && isfield(B,'names') && isfield(B,'dist')
    N = B.names;
    D = B.dist;
    B = B.tree;
    tr.tree = B;
    tr.dist = D;
    tr.names = N;
    justVerifyValidity = true;
end

if nargin
    if isnumeric(B)
        switch size(B,2)
            case 2
                % ok
            case 3
                D = B(:,3);
                B(:,3)=[];
            otherwise
                error(message('bioinfo:phytree:phytree:IncorrectSize'))
        end
    else
        error(message('bioinfo:phytree:phytree:IncorrectType'))
    end

    % test B
    if sum(diff(sort(B(:)))~=1) || (min(B(:))~=1)
        error(message('bioinfo:phytree:phytree:IncompleteTree'))
    end
    numBranches = size(B,1);
    numLeaves = numBranches + 1;
    numLabels = numBranches + numLeaves;
    h=all(B'>=repmat(numLeaves+1:numLabels,2,1));
    if any(h)
        error(message('bioinfo:phytree:phytree:NonChronologicalTree', num2str( find( h ) )))
    end

    if exist('D','var')
        if ~isnumeric(D) || any(D(:)<0) || ~all(isreal(D)) || size(D,2)~=1
            error(message('bioinfo:phytree:phytree:DistancesNotValid'))
        end
        switch size(D,1)
            case numBranches
                D = [zeros(numLeaves,1); D]; % add ultrametric distances of leaves
                D(B) = D((numLeaves+(1:numBranches))'*[1 1])-D(B);  %dist of edges
                D(end) = 0; % set root at zero
            case numLabels
                % ok
            otherwise
                error(message('bioinfo:phytree:phytree:DistancesBadSize'))
        end
    else % set defaut D
        % look for parents
        P = zeros(numLabels,1);
        P(B) = repmat((1:numBranches)',1,2);
        P(end) = numBranches;
        % look at which level is every branch
        L = zeros(numLabels,1);
        for ind = 1:numBranches
            L(ind+numLeaves) = max(L(B(ind,:))+1);
        end
        D = L(P+numLeaves)-L;
    end

    % set default names
    for ind = 1:numLeaves
        names{ind}=['Leaf ' num2str(ind)]; %#ok
    end
    for ind = 1:numBranches
        names{ind+numLeaves}=['Branch ' num2str(ind)];
    end

    if exist('N','var')
        if ~iscell(N)
            error(message('bioinfo:phytree:phytree:NamesNotValid'))
        end
        switch numel(N)
            case numLabels
                h = 1:numLabels;
            case numLeaves
                h = 1:numLeaves;
            case numBranches
                h = numLeaves+1:numLabels;
            otherwise
                error(message('bioinfo:phytree:phytree:NamesBadSize'))
        end
        for ind = 1:length(h);
            str = N{ind};
            if ~ischar(str)
                error(message('bioinfo:phytree:phytree:NamesNotStrings'))
            end
            names{h(ind)}=str;
        end
        % check that none of the names is empty
        for ind = 1:numLabels
            if isempty(names{ind})
                if ind > numLeaves
                    names{ind} = ['Branch ' num2str(ind-numLeaves)];
                else
                    names{ind} = ['Leaf ' num2str(ind)];
                end
            end
        end
        if numel(unique(names))~=numLabels
            error(message('bioinfo:phytree:phytree:NamesNotUnique'))
        end
    end

    % check and corrects a non-monotonic tree
    monotonicWarning = false;
    for ind = 1:numBranches
        if any(D(B(ind,:))<0)
            monotonicWarning = true;
            tmp = min(D(B(ind,:)));
            D(B(ind,:)) = D(B(ind,:)) - tmp;
            D(numLeaves+ind) = D(numLeaves+ind) + tmp;
        end
    end
    if monotonicWarning
        D(end) = 0; %in the monotonic correction the root might have been shifted to negative numbers G335257
        warning(message('bioinfo:phytree:phytree:NonMonotonicTree'))
    end

    if justVerifyValidity
        tr = class(tr,'phytree');
        return
    end

    tr.tree = B;
    tr.dist = D;
    tr.names = names(:);
    
    % check if the current tree structure leads to a tree with non-divided
    % clusters (no crossing branches when displayed):
    needsReorder = false;
    mi = [1:numLeaves,zeros(1,numBranches)];
    ma = [1:numLeaves,zeros(1,numBranches)];
    sz = ones(numBranches+numLeaves,1);
    for i = 1:numBranches
        j = i+numLeaves;
        mi(j) = min(mi(tr.tree(i,:)));
        ma(j) = max(ma(tr.tree(i,:)));
        sz(j) = sum(sz(tr.tree(i,:)));
        if ma(j)-mi(j)+1 ~= sz(j)
            needsReorder = true;
            break;
        end
    end

    % reorder such that there will be no crossings in the displayed tree
    if nargout > 1
        if needsReorder
            [tr,order] = prettyorder(tr);
        else
            order = 1:numLeaves;
        end
    elseif needsReorder
        tr = prettyorder(tr);
    end

end %if nargin

% for trees of only one branch correct dimensions
% if size(tr.tree,2) <2 tr.tree = tr.tree'; end

% Makes the tree a class
tr = class(tr,'phytree');


