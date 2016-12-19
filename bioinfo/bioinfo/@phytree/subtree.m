function subtr = subtree(tr,nodes)
%SUBTREE Extracts a subtree.
%
%   T2 = SUBTREE(T1,NODES) Extracts a new subtree T2 in which the new root
%   is the first common ancestor of the NODES vector from T1. NODES in the
%   tree are indexed as [1:NUMLEAVES] for the leaves and as [NUMLEAVES+1 :
%   NUMLEAVES+NUMBRANCHES] for the branches. NODES can also be a logical
%   array of following sizes: [NUMLEAVES+NUMBRANCHES x 1], [NUMLEAVES x 1]
%   or [NUMBRANCHES x 1]. 
%
%   Example:
%
%      % Load a phylogenetic tree created from a protein family:
%      tr = phytreeread('pf00002.tree')
%       
%      % Get the subtree that contains the VIPR2 and GLR human proteins:
%      sel = getbyname(tr,{'vipr2_human','glr_human'});
%      sel = any(sel,2);
%      tr =  subtree(tr,sel)
%      view(tr);
%       
%   See also PHYTREE, PHYTREE/PRUNE, PHYTREE/SELECT, PHYTREE/GET,
%   PHYTREE/GETBYNAME. 

% Copyright 2003-2005 The MathWorks, Inc.


if numel(tr)~=1
     error(message('bioinfo:phytree:subtree:NoMultielementArrays'));
end
numBranches = size(tr.tree,1);
numLeaves = numBranches + 1;
numLabels = numBranches + numLeaves; 

% validate nodes
if islogical(nodes)
    if numel(nodes)==numLabels 
        nodes = nodes(:)==true;
    elseif numel(nodes)==numLeaves
        nodes = [nodes(:);false(numBranches,1)];
    elseif numel(nodes)==numBranches
        nodes = [false(numLeaves,1);nodes(:)];
    else
        error(message('bioinfo:phytree:subtree:IncorrectSizeInputVector'));
    end
elseif isnumeric(nodes) && isreal(nodes) && all(nodes(:)>=1) && all(nodes(:)<=numLabels)
    tem = false(numLabels,1);
    tem(floor(nodes(:))) = true;
    nodes=tem(:);
else
    error(message('bioinfo:phytree:subtree:InvalidInputNodes'));
end

% at this point NODES should only be a logical vector

if (~any(nodes(numLeaves+1:numLabels)) && (sum(nodes(1:numLeaves))<2))
    error(message('bioinfo:phytree:subtree:InvalidSubtree'));
end


% look for the first common ancestor that contains all selected nodes,
% accumulating the selected nodes towards the root, the common ancestor
% will be the first sum equal to the number of selected nodes
branchWidth = double(nodes);
for ind = 1:numBranches
    branchWidth(numLeaves+ind) = branchWidth(numLeaves+ind) + ...
                                 sum(branchWidth(tr.tree(ind,:)));
end
commonAncestor = find(branchWidth==sum(nodes),1);

% now propagate the ancestor new) root towards the leaves to find all the
% nodes that should stay for the subtree
sel = false(1,numLabels);
sel(commonAncestor) = true;
for ind = commonAncestor:-1:numLeaves+1
    sel(tr.tree(ind-numLeaves,:)) = sel(ind);
end
    
% extract the subtree
permuta = find(sel);
subtr=phytree;
ipermuta(permuta) = 1:length(permuta);
subtrNumLeaves = (ipermuta(end) + 1)/2;
subtr.tree = ipermuta(tr.tree(permuta(subtrNumLeaves+1:end)-numLeaves,:));
subtr.dist = tr.dist(permuta);
subtr.names = tr.names(permuta);
subtr.dist(end) = 0;
