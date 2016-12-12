%% Building a Phylogenetic Tree for the Hominidae Species
% This example shows how to construct a phylogenetic trees from mtDNA
% sequences for the Hominidae taxa (also known as pongidae). This family
% embraces the gorillas, chimpanzees, orangutans and the humans.

%   Copyright 2003-2012 The MathWorks, Inc.


%% Introduction
% The mitochondrial D-loop is one of the fastest mutating sequence regions
% in animal DNA. Therefore, useful for comparing closely related organisms.
% The origin of modern man is a highly debated issue that has recently been 
% tackled by using mtDNA sequences. The limited genetic variability of
% human mtDNA has been explained in terms of a recent common genetic
% ancestry, thus implying that all modern-population mtDNAs originated from
% a single woman who lived in Africa less than 200,000 years.

%% Gather Sequence Data from GenBank(R)
% Get some data from GenBank(R), these are the accession codes for
% mitochondrial D-loop sequences that have been isolated for different
% hominidae species.  

%        Species Description      GenBank Accession
data = {'German_Neanderthal'      'AF011222';
        'Russian_Neanderthal'     'AF254446';
        'European_Human'          'X90314'  ;
        'Mountain_Gorilla_Rwanda' 'AF089820';
        'Chimp_Troglodytes'       'AF176766';
        'Puti_Orangutan'          'AF451972';
        'Jari_Orangutan'          'AF451964';
        'Western_Lowland_Gorilla' 'AY079510';
        'Eastern_Lowland_Gorilla' 'AF050738';
        'Chimp_Schweinfurthii'    'AF176722';
        'Chimp_Vellerosus'        'AF315498';
        'Chimp_Verus'             'AF176731';
       };
   
%%
% You can use the |getgenbank| function inside a for-loop to retrieve the
% sequences from the NCBI data repository and load them into MATLAB(R). 
%
%   for ind = 1:length(data)       
%       primates(ind).Header   = data{ind,1};
%       primates(ind).Sequence = getgenbank(data{ind,2},'sequenceonly','true');
%   end

%%
% For your convenience, previously downloaded sequences are included in a
% MAT-file. Note that data in public repositories is frequently curated and 
% updated; therefore the results of this example might be slightly different
% when you use up-to-date datasets.

load('primates.mat')

%% Build a UPGMA Phylogenetic Tree using Distance Methods 
% Compute pairwise distances using the 'Jukes-Cantor' formula and the
% phylogenetic tree with the 'UPGMA' distance method. Since the sequences
% are not pre-aligned, |seqpdist| will pairwise align them before computing
% the distances. 

distances = seqpdist(primates,'Method','Jukes-Cantor','Alpha','DNA');
UPGMAtree = seqlinkage(distances,'UPGMA',primates)

%%
% Render the UPGMA phylogenetic tree 
h = plot(UPGMAtree,'orient','top');
title('UPGMA Distance Tree of Primates using Jukes-Cantor model');
ylabel('Evolutionary distance')

%% Build a Neighbor-Joining Phylogenetic Tree using Distance Methods 
% Alternate tree topologies are important to consider when analyzing
% homologous sequences between species. A neighbor-joining tree can be 
% built using the |seqneighjoin| function. Neighbor-joining trees use the 
% pairwise distance calculated above (using the |seqpdist| function) to
% construct the tree. This method clusters using the minimum evolution
% method

NJtree = seqneighjoin(distances,'equivar',primates)

%%
% Render the NJ phylogenetic tree 
h = plot(NJtree,'orient','top');
title('Neighbor-Joining Distance Tree of Primates using Jukes-Cantor model');
ylabel('Evolutionary distance')

%% Comparing Tree Topologies
% Notice that the trees that are created have different topologies. The
% neighbor-joining tree groups Chimp Vellerosus in a clade with the 
% gorillas and the UPGMA tree groups it near chimps and orangutans. The
% |getcanonical| function can be used to compare these isomorphic trees.

sametree = isequal(getcanonical(UPGMAtree), getcanonical(NJtree))

%% Exploring the UPGMA Phylogenetic Tree
% Find the closest species to the 'European Human' entry (3).
names = get(UPGMAtree,'LeafNames')
[h_all,h_leaves] = select(UPGMAtree,'reference',3,'criteria','distance','threshold',0.6);

%%
% The logical index |h_all| selects the nodes (leaves and branches) within
% 0.6 of patristic distance to the 'European Human' leaf. The logical
% index |h_leaves| selects the leaf nodes within 0.6 of patristic distance
% to to the 'European Human' leaf.
subtree_names = names(h_leaves)

%%
% Reduce the tree to the sub-branch of interest
leaves_to_prune = ~h_leaves;
pruned_tree = prune(UPGMAtree,leaves_to_prune)
h = plot(pruned_tree,'orient','top');
title('Pruned UPGMA Distance Tree of Primates using Jukes-Cantor model');
ylabel('Evolutionary distance')

%%
% With |view| you can further explore/edit the phylogenetic tree using an
% interactive tool. See also |phytreeviewer|.
view(UPGMAtree,h_leaves)

%% References:
% [1] I.V. Ovchinnikov, et al., "Molecular analysis of Neanderthal DNA from
%     the northern Caucasus", Nature 404(6777), pp. 490-493, 2000.
%%
% [2] A. Sajantila, et al., "Genes and languages in Europe: an analysis of
%     mitochondrial lineages", Genome Research 5(1), pp. 42-52, 1995.
%%
% [3] M. Krings, et al., "Neandertal DNA sequences and the origin of modern
%     humans", Cell 90(1), pp. 19-30, 1997.
%%
% [4] M.I. Jensen-Seaman and K.K. Kidd, "Mitochondrial DNA variation and
%     biogeography of eastern gorillas", Molecular Ecology 10(9), pp.
%     2241-2247, 2001.

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20PRIMATESDEMO%20in%20Bioinformatics%20Toolbox%204.2 Suggest an enhancement for this example.>*

displayEndOfDemoMessage(mfilename)
