%% Analyzing the Human Distal Gut Microbiome
% This example shows several ways of visualizing the results of functional
% metagenomic analyses. The discussion is based on two studies focusing on
% the metagenomic analysis of the human distal gut microbiome.

% Copyright 2008-2012 The MathWorks, Inc. 


%% Introduction
% The human distal gut is the highest density, natural bacterial ecosystem
% known to date. Its size - up to 100 trillions cells - far exceeds the
% size of all the human body's other microbial communities. Recent studies
% have shown that the gut microbiota helps regulate energy balance, both by
% extracting calories from otherwise indigestible components, and by
% controlling the storage of energy in adipocytes. Furthermore, the gut
% microbiota is involved in a myriad of bioprocesses ranging from the
% synthesis of essential vitamins to the metabolism of carbohydrates,
% lipids and other xenobiotics that we ingest.
% 
% For this example, we will use two data sets. The first data set consists
% of data resulting from the analysis of the distal gut microbiome of two
% adult American subjects [1]. It comprises a phylogenetic survey of the
% microbial communities and a functional analysis of the metabolic
% functions represented by the identified gene pool. The variables included
% in |dataset1| are described below. Note that the taxonomic assignments
% are represented as a nominal categorical array.

load gutmicrobiomedata.mat

%=== first data set variables
rank1 = dataset1.rank1; % superkingdom assignments of each hit
rank2 = dataset1.rank2; % phylum assignments of each hit
rank3 = dataset1.rank3; % class assignments of each hit
subjF = dataset1.subjF; % number of hits in female subject
subjM = dataset1.subjM; % number of hits in male subject

%% Taxonomic Profiling of Adult Human Distal Gut Microbiome 
% We perform a taxonomic profiling of the first data set by considering the
% taxonomic assignment of the contigs according to the best BLASTX hit.
%
% We start by computing the number of assigned hits that belong to each
% superkingdom for each subject.

l1 = getlabels(rank1); % superkingdom labels
n1 = numel(l1); 
count1 = zeros(n1,1); % number of hits for each superkingdom and subject

for i = 1:n1
    obs = rank1 == l1{i};
    count1(i,1) = sum(subjF(obs)); count1(i,2) = sum(subjM(obs));
end

%=== plot
figure(); barh(count1);
colormap(summer);
set(gca, 'yTicklabel', l1);
xlabel('Number of hits');
title('Superkingdom assignment of best-BLASTX-hits');
legend('SubjF', 'SubjM');

%% 
% As you can see from the bar plot, the microbial community living in the
% distal gut microbiome is prevalently of bacterial nature. The differences
% observed between the male and female subjects cannot be addressed due to
% the limited subject sample size and the possibility that these
% differences might be related to host genotype or lifestyle.
%
% We now repeat the analysis at the phylum level.

l2 = getlabels(rank2); % phylum labels
n2 = numel(l2); 
count2 = zeros(n2,1); % number of hits for each phylum and subject

for i = 1:n2
    obs = rank2 == l2{i};
    count2(i,1) = sum(subjF(obs)); count2(i,2) = sum(subjM(obs));
end

%=== plot
figure(); barh(count2);
colormap(summer);
set(gca, 'ytick', 1:n2, 'yTicklabel', l2);
xlabel('Number of hits');
title('Phylum assignment of best-BLASTX-hits');
legend('SubjF', 'SubjM');

%%
% The bacterial phylotypes are assigned mostly to two divisions, the
% Firmicutes and the Actinobacteria. The relative paucity of Bacteroidetes
% assignments conflicts with data from other studies, but the discrepancy
% might be caused by the known biases of fecal lysis and DNA extraction
% methods used.
%
% Finally, we perform the same steps on the assignments at the class level.

l3 = getlabels(rank3); % class labels
n3 = numel(l3);
count3 = zeros(n3,1); % number of hits for each class and subject

for i = 1:n3
    obs = rank3 == l3{i};
    count3(i,1) = sum(subjF(obs)); count3(i,2) = sum(subjM(obs));
end

%=== plot
figure(); barh(count3);
colormap(summer);
set(gca, 'ytick', 1:n3, 'yTicklabel', l3);
xlabel('Number of hits');
title('Class assignment of best-BLASTX-hits');
legend('SubjF', 'SubjM');

%%
% The taxonomic distribution at the class level reveals an abundance of
% bacterial phylotypes in the Clostridia and Bacilli groups, and also
% Actinobacteria and Methanobacteria. 

%% Combining Taxonomic Distribution and the Underlying Classification
% You can combine the taxonomic distribution and the underlying taxonomic
% classification into a single representation by using a graph where each
% leaf node represents a class, and each internal node represents a phylum
% or a superkingdom.
%
% To construct such a graph, we need to determine the connectivity matrix
% |CM| representing the parent-child relationships among the nodes. We
% identify the phyla (children) belonging to each superkingdom (parent),
% and in turn the classes (children) belonging to each phylum (parent).

L = nominal([l1 l2 l3]);
N = n1 + n2 + n3;
CM = zeros(N, N); % connectivity matrix

%=== populate CM with relationships between superkingdoms and phyla
for i = 1:n1
    obs = rank1 == l1{i}; % entries classified in a given superkingdom 
    from = find(L == l1{i});  % parent node
    subobs =  unique(rank2(obs)); % phyla in a given superkigdom
    for j = 1:numel(subobs)
        to = find(L == subobs(j)); % child node
        CM(from, to) = 1;
    end
end

%=== populate CM with relationships between phyla and classes
for i = 1:n2
    obs = rank2 == l2{i}; % entries classified in a given phylum 
    from = find(L == l2{i});  
    from = from(end);
    subobs =  unique(rank3(obs)); % classes in a given phylum
    for j = 1:numel(subobs)
        to = find(L == subobs(j)); 
        to = to(end);
        CM(from, to) = 1;
    end
end

%=== create biograph object
bg = biograph(CM-diag(diag(CM)),[],'NodeAutoSize','off','ShowTextInNodes','Label');

%%
% The resulting graph has 60 nodes and 58 edges. Each level in the graph is
% associated with a given taxonomic rank, and the edges represent the
% underlying taxonomic classification. We can now label each node with the
% corresponding taxonomic assignment and rotate the entire graph
% counterclockwise by 90 degrees.

%=== label each node
set(bg.Nodes,'Size',[10 100]);
for i = 1:numel(bg.Nodes)
   bg.Nodes(i).Label = char(L(i));
end
dolayout(bg);

%=== rotate counterclockwise by 90 degrees
for i = 1:numel(bg.Nodes)
    bg.Nodes(i).Position = fliplr(bg.Nodes(i).Position).*[-1 1];
    bg.Nodes(i).Size = [100 15];
end

%=== redraw edges without changing node positions
bg.LayoutType = 'equilibrium';
dolayout(bg,'PathsOnly',true);
view(bg)


%%
% To include the distribution data in the graph, for each assignment we
% consider the average number of hits between the two subjects and the
% corresponding percentage. We then customize the color and size of each
% node. In particular, leaf nodes are represented with boxes, while
% internal nodes are represented with circles, and the size of each node is
% proportional to the number of hits (percentage) that fall within a given
% taxonomic assignment.

%=== compute distribution among ranks
count = [count1; count2; count3];
count = sum(count,2)/2; % avg between subjects
pct = (count + 1)/sum(count + 1) * 100; % add pseudocounts

%=== determine color schema 
t = accumarray(round(pct+1),1);
t(t>0) = 1:nnz(t);
colors = flipud(summer(nnz(t)));
cindex = t(round(pct+1));

%=== customize color of nodes according to distribution
for i = 1:numel(bg.Nodes)
    mynode = bg.Nodes(i);
    if (numel(getdescendants(mynode))~= 1) % leaf
        mynode.Shape = 'circle';
    end
    mynode.Color = colors(cindex(i),:);
end

view(bg)

%%
% From this representation, you can immediately see how the majority of the
% microbial communities are composed of Bacteria, in particular Firmicutes,
% including Clostridia and Bacilli.

%% Comparative Functional Analysis Using KEGG Categories
% Phylogenetic assessments of microbial communities provide a starting
% point for interpreting the functional predictions from metagenomic data.
% The metabolic potential of the microbiota is studied to understand how
% the human distal gut microbiome provides us with physiological properties
% that we have not had to evolve on our own.
%
% Here we consider the metabolic functions associated with the human distal
% gut microbiome through KEGG pathways assignments. We use odds ratios to
% rank the enrichment or depletion of KEGG categories with respect to
% reference genomic data sets, namely the _Homo sapiens_ genome, a
% collection of sequenced bacterial genomes, and a collection of the
% sequenced archaeal genomes.

genome = dataset1.genome; % reference genomes considered
keggCat = dataset1.keggCat; % KEGG category assignment
keggData = dataset1.keggData; % odds ratio for each KEGG category relative to reference genomes

%%
% An odds ratio of one (corresponding to a log of zero) indicates that the
% microbial community had the same proportion of hits to a given category
% as the reference data set. An odds ratio greater than one (corresponding
% to a log greater than zero) indicates enrichment, whereas an odds ratio
% less than one (corresponding to a log less than zero) indicates
% under-representation with respect to the reference data set. 

hi = imagesc(log(keggData));
colormap(redbluecmap);
colorbar;
ha = get(hi, 'Parent');
set(ha, 'XTick', 1:3, 'XTickLabel', genome);
set(ha, 'YTickLabel', keggCat);

%% 
% From the heat map above, we notice that the human gut microbiome is
% highly enriched relative to the human genome, similar to the sequenced
% bacteria, and moderately enriched relative to the sequenced archaea. 

%% Comparative Functional Analysis Using COG Categories
% COG categories, which use evolutionary relationships to group
% functionally related genes, can be used to perform functional analysis
% instead of KEGG categories, which map enzymes onto known metabolic
% pathways. The DataMatrix object |dm2| consists of data resulting from a
% comparative metagenomic analysis of the human distal gut microbiome of
% several Japanese subjects, including infants, children and adults [2].
% For reference, the data of American subjects considered above as well as
% other metagenomic data sets are reported. The rows represent the various
% COG observations, whereas the columns represent the various subject
% groups. The numeric data consists of normalized percentages of hits in a
% given COG category for a given subject group.

get(dm2)

%%
% For each main COG category, we compute a cumulative normalized percentage
% and store the results in a new DataMatrix object named |dm2Count|.

codes = {'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', ...
    'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V'}; % COG code to consider

n = numel(codes);
N = size(dm2,2);
count = zeros(n,N);

for i = 1:n
    try
        count(i,:) = sum(dm2.(codes{i}));
    catch
        sprintf('COG code %s is not found in the data set.',codes{i});
	end
end

dm2Count = bioma.data.DataMatrix(count, codes, dm2.ColNames);    

%%
% To investigate whether the COG enrichment patterns are different among
% the three age-related groups, we first consider the data associated with
% the adult, children and infant subjects.

group1 = {'Adult', 'Child', 'Infant'};
figure(); plot(dm2Count.(':')(group1), '.-', 'LineWidth', 2);
haxis = gca;
set(haxis, 'XTick', 1:n);
set(haxis, 'XTickLabel', codes);
legend(group1, 'location', 'northwest')
xlabel('COG categories');
ylabel('Normalized percentage of assigned genes');

%%
% We observe from this plot that adult subjects and children appear to have
% a similar pattern of enrichment in terms of COG categories. The infant
% subjects, on the other hand, display some singularities for categories G,
% K and L, corresponding to carbohydrate transport and metabolism,
% transcription, and replication respectively.
%
% In light of this affinity between adult and child microbiome functional
% patterns, we consider a combination of the two samples (Adult+Child) when
% performing a comparison against other environmental sample microbiomes.

group2 = {'Adult+Child', 'Soil', 'Whale Fall Ave.', 'Sargasso'};
figure(); plot(dm2Count.(':')(group2), '.-', 'LineWidth', 2);
haxis = gca;
set(haxis, 'XTick', 1:n, 'XTickLabel', codes);
legend(group2, 'location', 'north');
xlabel('COG categories');
ylabel('Normalized percentage of assigned genes');

%%
% The most striking differences between the human microbiome enrichment
% pattern and those of other environmental microbial communities is related
% to COG category G (carbohydrate metabolism). This is perhaps related to
% the notion that the colonic microbiota utilizes otherwise indigestible
% polysaccharides and peptides as major resource for energy production and
% biosynthesis of cellular components. The enrichment of several enzymes
% for DNA repair is also noteworthy (COG category L).

%%
% A more effective way of visualizing the distribution of patterns of
% COG-assigned genes between each type of microbiome consists of plotting
% the enrichment values for each COG category along a circumference. For
% each data point, the distance from the center of the circle is
% proportional to the enrichment value.

r = dm2Count.(':')(group2);
%colors = hsv(numel(group2));
colors = {'b', 'g', 'r', 'k'};
theta = (linspace(0,2*pi,n+1))';
figure();
hold on;

for i = 1:numel(group2)
    rho = [r(:,i); r(1,i)];
    plot(rho .* cos(theta), rho .* sin(theta), '-', 'Color', colors{i}, 'LineWidth', 2);
end 
legend(group2, 'location', 'NorthEastOutside')

%=== plot outside circle and labels
m = max(max(r));
for i = 1:n
    text( (m + .5) * cos(theta(i)), (m + .5) * sin(theta(i)), codes{i}, ...
        'HorizontalAlignment', 'center');
end

theta = (linspace(0,2*pi,100))';
plot(m * cos(theta), m * sin(theta), 'k-');

axis('equal')
axis([-1 1 -1 1] * (m+1))
axis('off')

%% Clustering Microbiomes Based on Their Functional Profiles
% We can examine the relationship between the human gut microbiomes and
% other environmental microbiomes using the enrichment values for each COG.
% We create a hierarchical cluster tree using the complete linkage
% algorithm and the distance matrix generated by considering the
% correlation between data points. The samples considered include: Adult,
% Child, Infant, American, Soil, Whale fall (1, 2, and 3) and Sargasso.

group3 = {'Adult', 'Child', 'Infant', 'American (SubjF)', 'American (SubjM)', ...
    'Soil', 'Whale Fall 1', 'Whale Fall 2', 'Whale Fall 3', 'Sargasso'};

z = linkage((dm2Count.(':')(group3))', 'complete', 'correlation');
dendrogram(z, 'orientation', 'left', 'labels', group3, 'colorthreshold', 'default')

%%
% The clustering analysis further shows that, while the adult and child
% microbiomes present similar profiles, those of infants have a distinct
% profile. Furthermore, some differences can be observed between the
% Japanese individuals and the American subjects. Finally, as expected, the
% human gut microbiome appears to be specific of the human species, and not
% related to the other environmental microbial communities.

%% References
% [1] Gill, S., Pop, M., DeBoy, R., Eckburg, P., Turnbaugh, P., Samuel, B.,
% Gordon, J., Relman, D., Fraser-Liggett, C., Nelson, K. (2006).
% Metagenomic Analysis of the Human Distal Gut Microbiome, Science, 312,
% 1355-1359.
%
% [2] Kurokawa, K., Itoh, T., Kuwahara, T., Oshima, K., Toh, H., Toyoda,
% A., Takami, H., Morita, H., Sharma, V., Shrivastava, T., Taylor, T.,
% Noguchi, H., Mori, H., Ogura, Y., Ehrlich, D., Itoh, K., Takagi, T.,
% Sakaki, Y., Hayashi, T., Hattori, M. (2007). Comparative Metagenomics
% Revealed Commonly Enriched Gene Sets in Human Gut Microbiomes, DNA
% Research, 14, 168-181.

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20GUTMICROBIOMEDEMO%20in%20Bioinformatics%20Toolbox%204.2 Provide feedback for this example.>*

displayEndOfDemoMessage(mfilename)
