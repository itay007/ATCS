%% Performing a Metagenomic Analysis of a Sargasso Sea Sample
% This example illustrates a simple metagenomic analysis on a sample data
% set from the Sargasso Sea. It requires the taxonomy information included
% in the files |gi_taxid_prot.dmp|, |names.dmp| and |nodes.dmp| (see the
% compressed file |taxdump|), which you can download from the
% <ftp://ftp.ncbi.nih.gov/pub/taxonomy/ NCBI taxonomy FTP site>. 

%   Copyright 2007-2012 The MathWorks, Inc.

%% Introduction
% Metagenomics is the study of the taxonomic composition of a sample of
% organisms obtained from a common habitat. It usually consists of the
% comparison of the sequence samples against databases of known sequences
% and the use of taxonomy information to classify the sample species. The
% main goals of a metagenomic analysis include the quantification of the
% relative abundance of known species and the identification of unknown
% sequences for which no relatives have yet been identified.

%% Reading BLASTX Hit Report
% In this example, we consider a small subset (100 reads) of the Sargasso
% Sea data set [1], which has been searched against the NCBI-NR database
% using BLASTX with default parameters. For convenience, the resulting
% BLAST report has been saved and compressed into the file
% |sargasso-sample1-100.rpt.gz|, and it is provided with 
% Bioinformatics Toolbox(TM). We read the report content and extract
% relevant information such as the high-scoring pairs, their score,
% expectation value and percent identity.

%=== open the blastx report
reportFilename = gunzip('sargasso-sample1-100.rpt.gz',tempdir);
fid = fopen(reportFilename{1}, 'rt');

%=== read all strings to be able to write into xls
blastInfo = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s %s'); 
fclose(fid);
delete(reportFilename{1});

%=== extract relevant information
queries = blastInfo{1};
hits = blastInfo{2};
ident = str2double(blastInfo{3});
evalue = str2double(blastInfo{11});
score = str2double(blastInfo{12});

numEntries = numel(queries)

%% Filtering BLAST Hits
% Because we are interested only in significant hits, we filter the results
% based on their score, expectation value and percent identity with the
% query sequences. By using this filtering process, we reduce the
% number of hits to approximately one quarter of the original hits.
 
%=== setup filter criteria
scoreThreshold = 100;
evalueThreshold = 10^-5;
identThreshold = 50;

%=== consider only hits satisfying the criteria
k = find(score > scoreThreshold & evalue < evalueThreshold & ident > identThreshold);
queries = queries(k);
hits = hits(k);
evalue = evalue(k);
score = score(k);

numEntries = length(k)

%=== clear report
clear blastInfo

%% Memory-Mapping the Taxonomy Data File
% The taxonomic classifications for all GenBank(R) sequences are
% stored in large files that are updated weekly as new sequences are
% submitted and the taxonomic information is refined. 
% To retrieve this information in a quick and efficient way,
% we create a map between any possible gi number in the GenBank database
% and its associated taxonomic identifier (taxid). Because currently there
% are more than 100 million live gi numbers, the memory requirements for
% loading such a large data set can be very demanding. Thus, using the
% provided helper function |mapTaxoFile|, we read the data in blocks of
% 1MB, save it as a binary file and then use the function |memmapfile| to
% map into memory the content of the file itself, so that the data can be
% accessed using standard indexing operations. See |memmapfile| help for
% more details.  

taxoFilenameIn = 'gi_taxid_prot.dmp';
taxoFilenameOut = 'gi_taxid_prot_map.dmp';

%=== create map so that gi --> taxid, taxid = -1 if no live gi
blockSize = 2^20; % block size (1MB)
mapTaxoFile(taxoFilenameIn, taxoFilenameOut, blockSize);

%=== map file into memory
mt = memmapfile(taxoFilenameOut, 'format', 'int32'); 

%%
% We can access the taxid of first ten live GenBank sequences as follows:

q = find(mt.Data(1:100)>0);
mt.Data(q(1:10))
clear q

%% Annotating the BLAST Report with Taxonomic Information
% We are now interested in performing a taxonomic annotation of each hit in
% the BLAST report. We extract the gi number of each hit and
% retrieve its associated taxid. 

%=== extract gi number for each hit
gi = zeros(1, numEntries);
for i = 1:numEntries
    g = str2double(regexpi(hits{i}, '(?<=gi\|)\d+', 'match', 'once'));
    if ~isempty(g)
        gi(i) = g;
    end
end 

%=== determine taxid for each hit
taxid = mt.Data(gi);

%%
% If you performed the BLAST search against a database that is outdated
% with respect to the taxonomy information included in the |nodes.dmp|
% file, some gi numbers might be superseded. Therefore, you need to exclude
% from the analysis those sequences associated with superseded entries.

%=== ignore dead gi numbers
livegi = (taxid > 0);
gi = gi(livegi); taxid = taxid(livegi);
queries = queries(livegi); hits = hits(livegi); 
evalue = evalue(livegi); score = score(livegi);


%%
% During the search against the NCBI-NR Database, the first query
% (SHAA001TR) hit  |n| sequences with significant expectation value and
% score. We can look at the taxonomic assignment of these hits using the
% array |taxid|.

SHAA001TR = strcmp('SHAA001TR', queries);
n = sum(SHAA001TR)
hits(SHAA001TR)
taxid(SHAA001TR)

%% Classifying BLAST Hits by Scientific Name
% Every taxid corresponds to a specific taxon, which has been given a
% scientific name and possibly various synonyms. For our classification
% purposes, we are interested in the scientific names only. Thus, we extract
% this information and annotate each BLAST hit in the report using the
% scientific names, rather than the taxids. 

%=== read taxonomy name file
taxonomyFilenameIn = 'names.dmp';
fid1 = fopen(taxonomyFilenameIn,'rt');
nameInfo = textscan(fid1, '%d%s%s%s', 'delimiter', '|'); 
fclose(fid1);

%=== preallocate space for SN
maxTaxid = max(double(nameInfo{1}));
SN = repmat({''}, maxTaxid, 1);

%=== populate array so that taxid --> scientific name
ind = strncmp('scientific',nameInfo{4},10); % indices of scientific names in the array
SN(nameInfo{1}(ind)) = strtrim(nameInfo{2}(ind)); 

%=== assign name to every hit
sciNames = SN(taxid);

%%
% We can look at the scientific names of the organisms whose sequences were hit by the
% first query by considering the first |n| elements in the array |sciNames|, as follows: 

sciNames(1:n)

%% Saving Annotated BLAST Report
% Once we determine the taxonomic classification for each hit, we can
% include the information in a text file as shown below: 

%=== create annotated report for first n hits
textFilename = 'sargasso-annotated-report.txt';
fid = fopen(textFilename, 'wt');
for i = 1:n
    fprintf(fid, '%s\t%s\t%d\t%d\t%s\n', queries{i}, hits{i}, evalue(i), taxid(i), sciNames{i});
end
fclose(fid);

type sargasso-annotated-report.txt

%% Determining the Taxonomic Distribution of BLAST Hits
% One reason to classify the sequence hits in a BLAST report is to study
% their taxonomic distribution. We can easily create a list of organisms
% that are represented in the report, their taxids and their frequency as
% follows:

%=== distribution by taxid
taxidList = unique(taxid);       % list of unique taxids
T = accumarray(taxid, 1);        % multiplicity of taxids
taxidCount = T(unique(taxid));   % number of hits for each taxon

%=== simple statistics of the hit distribution
numTaxa = length(taxidList)          % number of distinct taxa
[maxCount,maxInd] = max(taxidCount); % most represented taxon
maxTaxid = taxidList(maxInd)         % taxid of the most represented taxon
maxSN = SN(maxTaxid)                 % name of the most represented taxon
maxCount

%%
% From the simple statistics on the taxonomic distribution, we observe that
% the most represented taxon is the _Burkholderia sp.383_ (taxid 269483).
% The over-representation of this bacterium, which is usually found in
% terrestrial settings, in the Sample 1 of the Sargasso Sea data set is
% discussed in [1].

%% Filtering Out Isolated Assignments
% Several taxa in the report appear to be isolated assignments because they
% are hit by only one sequence. These taxa are rarely true members of the
% environmental community under investigation. Thus it is useful to
% identify them and discard them if needed.

t1 = taxidCount == 1;
isolated = length(find(t1))

taxidList = taxidList(~t1);
taxidCount = taxidCount(~t1);

numTaxaFiltered = length(taxidCount)

%% Plotting Taxonomic Distribution of BLAST Hits
% If we plot the taxonomic distribution of the hits on a bar chart, we
% observe that the majority of taxa has a low number of occurrences. 

%=== plot by sorting the counts
hFig = figure(); bar(sort(taxidCount));
xlabel('Distinct taxonomic assignments');
ylabel('Number of hits');
title('Taxonomic distribution of filtered hits');
set(gca, 'XTickLabel','');

%% Limiting the Analysis to the Best Hit for Each Query
% We can repeat the above procedure by limiting the analysis to only the
% best-scoring hit for each query sequence. Even though analyses limited to
% the best-scoring hits cannot depict a complete and accurate picture of
% the situation, they can be useful as a first approximation and overcome
% the difficulty inherent with large data sets.

%== get only best hits
[queriesUnique, idx] = unique(queries, 'first');      % best hits rows
bestHitTaxid = taxid(idx);
bestHitSciName = sciNames(idx);

%=== count occurrences
T = accumarray(bestHitTaxid, 1);        % multiplicity of taxids
bestCount = T(unique(bestHitTaxid));    % number of hits for each taxon
bestCountNames = SN(unique(bestHitTaxid)); 

%=== five most represented taxa
[bestCountSorted, idx] = sort(bestCount, 'descend');
bestCountSorted(1:5)
bestCountNames(idx(1:5))

%%
% In our example, when only the best-scoring hits are considered,
% _Burkholderia_, _Candidatus pelagibacter ubique_ and _Shewanella_
% appear to be the most represented taxa in the report. While finding
% _Candidatus pelagibacter ubique_ is not surprising, because it is a
% dominant form of life in the Sargasso Sea, _Burkholderia_ and
% _Shewanella_ are not expected to be present in this marine sample where
% nutrients and resources are low, because they live either in terrestrial
% settings or in aquatic, nutrient-rich environments respectively. For a
% detailed discussion regarding the presence of these bacteria in the
% Sargasso Sea, see [1].
%

%% Memory-Mapping the Taxon Node Information
% Oftentimes, to gain a clear vision of the taxonomic distribution
% of a sequence set, Linnaean categories higher than species are
% considered. To perform this analysis, we need to create a map between
% each taxid and its assigned rank, as well as a map between each taxid and
% the taxid of its parent node, according to the NCBI Taxonomy Database
% schema. Files containing this information can be created with the helper
% function |mapNodeFile|.

nodeFilename   = 'nodes.dmp';
parentFilename = 'nodes_parent_map.dmp';
rankFilename   = 'nodes_rank_map.dmp';

%=== create a map 
mapNodeFile(nodeFilename, parentFilename, rankFilename, blockSize);

%=== map the files into memory
mmParentObj = memmapfile(parentFilename, 'format', 'int32'); % taxid --> taxid_parent
mmRankObj = memmapfile(rankFilename, 'format', 'int32');     % taxid --> rank
           
%% Classifying BLAST Hits by Higher Taxonomic Rank
% After the maps are created, for every hit that is associated with a taxid
% corresponding to a Linnaean category more specific than the target rank,
% we determine the parental taxid and its rank until the target is reached.
% Then we annotate the hit with the taxid of its more distant ancestor.
% Synthetic constructs or nodes with no rank are considered descendants of
% the root. This procedure of walking up the taxonomic hierarchy is
% performed by the helper function |findTaxoRank|.

%%
% Suppose we are interested in classifying our hits according to the
% superkingdom to which they belong. After assigning the superkigdom taxid
% to each hit, we group and count the occurrences as follows:
 
%=== find superkingdom assignments
skRank = findTaxoRank(taxidList, mmRankObj, mmParentObj, 1);
sk = accumarray(skRank, 1);
skCount = sk(unique(skRank));
skNames = SN(unique(skRank));

%=== plot pie chart
hFig = figure(); 
pie(skCount); colormap(summer)
legend(skNames, 'location', 'EastOutside');

%% 
% As expected, the majority of the hits are bacteria. 
% Similarly, we can determine the taxonomic distribution at the level of phylum, class,
% order and family as shown below:

rTargetString = {'phylum', 'class', 'order', 'family'}
rTarget = [5 8 11 14];

numTarget = numel(rTarget);
rank = cell(1,numTarget);

%=== annotate hits with the taxid at the target level
for i = 1:numTarget
    rank{i} = findTaxoRank(taxidList, mmRankObj, mmParentObj, rTarget(i));
end

%=== determine the distribution
count = cell(1,numTarget);
names = cell(1,numTarget);

for i = 1:numTarget
    list = unique(rank{i});
    T = accumarray(rank{i}, 1);
    count{i} = T(list);
    names{i} = SN(list);
end

%=== plot the first two classifications
for i = 1:2
    figure(); barh(count{i});
    set(gca, 'yTick', 1:numel(names{i}));
    set(gca, 'yTicklabel', names{i});
    xlabel('Occurrences');
    title(['Taxonomic distribution at the ' rTargetString{i} ' level'])
end

%=== Draw a Pareto chart for the phyla
pnames = names{1}; pcount = count{1};
[ppeaks, pind] = sort(pcount, 'descend'); 
plabels = pnames(pind);
figure(); pareto(pcount, pnames);
ylabel('Occurrences');
pticks = get(gca, 'XTick'); 
ht = text(pticks', ppeaks(pticks)+10, plabels(pticks), 'rotation', 45);
title('Pareto chart for distribution at the phylum level');
set(gca, 'XtickLabel', '')

%% Representing the Taxonomic Distribution on a Graph
% The taxonomic distributions at different levels are related to each other
% by hierarchy. Suppose we want to look at the distribution of hits across
% phyla and visualize them on a graph. After filtering out the counts of the
% low represented phyla (<5 counts), we create a connectivity matrix where
% all phyla are direct children of the root.

k = count{1} > 5;
phylaNames = names{1}(k);
n1 = length(phylaNames);
CM = zeros(n1); CM(1,2:end) = 1;
bg = biograph(CM, phylaNames);
view(bg)

%% 
% We can now consider all the hits classified as Proteobacteria (taxid
% 1224), and perform the same distribution analysis at the level of classes.

%=== consider only Proteobacteria
pb = taxidList(rank{1} == 1224);
pbRank = findTaxoRank(pb, mmRankObj, mmParentObj, 8);
pbList = unique(pbRank);
pbT = accumarray(pbRank, 1);
pbCount = pbT(pbList);
pbNames = SN(pbList);
 
%=== filter out if less than 5 counts
h = pbCount > 5;
pbCount = pbCount(h);
n2 = length(pbCount);
pbNames = pbNames(h)

%%
% To represent the various phyla and the Proteobacteria class in the same graph,
% we need to create a connectivity matrix such that all the phyla are
% children of the root and all the Proteobacteria classes are children of
% the Proteobacteria node. In the graph, the labels include the class names
% and the number of occurrences in the BLAST report.

%=== find Proteobacteria node
x = strcmp('Proteobacteria', phylaNames);

%=== combine names and counts
numNodes = n1 + n2;

allNames(1:n1) = phylaNames;
allNames(n1+1:numNodes) = pbNames;

allCount(1:n1) = count{1}(k);
allCount(n1+1:numNodes) = pbCount;

%=== create labels for nodes (scientific name and count)
labels = cell(1,numNodes);
for node = 1:numNodes
    labels{node} = [allNames{node} ' (' num2str(allCount(node)) ')'];
end

%=== create graph
CM = zeros(numNodes);
CM(1,2:n1) = 1; CM(x, n1+1:numNodes) = 1; CM(x,x) = 0;
bg = biograph(CM, labels, 'showArrows', 'off');
bg.view

%=== clear memory mapped variables
clear mmParentObj mmRankObj mt

%% References
% [1] Venter, J.C., et al. "Environmental genome shotgun sequencing of the
%     Sargasso sea", Science, 304, 66-74, 2004.
%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20METAGENOMICDEMO%20in%20Bioinformatics%20Toolbox%204.2 Provide feedback for this example.>*

displayEndOfDemoMessage(mfilename)
