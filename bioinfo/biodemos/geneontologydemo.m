%% Gene Ontology Enrichment in Microarray Data
% This example shows how to enrich microarray gene expression data using 
% the Gene Ontology relationships.

%   Copyright 2005-2013 The MathWorks, Inc. 


%% Introduction
% Gene Ontology is a controlled method for describing terms related to
% genes in any organism. As more gene data is obtained from organisms, it
% is annotated using Gene Ontology. Gene Ontology is made of three smaller
% ontologies or aspects: Molecular Function, Biological Process, and
% Cellular Component. Each of these ontologies contains terms that are
% organized in a directed acyclic graph with these three terms as the
% roots. The roots are the broadest terms relating to genes. Terms further
% away from the roots get more specific. For this example you will use
% microarray data from the <yeastdemo.html Gene Expression Profile
% Analysis> example to look at the significance of interesting genes and
% Gene Ontology terms that are associated with the microarray probes. More
% specifically, you will further investigate to determine if a set of genes
% that cluster together are also involved in a common molecular function.

%% Examples Using Gene Ontology Functions
% The Gene Ontology database is loaded into a MATLAB(R) object using the
% Bioinformatics Toolbox |geneont| function.

GO = geneont('live',true); % this step takes a while 
get(GO)

%%
% Every Gene Ontology term has an accession number which is a seven digit
% number preceded by the prefix 'GO:'. To look at a specific term you first
% create a sub-ontology by sub-scripting the object and then inspect the
% 'terms' property. A hash-table is implemented in the GO object for
% efficiently looking up term IDs.

GO(5840).terms % the ribosome Gene Ontology term

%%
% This term represents the 'ribosome', and it has fields for |is_a| and
% |part_of|. These fields represent the relationships between Gene Ontology
% terms. Gene Ontology terms can be seen as nodes in an acyclic graph. You
% can traverse such relationships with the methods |getancestors|,
% |getdescendants|, |getrelatives|, and |getmatrix|. For example, the
% |getancestors| method returns any less specific term than 'ribosome'
% (i.e., its parents in the graph). 

ancestors = getancestors(GO,5840)
riboanc = GO(ancestors)

%% 
% To visualize these relationships we use the |biograph| function and the
% |getmatrix| method from the Bioinformatics Toolbox. The |getmatrix|
% method returns a square matrix of relationships of the given Gene
% Ontology object. This graph is sometimes called an 'induced' graph. 

cm = getmatrix(riboanc);
BG = biograph(cm,get(riboanc.Terms,'name'))
view(BG)

%% Using Clustering to Select an Interesting Subset of Genes
% To show how Gene Ontology information is useful, we will look at
% microarray data from the <yeastdemo.html Gene Expression Profile
% Analysis> example. This data has 6400 genes on the microarray that are
% involved with many different aspects of yeast gene expression. A small
% portion of these might show interesting behavior for this microarray
% experiment. We will use the Gene Ontology to better understand if and how
% these genes are related in the cell. The full yeast data can be found at
% the NCBI <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28
% Website>. 

load yeastdata
whos yeastvalues genes

%% 
% The example "Gene Expression Profile Analysis" shows several ways to
% cluster the data from the experiment. In this example, K-means clustering
% is used to select a group of about 240 genes for study.

%%
% First, the data needs cleaning up. There are some empty spots on the
% chip, and some genes have missing values. In this example the empty spots
% are removed from the data set and the |knnimpute| function imputes the
% missing values (marked with NaNs).

% Remove data for empty spots
emptySpots = strcmp('EMPTY',genes);
yeastvalues = yeastvalues(~emptySpots,:);
genes = genes(~emptySpots);
fprintf('Number of genes after removing empty spots is %d.\n',numel(genes))

% Impute missing values
yeastvalues = knnimpute(yeastvalues);

%%
% Next, the function |genelowvalfilter| removes genes with low overall
% expression. In this example a fairly large value is used to filter the
% genes. The <yeastdemo.html Gene Expression Profile Analysis> example
% shows alternative filtering techniques.

mask = genelowvalfilter(yeastvalues,genes,'absval',log2(3.5));
highexpIdx = find(mask);
yeastvalueshighexp = yeastvalues(highexpIdx,:);
fprintf('Number of genes with high expression is %d.\n',numel(highexpIdx))

%%
% The |kmeans| function from the Statistics Toolbox(TM) groups the data into
% four clusters using correlation as the distance metric. 

rand('state',0);
[cidx, ctrs] = kmeans(yeastvalueshighexp, 4, 'dist','corr', 'rep',20);
for c = 1:4
    subplot(2,2,c);
    plot(times,yeastvalueshighexp((cidx == c),:)'); 
    axis('tight');
end
suptitle('K-Means Clustering of Profiles');

%%
% The plots show four fairly different clusters. By looking at the
% centroids of the clusters you can see clearly how they differ.   
figure
for c = 1:4
    subplot(2,2,c);
    plot(times,ctrs(c,:)','linewidth',4,'color','k');
    axis tight
    axis off    % turn off the axis
end
suptitle('Centroids of K-Means Clustering of Profiles');

%%
% The first cluster in the top left corner represents the genes that are
% up-regulated with their expression levels falling off a little in the
% final chip. The genes in this cluster will be the subset used for the
% remainder of this experiment.

clusterIdx = highexpIdx(cidx==1);
fprintf('Number of genes in the first cluster is %d.\n',numel(clusterIdx))

%% Getting Annotated Genes from the Saccharomyces Genome Database
% Many Genome Projects interact with the Gene Ontology Consortium when
% annotating genes. Gene annotations for several organisms can be found at
% the <http://www.geneontology.org/GO.current.annotations.shtml Gene
% Ontology Website>. In addition, annotations for individual organisms can
% be found at their respective websites (such as the
% <http://www.yeastgenome.org/ Yeast Genome database>). These annotations
% are updated frequently and are usually curated by members of the genome
% group for each organism. NCBI also has a collective list of gene
% annotations that relate to their
% <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene Entrez Gene 
% database>. These annotation files consist of large lists of genes and
% their associated Gene Ontology terms. These files follow the structure
% defined by the <http://www.geneontology.org/GO.annotation.shtml#file Gene
% Ontology Consortium>. The function |goannotread| will parse these
% uncompressed files and put the information into a MATLAB structure. The
% file |yeastgenes.sgd| was obtained from the Gene Ontology Annotation
% site.  

%%
% For this analysis you will look at genes that are annotated as molecular
% function (i.e., the 'Aspect' field is set to 'F'). However, you could
% extend this analysis to see if the clustered genes are involved in common
% processes ('P') or are co-located in the same cell compartment ('C'). The
% fields that are of interest are the gene symbol and associated ID. In GO
% Annotation files, these have field names DB_Object_Symbol and GOid,
% respectively. Create a map for efficient search of the Gene Symbol.
% Observe that not every gene from the 6314 genes on the microarray is
% annotated.

SGDann = goannotread('yeastgenes.sgd','Aspect','F',...
                     'Fields',{'DB_Object_Symbol','GOid'}); 

SGDmap = containers.Map();
for i=1:numel(SGDann)
    key = SGDann(i).DB_Object_Symbol;
    if isKey(SGDmap,key)
        SGDmap(key) = [SGDmap(key) SGDann(i).GOid];
    else
        SGDmap(key) = SGDann(i).GOid;
    end
end

fprintf('Number of annotated genes related to molecular function is %d.\n',SGDmap.Count)
fprintf('Number of unique GO terms associated to annotated genes is %d.\n',numel(unique([SGDann.GOid])))
fprintf('Number of gene-GO term associations is %d.\n',numel(SGDann))


%% Counting Annotated Genes From the Microarray
% We will keep a tally of the number of genes that are annotated for every
% Gene Ontology term. At the same time, we will keep track of how many Gene
% Ontology terms have interesting genes associated with them from above.
% You can use this information to determine how often Gene Ontology terms are
% represented in the microarray experiment.

m = GO.Terms(end).id;           % gets the last term id
geneschipcount = zeros(m,1);    % a vector of GO term counts for the entire chip.
genesclustercount = zeros(m,1); % a vector of GO term counts for interesting genes.
for i = 1:numel(genes)
    if isKey(SGDmap,genes{i})
        goid = getrelatives(GO,SGDmap(genes{i}));
        % update vector counts
        geneschipcount(goid) = geneschipcount(goid) + 1;
        if (any(i == clusterIdx))
           genesclustercount(goid) = genesclustercount(goid) + 1;
        end
    end
end

%% Looking at Probability of Gene Ontology Annotation  
% You can find the most significant annotated terms by looking at
% the probabilities that the terms are counted by chance. For this you
% can use the hypergeometric probability distribution function (|hygepdf|).
% This function returns the p-value associated to each term, you can create
% a list of the most significant GO terms by ordering the p-values. 

pvalues = hygepdf(genesclustercount,max(geneschipcount),...
                  max(genesclustercount),geneschipcount);
[dummy idx] = sort(pvalues); 

% create a report
report = sprintf('GO Term      p-val  counts  definition\n');
for i = 1:10
    term = idx(i);
    report = sprintf('%s%s\t%-1.4f\t%-d / %-d\t%s...\n', report, ...
                    char(num2goid(term)), pvalues(term),...
                    genesclustercount(term),geneschipcount(term),...
                    GO(term).Term.definition(2:min(end,60)));
end
disp(report);

%% Further Analysis of the Most Significant Terms
% You can use the methods described earlier in this example to find out
% more about the terms that appear high on this list.

%%
% First look at the ancestors of the top item on the list.

topItem = idx(1);
GO(topItem).terms % the most significant gene
topItemAncestors = getancestors(GO,topItem)

%%
% Note that the term 16616 appears as one of the ancestors. If you look at
% the list for the second item, you will see many of the same ancestors.

secondItem = idx(2);
GO(secondItem).terms % the second most significant gene
secondItemAncestors = getancestors(GO,secondItem)

%%
% You can now build a sub-ontology that includes the ancestors of the ten
% (for example) most significant terms and visualize this using the
% |biograph| function.  

subGO = GO(getancestors(GO,idx(1:10)))
[cm acc rels] = getmatrix(subGO);
BG = biograph(cm,get(subGO.Terms,'name'))

%%
% Use the p-values, calculated before, to assign a color to the graph
% nodes. In this example an arbitrary color map is used, where bright red
% is the most significant and bright green is the least significant.

for i=1:numel(acc) 
    pval = pvalues(acc(i));
    color = [(1-pval).^(10),pval.^(1/10),0.3];
    set(BG.Nodes(i),'Color',color); 
    set(BG.Nodes(i),'Label',num2str(acc(i))) % add info to datatips
end

view(BG);

%% References
% [1] http://www.geneontology.org/
%%
% [2] http://www.yeastgenome.org/
%%
% [3] Gentleman, R. 'Basic GO Usage'. Bioconductor vignette May 16, 2005
% http://bioconductor.org/docs/vignettes.html
%%
% [4] Gentleman, R. 'Using GO for Statistical Analyses'. Bioconductor 
% vignette May 16, 2005 http://bioconductor.org/docs/vignettes.html

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20GENEONTOLOGYDEMO%20in%20Bioinformatics%20Toolbox%204.2 Provide feedback for this example.>*

displayEndOfDemoMessage(mfilename)
