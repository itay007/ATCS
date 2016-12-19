%% Exploring Gene Expression Data
% This example shows how to identify differentially expressed genes. Then
% it uses Gene Ontology to determine significant biological functions that
% are associated to the down- and up-regulated genes.
 
%   Copyright 2006-2012 The MathWorks, Inc.


%% Introduction
% Microarrays contain oligonucleotide or cDNA probes for comparing the
% expression profile of genes on a genomic scale. Determining if changes in
% gene expression are statistically significant between different
% conditions, e.g. two different tumor types, and determining the
% biological function of the differentially expressed genes, are important
% aims in a microarray experiment.
%
% A publicly available dataset containing gene expression data of 42 tumor
% tissues of the embryonal central nervous system (CNS, Pomeroy et al.
% 2002) is used for this example. The samples were hybridized on
% Affymetrix(R) HuGeneFL GeneChip(R) arrays.
% 
% The CNS dataset (CEL files) is available at the
% <http://www.broad.mit.edu/mpr/CNS/ CNS experiment> web site. The 42 tumor
% tissue samples include 10 medulloblastomas, 10 rhabdoid, 10 malignant
% glioma, 8 supratentorial PNETS, and 4 normal human cerebella. The CNS raw
% dataset was preprocessed with the Robust Multi-array Average (RMA) and GC
% Robust Multi-array Average (GCRMA) procedures. For further information on
% Affymetrix oligonucleotide microarray preprocessing, see
% <affypreprocessdemo.html Preprocessing Affymetrix Microarray Data at the
% Probe Level>.
% 
% You will use the t-test and false discovery rate to detect differentially
% expressed genes between two of the tumor types. Additionally, you will
% look at Gene Ontology terms related to the significantly up-regulated
% genes.

%% Loading the Expression Data
% Load the MAT file |cnsexpressiondata| containing three DataMatrix
% objects. Gene expression values preprocessed by RMA and GCRMA (MLE and
% EB) procedures are stored in the DataMatrix objects |expr_cns_rma|,
% |expr_cns_gcrma_mle|, and |expr_cns_gcrma_eb| respectively.  
load cnsexpressiondata

%%
% In each DataMatrix object, each row corresponds to a probe set on the
% HuGeneFl array, and each column corresponds to a sample. The row names
% are the probe set IDs and column names are the sample names. The
% DataMatrix object |expr_cns_gcrma_eb| will be used in this example. You
% can use either of the other two expression variables as well.
%%
% You can get the properties of the DataMatrix object |expr_cns_gcrma_eb|
% using the |get| command.
get(expr_cns_gcrma_eb)

%%
% Determine the number of genes and number of samples in the DataMatrix
% object |expr_cns_gcrma_eb|.
[nGenes, nSamples] = size(expr_cns_gcrma_eb)

%%
% You can use gene symbols instead of the probe set IDs to label the
% expression values. The gene symbols for the HuGeneFl array are provided
% in a MAT file containing a |Map| object. 
load HuGeneFL_GeneSymbol_Map;

%%
% Create a cell array of gene symbols for the expression values from the
% |hu6800GeneSymbolMap| object.
huGenes = values(hu6800GeneSymbolMap, expr_cns_gcrma_eb.RowNames);

%%
% Set the row names of the |exprs_cns_gcrma_eb| to gene symbols using the
% |rownames| method of the DataMatrix object.
expr_cns_gcrma_eb = rownames(expr_cns_gcrma_eb, ':', huGenes);

%% Filtering the Expression Data
% Remove gene expression data with empty gene symbols. In the example, the
% empty symbols are labeled as |'---'|.
expr_cns_gcrma_eb('---', :) = [];

%% 
% Many of the genes in this study are not expressed, or have only small
% variability across the samples. Remove these genes using non-specific
% filtering.

%%
% Use |genelowvalfilter| to filter out genes with very low absolute
% expression values.
[mask, expr_cns_gcrma_eb] = genelowvalfilter(expr_cns_gcrma_eb);

%%
% Use |genevarfilter| to filter out genes with a small variance across
% samples.
[mask, expr_cns_gcrma_eb] = genevarfilter(expr_cns_gcrma_eb);
%%
% Determine the number of genes after filtering.
nGenes = expr_cns_gcrma_eb.NRows 

%% Identifying Differential Gene Expression
% You can now compare the gene expression values between two groups of
% data: CNS medulloblastomas (MD) and non-neuronal origin malignant gliomas
% (Mglio) tumor.

%%
% From the expression data of all 42 samples, extract the data of the 10
% MD samples and the 10 Mglio samples.
MDs = strncmp(expr_cns_gcrma_eb.ColNames,'Brain_MD', 8);
Mglios = strncmp(expr_cns_gcrma_eb.ColNames,'Brain_MGlio', 11);

%%
MDData = expr_cns_gcrma_eb(:, MDs);
get(MDData)
%%
MglioData = expr_cns_gcrma_eb(:, Mglios);
get(MglioData)

%%
% A standard statistical test for detecting significant changes between the
% measurement of a variable in two groups is the t-test. Conduct a t-test
% for each gene to identify significant changes in expression values
% between the MD samples and Mglio samples. You can inspect the test
% results from the normal quantile plot of t-scores and the histograms of
% t-scores and _p-values_ of the t-tests. 
[pvalues, tscores] = mattest(MDData, MglioData,...
                                'Showhist', true', 'Showplot', true);

%%
% In any test situation, two types of errors can occur, a false positive by
% declaring that a gene is differentially expressed when it is not, and a
% false negative when the test fails to identify a truly differentially
% expressed gene. In multiple hypothesis testing, which simultaneously
% tests the null hypothesis of thousands of genes using microarray
% expression data, each test has a specific false positive rate, or a false
% discovery rate (FDR). False discovery rate is defined as the expected
% ratio of the number of false positives to the total number of positive
% calls in a differential expression analysis between two groups of samples
% (Storey et al., 2003). 

%%
% In this example, you will compute the FDR using the Storey-Tibshirani
% procedure (Storey et al., 2003). The procedure also computes the q-value
% of a test, which measures the minimum FDR that occurs when calling the
% test significant. The estimation of FDR depends on the truly null
% distribution of the multiple tests, which is unknown. Permutation methods
% can be used to estimate the truly null distribution of the test
% statistics by permuting the columns of the gene expression data matrix
% (Storey et al., 2003, Dudoit et al., 2003). Depending on the sample size,
% it may not be feasible to consider all possible permutations. Usually a
% random subset of permutations are considered in the case of large sample
% size. Use the |nchoosek| function in Statistics Toolbox(TM) to find out
% the number of all possible permutations of the samples in this example.
all_possible_perms = nchoosek(1:MDData.NCols+MglioData.NCols, MDData.NCols);
size(all_possible_perms, 1)

%%
% Perform a permutation t-test using |mattest| and the |PERMUTE| option to
% compute the _p-values_ of _10,000_ permutations by permuting the columns
% of the gene expression data matrix of MDData and MglioData (Dudoit et
% al., 2003).
pvaluesCorr = mattest(MDData, MglioData, 'Permute', 10000);

%%
% Determine the number of genes considered to have statistical significance
% at the p-value cutoff of 0.05. Note: You may get a different number of
% genes due to the permutation test outcome.
cutoff = 0.05;
sum(pvaluesCorr < cutoff)

%%
% Estimate the FDR and q-values for each test using |mafdr|. The quantity
% _pi0_ is the overall proportion of true null hypotheses in the study. It
% is estimated from the simulated null distribution via bootstrap or the
% cubic polynomial fit. Note: You can also manually set the value of
% lambda for estimating _pi0_.
figure;
[pFDR, qvalues] = mafdr(pvaluesCorr, 'showplot', true);

%%
% Determine the number of genes that have q-values less than the cutoff
% value. Note: You may get a different number of genes due to the
% permutation test and the bootstrap outcomes.
sum(qvalues < cutoff)

%%
% Many genes with low FDR implies that the two groups, MD and Mglio, are
% biologically distinct.

%%
% You can also empirically estimate the FDR adjusted _p-values_ using the
% Benjamini-Hochberg (BH) procedure (Benjamini et al, 1995) by setting the
% |mafdr| input parameter |BHFDR| to true. 
pvaluesBH = mafdr(pvaluesCorr, 'BHFDR', true);
sum(pvaluesBH < cutoff)

%%
% You can store the t-scores, _p-values_, pFDRs, q-values and BH FDR
% corrected _p-values_ together as a DataMatrix object. 
testResults = [tscores pvaluesCorr pFDR qvalues pvaluesBH];
%%
% Update the column name for BH FDR corrected _p-values_ using the
% |colnames| method of DataMatrix object.
testResults = colnames(testResults, 5, {'FDR_BH'});
%%
% You can sort by _p-values_ |pvaluesCorr| using the |sortrows| mathod.
testResults = sortrows(testResults, 2);

%%
% Display the first 23 genes in |testResults|. Note: Your results may be
% different from those shown below due to the permutation test and the
% bootstrap outcomes.
testResults(1:23, :)

%% 
% A gene is considered to be differentially expressed between the two
% groups of samples if it shows both statistical and biological
% significance. In this example, you will compare the gene expression ratio
% of MD over Mglio tumor samples. Therefore an up-regulated gene in this
% example has higher expression in MD and down-regulate gene has higher
% expression in Mglio.

%%
% Plot the -log10 of _p-values_ against the biological effect in a volcano
% plot. Note: From the volcano plot UI, you can interactively change the
% p-value cutoff and fold change limit, and export differentially expressed
% genes. 
diffStruct = mavolcanoplot(MDData, MglioData, pvaluesCorr)

%%
% Ctrl-click genes in the gene lists to label the genes in the plot. As
% seen in the volcano plot, genes specific for neuronal based cerebella
% granule cells, such as _ZIC_ and _NEUROD_, are found in the up-regulated
% gene list, while genes typical of the astrocytic and oligodendrocytic
% lineage and cell differentiation, such as _SOX2_, _PEA15_, and _ID2B_,
% are found in the down-regulated list.  

%%
% Determine the number of differentially expressed genes.
nDiffGenes = diffStruct.PValues.NRows

%%
% Get the list of up-regulated genes for MD compared to Mglio. 
up_geneidx = find(diffStruct.FoldChanges > 0);
up_genes = rownames(diffStruct.FoldChanges, up_geneidx); 
nUpGenes = length(up_geneidx)

%%
% Determine the number of down-regulated genes for MD compared to Mglio.
nDownGenes = sum(diffStruct.FoldChanges < 0)

%% Annotating Up-Regulated Genes Using Gene Ontology
% Use Gene Ontology (GO) to annotate the differentially expressed genes.
% You can look at the up-regulated genes from the analysis above. Download
% the _Homo sapiens_ annotations (|gene_association.goa_human.gz|  file)
% from <http://www.geneontology.org/GO.current.annotations.shtml Gene
% Ontology Current Annotations>, unzip, and store it in your the current
% directory.

%%
% Find the indices of the up-regulated genes for Gene Ontology analysis.
huGenes = rownames(expr_cns_gcrma_eb);
for i = 1:nUpGenes
    up_geneidx(i) = find(strncmpi(huGenes, up_genes{i}, length(up_genes{i})), 1);
end

%%
% Load the Gene Ontology database into a MATLAB object using the |geneont|
% function.
GO = geneont('live',true);

%%
% Read the _Homo sapiens_ gene annotation file. For this example, you will
% look only at genes that are related to molecular function, so you only
% need to read the information where the |Aspect| field is set to 'F'. The
% fields that are of interest are the gene symbol and associated ID. In GO
% Annotation files these have field names |DB_Object_Symbol| and |GOid|
% respectively.
HGann = goannotread('gene_association.goa_human',...
    'Aspect','F','Fields',{'DB_Object_Symbol','GOid'});

%%
% Create a map between annotated genes and GO terms.
HGmap = containers.Map();
for i=1:numel(HGann)
    key = HGann(i).DB_Object_Symbol;
    if isKey(HGmap,key)
        HGmap(key) = [HGmap(key) HGann(i).GOid];
    else
        HGmap(key) = HGann(i).GOid;
    end
end

%%
% Determine the number of _Homo sapiens_ annotated genes related to
% molecular function.
HGmap.Count

%%
% Not all of the 5758 genes on the HuGeneFL chip are annotated. For every
% gene on the chip, see if it is annotated by comparing its gene symbol to
% the list of gene symbols from GO. Track the number of annotated genes and
% the number of up-regulated genes associated with each GO term. Note: You
% might get warnings about invalid or obsolete IDs due to the frequent
% update to the _Homo sapiens_ gene annotation file.
m = GO.Terms(end).id;        % gets the last term id
chipgenesCount = zeros(m,1); % a vector of GO term counts for the entire chip.
upgenesCount  = zeros(m,1);  % a vector of GO term counts for up-regulated genes.
for i = 1:length(huGenes)
    if isKey(HGmap,huGenes{i})
        goid = getrelatives(GO,HGmap(huGenes{i}));
        % Update the tally
        chipgenesCount(goid) = chipgenesCount(goid) + 1;
        if (any(i == up_geneidx))
            upgenesCount(goid) = upgenesCount(goid) +1;
        end
    end
end

%% 
% Determine the statistically significant GO terms using the hypergeometric
% probability distribution. For each GO term, a p-value is calculated
% representing the probability that the number of annotated genes
% associated with it could have been found by chance.  
gopvalues = hygepdf(upgenesCount,max(chipgenesCount),...
                        max(upgenesCount),chipgenesCount);
[dummy, idx] = sort(gopvalues);                   

report = sprintf('GO Term     p-value     counts      definition\n');
for i = 1:10
    term = idx(i);
    report = sprintf('%s%s\t%-1.5f\t%3d / %3d\t%s...\n',... 
             report, char(num2goid(term)), gopvalues(term),...
             upgenesCount(term), chipgenesCount(term),...
            GO(term).Term.definition(2:min(50,end)));
end
disp(report);

%%
% Inspect the significant GO terms and select the terms related to specific
% molecule functions to build a sub-ontology that includes the ancestors of
% the terms. Visualize this ontology using the |biograph| function. You can
% also color the graphs nodes. In this example, the red nodes are the most
% significant, while the blue nodes are the least significant gene ontology
% terms. Note: The GO terms returned may differ from those shown due to the
% frequent update to the _Homo sapiens_ gene annotation file.
fcnAncestors = GO(getancestors(GO,idx(1:5)))
[cm acc rels] = getmatrix(fcnAncestors);
BG = biograph(cm,get(fcnAncestors.Terms,'name'))

for i=1:numel(acc) 
    pval = gopvalues(acc(i));
    color = [(1-pval).^(1) pval.^(1/8) pval.^(1/8)];  
    set(BG.Nodes(i),'Color',color); 
end
view(BG)

%% Finding the Differentially Expressed Genes in Pathways
% You can query the pathway information of the differentially expressed
% genes from the KEGG pathway database through <http://www.genome.jp KEGG's
% Web Service>.
%% 
% Following are a few pathway maps with the genes in the up-regulated gene
% list highlighted:
%%
% <http://www.genome.ad.jp/dbget-bin/show_pathway?hsa04110+3065+5111+4171
% Cell Cycle>
%%
% <http://www.genome.ad.jp/dbget-bin/show_pathway?hsa04340+122011 Hedgehog
% Signaling pathway>
%%
% <http://www.genome.ad.jp/dbget-bin/show_pathway?hsa04150+1975+6194 mTor
% Signaling pathway>

%% References
% [1] Pomeroy, S.L., Tamayo, P., Gaasenbeek, M., Sturla, L.M., Angelo, M.,
%     McLaughlin, M.E., Kim, J.Y., Goumnerova, L.C., Black, P.M., Lau, C.,
%     Allen, J.C., Zagzag, D., Olson, J.M., Curran, T., Wetmore, C.,
%     Biegel, J.A., Poggio, T., Mukherjee, S., Rifkin, R., Califano, A.,
%     Stolovitzky, G., Louis, DN, Mesirov, J.P., Lander, E.S., and Golub,
%     T.R. (2002). Prediction of central nervous system embryonal tumour
%     outcome based on gene expression. Nature, 415(6870), 436-442.
%%
% [2] Storey, J.D., and Tibshirani, R. (2003). Statistical significance for
%     genomewide studies. Proc.Nat.Acad.Sci., 100(16), 9440-9445.
%%
% [3] Dudoit, S., Shaffer, J.P., and Boldrick, J.C. (2003). Multiple
%     hypothesis testing in microarray experiment. Statistical Science, 18,
%     71-103.
%%
% [4] Benjamini, Y., and Hochberg, Y. (1995). Controlling the false
%     discovery rate: a practical and powerful approach to multiple
%     testing. J. Royal Stat. Soc., B 57, 289-300.

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20CNSGENEEXPDEMO%20in%20Bioinformatics%20Toolbox%204.2 Suggest an enhancement for this example.>*

displayEndOfDemoMessage(mfilename)
