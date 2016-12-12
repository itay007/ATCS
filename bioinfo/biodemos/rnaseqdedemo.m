%% Identifying Differentially Expressed Genes from RNA-Seq Data
% This example shows how to load RNA-seq data and test for differential
% expression using a negative binomial model.
 
%   Copyright 2010-2012 The MathWorks, Inc.

%% Introduction
% RNA-seq is an emerging technology for surveying gene expression and
% transcriptome content by directly sequencing the mRNA molecules in a
% sample. RNA-seq can provide gene expression measurements and is regarded
% as an attractive approach to analyze a transcriptome in an unbiased and
% comprehensive manner.
% 
% In this example, you will use Bioinformatics Toolbox(TM) and Statistics
% Toolbox(TM) functions to load publicly available transcriptional profiling
% sequencing data into MATLAB(R), compute the digital gene expression, and
% then identify differentially expressed genes in RNA-seq data from hormone
% treated prostate cancer cell line samples [1].

%% The Prostate Cancer Data Set
% In the prostate cancer study, the prostate cancer cell line LNCap was
% treated with androgen/DHT. Mock-treated and androgen-stimulated LNCap
% cells were sequenced using the Illumina(R) 1G Genome Analyzer [1]. For
% the mock-treated cells, there were four lanes totaling ~10 million reads.
% For the DHT-treated cells, there were three lanes totaling ~7 million
% reads. All replicates were technical replicates. Samples labeled s1
% through s4 are from mock-treated cells. Samples labeled s5, s6, and s8
% are from DHT-treated cells. The read sequences are stored in FASTA files.
% The sequence IDs break down as follows: seq_(unique sequence id)_(number
% of times this sequence was seen in this lane).
%
% This example assumes that you:
% 
% (1) Downloaded and uncompressed the seven FASTA files (|s1.fa|,
%     |s2.fa|, |s3.fa|, |s4.fa|, |s5.fa|, |s6.fa| and |s8.fa|) containing
%     the raw, 35bp, unmapped short reads from
%     <http://yeolab.ucsd.edu/yeolab/Papers.html the author's Web Site>.
% 
% (2) Produced a SAM-formatted file for each of the seven FASTA files by
%     mapping the short reads to the NCBI version 37 of the human genome
%     using a mapper such as Bowtie [2], 
% 
% (3) Ordered the SAM-formatted files by reference name first, then by
%     genomic position.
% 
% For the published version of this example, 4,388,997 short reads
% were mapped using the Bowtie aligner [2]. The aligner was instructed to
% report one best valid alignment. No more than two mismatches were allowed
% for alignment. Reads with more than one reportable alignment were
% suppressed, i.e. any read that mapped to multiple locations was
% discarded. The alignment was output to seven SAM files (|s1.sam|,
% |s2.sam|, |s3.sam|, |s4.sam|, |s5.sam|, |s6.sam| and |s8.sam|). Because
% the input files were FASTA files, all quality values were assumed to be
% 40 on the Phred quality scale [2]. We then used SAMtools [3] to sort the
% mapped reads in the seven SAM files, one for each replicate.

%% Creating an Annotation Object of Target Genes
% Download from Ensembl a tab-separated-value (TSV) table with all
% protein encoding genes to a text file, |ensemblmart_genes_hum37.txt|. For
% this example, we are using Ensamble release 64. Using Ensembl's
% <http://www.ensembl.org/biomart/martview/ BioMart> service, you can
% select a table with the following attributes: chromosome name, gene
% biotype, gene name, gene start/end, and strand direction.
%
% Use the provided helper function |ensemblmart2gff| to convert the
% downloaded TSV file to a GFF formatted file. Then use |GFFAnnotation| to
% load the file into MATLAB.

GFFfilename = ensemblmart2gff('ensemblmart_genes_hum37.txt');
genes = GFFAnnotation(GFFfilename)

%%
% Create a subset with the genes present in chromosomes only (without
% contigs). The |GFFAnnotation| object contais 20012 annotated
% protein-coding genes in the Ensembl database.

chrs = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14',...
        '15','16','17','18','19','20','21','22','X','Y','MT'};  
genes = getSubset(genes,'reference',chrs)   

%% 
% Copy the gene information into a structure and display the first entry.
getData(genes,1)

%% Importing Mapped Short Read Alignment Data
% The size of the sorted SAM files in this data set are in the order of
% 250-360MB. You can access the mapped reads in |s1.sam| by creating a
% |BioMap|. |BioMap| has an interface that provides direct access to the
% mapped short reads stored in the SAM-formatted file, thus minimizing the
% amount of data that is actually loaded into memory.

bm = BioMap('s1.sam')

%%
% Use the |getSummary| method to obtain a list of the existing references
% and the actual number of short read mapped to each one. Observe that the
% order of the references is equivalent to the previously created cell
% string |chrs|.

getSummary(bm)

%%
% You can access the alignments, and perform operations like getting counts
% and coverage from |bm|. For more examples of getting read coverage at the
% chromosome level, see <chipseqpedemo.html Exploring Protein-DNA Binding
% Sites from Paired-End ChIP-Seq Data>. 

%% Determining Digital Gene Expression
% Next, you will determine the mapped reads associated with each Ensembl
% gene. Because the strings used in the SAM files to denote the reference
% names are different to those provided in the annotations, we find a
% vector with the reference index for each gene: 

geneReference =  seqmatch(genes.Reference,chrs,'exact',true);

%%
% For each gene, count the mapped reads that overlap any part of the
% gene. The read counts for each gene are the digital gene expression of
% that gene. Use the |getCounts| method of a |BioMap| to compute the read
% count within a specified range.

counts = getCounts(bm,genes.Start,genes.Stop,1:genes.NumEntries,geneReference);

%%
% Gene expression levels can be best respresented by a |DataMatrix|, with
% each row representing a gene and each column representing a sample.
% Create a |DataMatrix| with seven columns, one for each sample. Copy the
% counts of the first sample to the first column.

filenames = {'s1.sam','s2.sam','s3.sam','s4.sam','s5.sam','s6.sam','s8.sam'};
colnames =  {'Mock_1','Mock_2','Mock_3','Mock_4','DHT_1','DHT_2','DHT_3'};

lncap_counts = bioma.data.DataMatrix(NaN([genes.NumEntries,7]),genes.Feature,colnames);
lncap_counts(:,1) = counts;

lncap_counts(1:10,:)

%%
% Determine the number of genes that have counts greater than or equal to
% 50 in chromosome 1.
lichr1 = geneReference == 1;  % logical index to genes in chromosome 1
sum(lncap_counts(:,1) >= 50 & lichr1)

%%
% Repeat this step for the other six samples (SAM files) in the data set to
% get their gene counts and copy the information to the previously
% created |DataMatrix|.

for i = 2:7 
    bm = BioMap(filenames{i});
    counts = getCounts(bm,genes.Start,genes.Stop,1:genes.NumEntries,geneReference);
    lncap_counts(:,i) = counts;
end

%% 
% Inspect the first 10 rows in the count table.
lncap_counts(1:10, :)

%%
% The DataMatrix |lncap_counts| contains counts for samples from two
% biological conditions: mock-treated and DHT-treated.
cond_Mock = logical([1 1 1 1 0 0 0]);
cond_DHT  = logical([0 0 0 0 1 1 1]);

%%
% You can plot the counts for a chromosome along the chromosome genome
% coordinate. For example, plot the counts for chromosome 1 for
% mock-treated sample |Mock_1| and DHT-treated sample |DHT_1|. Add the
% ideogram for chromosome 1 to the plot using the |chromosomeplot|
% function.

ichr1 = find(lichr1);  % linear index to genes in chromosome 1
[~,h] = sort(genes.Start(ichr1));
ichr1 = ichr1(h);      % linear index to genes in chromosome 1 sorted by 
                       % genomic position

figure
plot(genes.Start(ichr1), lncap_counts(ichr1,'Mock_1'), '.-r',...
     genes.Start(ichr1), lncap_counts(ichr1,'DHT_1'), '.-b');
ylabel('Gene Counts')
title('Gene Counts on Chromosome 1')
fixGenomicPositionLabels(gca)  % formats tick labels and adds datacursors
chromosomeplot('hs_cytoBand.txt', 1, 'AddToPlot', gca)

%% Inference of Differential Signal in RNA Expression
% For RNA-seq experiments, the read counts have been found to be linearly
% related to the abundance of the target transcripts [4]. The interest lies
% in comparing the read counts between different biological conditions.
% Current observations suggest that typical RNA-seq experiments have low
% background noise, and the gene counts are discrete and could follow the
% Poisson distribution. While it has been noted that the assumption of the
% Poisson distribution often predicts smaller variation in count data by
% ignoring the extra variation due to the actual differences between
% replicate samples [5]. Anders _et.al._,(2010) proposed an error model for
% statistical inference of differential signal in RNA-seq expression data
% that could address the overdispersion problem [6]. Their approach uses
% the negative binomial distribution to model the null distribution of the
% read counts. The mean and variance of the negative binomial distribution
% are linked by local regression, and these two parameters can be reliably
% estimated even when the number of replicates is small [6].

%%
% In this example, you will apply this statistical model to process the
% count data and test for differential expression. The details of the
% algorithm can be found in reference [6]. The model of Anders _et.al._,
% (2010) has three sets of parameters that need to be estimated from the
% data:
% 
% 1. Library size parameters; 
% 
% 2. Gene abundance parameters under each experimental condition;
% 
% 3. The smooth functions that model the dependence of the raw
%          variance on the expected mean. 

%% Estimating Library Size Factor
% The expectation values of all gene counts from a sample are proportional
% to the sample's library size. The effective library size can be estimated
% from the count data.
%%
% Compute the geometric mean of the gene counts (rows in |lncap_counts|)
% across all samples in the experiment as a pseudo-reference sample.
geoMeans = exp(mean(log(lncap_counts), 2));

%%
% Each library size parameter is computed as the median of the ratio of
% the sample's counts to those of the pseudo-reference sample.
ratios = dmbsxfun(@rdivide, lncap_counts(geoMeans >0, :), geoMeans(geoMeans >0));
sizeFactors = median(ratios, 1);
                               
%% 
% The counts can be transformed to a common scale using size factor
% adjustment.
base_counts = dmbsxfun(@rdivide, lncap_counts, sizeFactors);

%%
% Use the |boxplot| function to inspect the count distribution of the
% mock-treated and DHT-treated samples and the size factor adjustment.
figure
subplot(2,1,1)
maboxplot(log2(lncap_counts), 'title','Raw Read Counts',...
                              'orientation', 'horizontal') 
subplot(2,1,2)
maboxplot(log2(base_counts), 'title','Size Factor Adjusted Read Counts',...
                             'orientation', 'horizontal') 
      
%% Estimating Negative Binomial Distribution Parameters
% The expectation value of counts for a gene is also proportional to the
% gene abundance parameter. You can estimate the gene abundance parameter
% from the average of counts from samples corresponding to an experimental
% condition. For example, compute the mean counts and sample variances from
% mock-treated samples.
base_mean_mock = mean(base_counts(:, cond_Mock), 2);
base_var_mock = var(base_counts(:, cond_Mock), 0, 2);

%%
% To avoid code duplication in the example for computing parameters for
% samples of different conditions, we provide a helper function,
% |estimateBaseParams|, to compute the mean, the variance, the smooth
% function fit data for raw variance estimation, and the diagnostic
% variance residual distribution from replicates under the same condition.
% For example, compute the base means and variances for DHT-treated
% samples.
[base_mean_dht, base_var_dht] = estimateBaseParams(lncap_counts(:, cond_DHT),...
                                                   sizeFactors(cond_DHT),...
                                                   'MeanAndVar');
%%
% In the model, the full variances of the negative binomial distribution of
% the counts of a gene are considered as the sum of a shot noise term and a
% raw variance term. The shot noise term is the read counts of the gene,
% while the raw variance can be predicted from the mean, i.e., genes with a
% similar expression level have similar variance across the replicates
% (samples of the same biological condition). A smooth function that models
% the dependence of the raw variance on the mean is obtained by fitting the
% sample mean and variance within replicates for each gene using the local
% regression function |malowess|. For example, get the smooth fit data from
% the sample mean and variance of the mock-treated samples.
[rawVarSmooth_X_mock, rawVarSmooth_Y_mock] = ...
                                estimateBaseParams(lncap_counts(:, cond_Mock),...
                                                   sizeFactors(cond_Mock),...
                                                   'SmoothFunc');

%%
% Find the raw variances for each gene from its base mean value by
% interpolation.
raw_var_mock_fit = interp1(rawVarSmooth_X_mock, rawVarSmooth_Y_mock,...
                           log(base_mean_mock), 'linear', 0);
                   
%%
% Add the bias correction term [6] to get the raw variances.
zConst = sum(1 ./sizeFactors(cond_Mock), 2) / length(sizeFactors(cond_Mock));
raw_var_mock = raw_var_mock_fit - base_mean_mock * zConst;

%%
% Plot the sample variance and the raw variance data to check the fit of
% the variance function.
[base_mean_mock_sort, sidx] = sort(log10(base_mean_mock));
raw_var_mock_sort = log10(raw_var_mock_fit(sidx));

figure
plot(log10(base_mean_mock), log10(base_var_mock), '*')
hold on
line(base_mean_mock_sort, real(raw_var_mock_sort), 'Color', 'r', 'LineWidth',2)
ylabel('log10(base variances) of mock-treated samples')
xlabel('log10(base means) of mock-treated samples')

%%
% The fit (red line) follows the single-gene estimates well, even though
% the spread of the latter is considerable, as one would expect, given that
% each raw variance value is estimated from only four values (four
% mock-treaded replicates).

%%
% As RNA-seq experiments typically have few replicates, the single-gene
% estimate of the base variance can deviate wildly from the fitted value.
% To see whether this might be too wild, the cumulative probability for the
% ratio of single-gene estimate of the base variance to the fitted value is
% calculated from the chi-square distribution, as explained in reference
% [6].

%%
% Compute the cumulative probabilities of the variance ratios of
% mock-treated samples.
df_mock = sum(cond_Mock) - 1;
varRatio_mock = base_var_mock ./ raw_var_mock_fit;
pchisq_mock = chi2cdf(df_mock * varRatio_mock, df_mock);

%%
% Compute the empirical cumulative density functions (ECDF) stratified by
% base count levels, and show the ECDFs curves. Group the counts into seven
% levels.
count_levels = [0 3; 3.1 12; 12.1 30; 30.1 65; 65.1 130; 130.1 310; 310.1 2500];
%%
figure;
hold on
cm = jet(7);
for i = 1:7
   [Y1,X1] = ecdf(pchisq_mock(base_mean_mock>count_levels(i, 1) &...
                              base_mean_mock<count_levels(i,2)));
   plot(X1,Y1,'LineWidth',2,'color',cm(i,:))
end
plot([0,1],[0,1] ,'k', 'linewidth', 2)
set(gca, 'Box', 'on')
legend('0-3', '3-12', '12-30', '31-65', '65-130', '131-310', '311-2500',...
       'Location','NorthWest')
xlabel('Chi-squared probability of residual')
ylabel('ECDF')
title('Residuals ECDF plot for mock-treated samples')

%%
% The ECDF curves of count levels greater than 3 and below 130 follows the
% diagonal well (black line). If the ECDF curves are below the black line,
% variance is underestimated. If the ECDF curves are above the black line,
% variance is overestimated [6]. For very low counts (below 3), the
% deviations become stronger, but at these levels, shot noise dominates.
% For the high count cases, the variance is overestimated. The reason might
% be there aren't enough genes with high counts. Get the number of genes in
% each of the count levels.
num_in_count_levels = zeros(1, 7);
for i = 1:7
    num_in_count_levels(i) = sum(base_mean_mock>count_levels(i, 1) & ...
                                 base_mean_mock<count_levels(i,2));
end
%%
num_in_count_levels
%%
% Increasing the sequence depth, which in turn increases the number of
% genes with higher counts, improves the variance estimation.
%%
% You can produce the same ECDF plot for the DHT-treated samples by
% following the same steps.
pchisq_dht = estimateBaseParams(lncap_counts(:, cond_DHT),...
                                                sizeFactors(1, cond_DHT),...
                                                'Diagnostic');
%%
figure;
hold on
for i = 1:7
   [Y1,X1] = ecdf(pchisq_dht(base_mean_dht>count_levels(i, 1) & ...
                             base_mean_dht<count_levels(i,2)));
   plot(X1,Y1,'LineWidth',2,'color',cm(i,:))
end
plot([0,1],[0,1] ,'k', 'linewidth', 2)
set(gca, 'Box', 'on')
legend('0-3', '3-12', '12-30', '31-65', '65-130', '131-310', '311-2500',...
       'Location','NorthWest')
xlabel('Chi-squared probability of residual')
ylabel('ECDF')
title('Residuals ECDF plot for DHT-treated samples')

%%
% In both cases, most of the ECDF curves follow the diagonal well. The fits
% are reasonably good.

%% Testing for Differential Expression
% Having estimated and verified the mean-variance dependence, you can test
% for differentially expressed genes between the samples from the mock- and
% DHT- treated conditions. For your convenience, we provide the helper
% function |estimateNBParams| to estimate the mean and full variance of the
% two-parametric negative binomial distribution for each gene from the
% three sets of parameters discussed above.
%%
[mu_mock, full_var_mock, mu_dht, full_var_dht] =...
          estimateNBParams(lncap_counts, sizeFactors, cond_DHT, cond_Mock);
%%
% Compute the p-values for the statistical significance of the change from
% DHT-treated condition to mock-treated condition. The helper function
% |computePVal| implements the numerical computation of the p-values
% presented in the reference [6]. We use the |nbinpdf| function to compute
% the negative binomial probability density. Note: The computation was not
% optimized, and it will take several minutes to run.
%%
% Get the gene counts for each condition:
k_mock = sum(lncap_counts(:, cond_Mock), 2);
k_dht = sum(lncap_counts(:, cond_DHT), 2); 
%%
pvals =  computePVal(k_dht, mu_dht, full_var_dht,...
                          k_mock, mu_mock, full_var_mock);
%%
% You can empirically adjust the p-values from the multiple tests for false
% discovery rate (FDR) with the Benjamini-Hochberg procedure [7] using the
% |mafdr| function.
p_fdr = mafdr(pvals, 'BHFDR', true);

%%
% Determine the fold change estimated from the DHT-treated to the
% mock-treated condition.
foldChange = base_mean_dht ./ base_mean_mock;
%% 
% Determine the base 2 logarithm of the fold change.
log2FoldChange = log2(foldChange);

%%
% Determine the mean expression level estimated from both conditions.
base_mean_com = estimateBaseParams(lncap_counts, sizeFactors, 'MeanAndVar');

%%
% Assume a p-value cutoff of 0.01.
de_idx = p_fdr < 0.01; 

%%
% Plot the log2 fold changes against the base means, and color those genes
% with p-values less than the cutoff value red.
figure;
plot(log2(base_mean_com(~de_idx, :)), log2FoldChange(~de_idx,:), 'b.')
hold on
plot(log2(base_mean_com(de_idx, :)), log2FoldChange(de_idx, :), 'r.')
xlabel('log2 Mean')
ylabel('log2 Fold Change')

%%
% You can identify up- or down- regulated genes for mean base count levels
% over 3.
up_idx = find(p_fdr < 0.01 & log2FoldChange >= 2 & base_mean_com > 3 );
numel(up_idx)
%%
down_idx = find(p_fdr < 0.01 & log2FoldChange <= -2 & base_mean_com > 3 );
numel(down_idx)

%%
% This analysis identified 462 genes (out of 20,012 genes) that were
% differentially up- or down- regulated by hormone treatment.

%% References
% [1] Li, H., Lovci, M.T., Kwon, Y-S., Rosenfeld, M.G., Fu, X-D., and Yeo,
%     G.W. "Determination of Tag Density Required for Digital Transcriptome
%     Analysis: Application to an Androgen-Sensitive Prostate Cancer
%     Model", PNAS, 105(51), pp 20179-20184, 2008.
%%
% [2] Langmead, B., Trapnell, C., Pop, M., and Salzberg, S.L. "Ultrafast
%     and Memory-efficient Alignment of Short DNA Sequences to the Human
%     Genome", Genome Biology, 10:R25, pp 1-10, 2009.
%%
% [3] Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N.,
%     Marth, G., Abecasis, G., Durbin, R. and 1000 Genome Project Data
%     Processing Subgroup, "The Sequence Alignment/map (SAM) Format and
%     SAMtools", Bioinformatics, 25, pp 2078-2079, 2009.
%%
% [4] Mortazavi, A., Williams, B.A., McCue, K., Schaeffer, L., and Wold, B.
%     "Mapping and quantifying mammalian transcriptomes by RNA-Seq", Nature
%     Methods, 5, pp 621-628, 2008.
%%
% [5] Robinson, M.D., and Oshlack, A. "A Scaling Normalization method for
%     differential Expression Analysis of RNA-seq Data", Genome Biology
%     11:R25, 1-9, 2010.
%%
% [6] Anders, S. and Huber W. "Differential Expression Analysis for
%     Sequence Count Data", Genome Biology, 11:R106, 2010.
%%
% [7] Benjamini, Y., and Hochberg, Y.  "Controlling the false discovery
%     rate: a practical and powerful approach to multiple testing", J.
%     Royal Stat. Soc., B 57, 289-300, 1995.
%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20RNASEQDEDEMO%20in%20Bioinformatics%20Toolbox%204.2 Suggest an enhancement for this example.>*

displayEndOfDemoMessage(mfilename)
