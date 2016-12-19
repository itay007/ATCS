%% Analyzing Illumina(R) Bead Summary Gene Expression Data
% This example shows how to analyze Illumina BeadChip gene expression
% summary data using MATLAB(R) and Bioinformatics Toolbox(TM) functions.
 
%   Copyright 2009-2012 The MathWorks, Inc.


%% Introduction
% This example shows a number of ways to import and analyze gene expression
% data from the Illumina BeadChip microarray platform. The data set in the
% example is from the study of gene expression profiles of human
% spermatogenesis by Platts et al., 2007. The expression levels were
% measured on Illumina Sentrix Human 6 (or WG6) BeadChips.
% 
% Each Illumina WG6 BeadChip contains six identical arrays with 47,293
% unique probes. Illumina's BeadStudio(TM) software outputs the summarized
% expression levels for each bead type on the arrays on the BeadChip. The
% output from the BeadStudio software can be either a text file or an
% Excel(R) file. 
% 
% Both raw and normalized Illumina expression data are available on the
% Gene Expression Omnibus (GEO) database Web site:
% http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6967.   
% 
% The data from most microarray gene expression experiments usually
% consists of four components: experiment data values, sample information,
% feature annotations, and information about the experiment. You will work
% with the data from the experiment, construct each of the four components,
% assemble them into an |ExpressionSet| object, and find the differentially
% expressed genes. For more examples about the |ExpressionSet| class, see
% <maexptdemo.html Working with Objects for Microarray Experiment Data>.

%% Importing Experiment Data from the GEO Database
% Samples were hybridized on three Illumina Sentrix Human 6 (or WG6)
% BeadChips in the experiment. Use the |getgeodata| function to retrieve
% the GEO Series data, _GSE6967_, and read it into the MATLAB Workspace as
% a structure, |TNGEOData|. You can also download the raw Illumina summary
% data files of the _GSE6967_ record via FTP from the GEO database. 

%%
TNGEOData = getgeodata('GSE6967')

%%
% The |TNGEOData| structure contains |Header| and |Data| fields. The
% |Header| field has two fields, |Series| and |Samples|, containing a
% description of the experiment and samples respectively. The |Data| field
% contains a |DataMatrix| of normalized summary expression levels from the
% experiment.

%%
% Determine the number of samples in the experiment.
nSamples = numel(TNGEOData.Header.Samples.geo_accession)

%%
% Inspect the sample titles from the |Header.Samples| field.
TNGEOData.Header.Samples.title'

%%
% For simplicity, extract the sample labels from the sample titles.
sampleLabels = cellfun(@(x) char(regexp(x, '\w\d+', 'match')),...
                TNGEOData.Header.Samples.title, 'UniformOutput',false)

%% Importing Expression Data from Illumina BeadStudio Summary Files 
% Raw, non-normalized summary data was deposited as supplementary data for
% the _GSE6967_ record in the GEO database. Download the supplementary file
% |GSE6967_RAW.tar|. Unzip the file to access the 13 text files produced by
% the BeadStudio software, which contain the raw, non-normalized bead
% summary values. 
% 
% There is one file for each array (or each sample). Inside each file is a
% line of column headers followed by a data matrix with 47,293 rows. Each
% row corresponds to a different target probe (gene) in the experiment. The
% matrix contains the summarized expression values (Avg_Signal), standard
% error of the bead replicates (BEADSTDEV), number of beads used
% (Avg_NBEADS) and a detection score, which estimates the confidence limit
% of detection of a target probe from the samples hybridized on the
% BeadChip. The column header information also includes the Sentrix chip
% IDs and sample IDs placed on each array.

%%
% The raw data text files are named with their GSM accession numbers. For
% this example, construct the file names of the text data files using the
% path where the text files are located. 
rawDataFiles = cell(1,nSamples);
for i = 1:nSamples
    rawDataFiles {i} = [TNGEOData.Header.Samples.geo_accession{i} '.txt'];
end

%%
% Modify this line to contain the path and directory to which you extracted
% the data files. 
rawDataPath = fullfile('C:', 'Examples', 'illuminagedemo', 'GSE6967')

%%
% Use the |ilmnbsread| function to read one of the summary files and
% inspect the returned structure.
rawData =ilmnbsread(fullfile(rawDataPath, rawDataFiles{1}))

%%
% Inspect the column names in the |rawData| structure.
rawData.ColumnNames'

%%
% Determine the number of target probes.
nTargets = size(rawData.Data, 1)
%%
% Read the non-normalized expression values (Avg_Signal), the detection
% confidence limits and the Sentrix chip IDs from the 13 summary data
% files. The gene expression values are identified with Illumina probe
% target IDs. 
rawMatrix = bioma.data.DataMatrix(zeros(nTargets, nSamples),...
                                  rawData.TargetID, sampleLabels);
detectionConf = bioma.data.DataMatrix(zeros(nTargets, nSamples),...
                                      rawData.TargetID, sampleLabels);
chipIDs = cell(1, nSamples);
%%
% You can specify the columns to read from the data file.
for i = 1:nSamples
    rawData =ilmnbsread(fullfile(rawDataPath, rawDataFiles{i}),...
                                'COLUMNS', {'AVG_Signal', 'Detection'});
    chipIDs(i) = regexp(rawData.ColumnNames(1), '\d*', 'match', 'once');
    rawMatrix(:, i) = rawData.Data(:, 1);
    detectionConf(:,i) = rawData.Data(:,2);
end

%%
% There are three Sentrix BeadChips used in the experiment. Inspect the
% Illumina Sentrix BeadChip IDs in |chipIDs| and determine the number of
% samples hybridized on each chip. 
[beadChipIDs, I] = unique(chipIDs);
beadChipIDs
diff([0 I])

%%
sampleLabels(1:I(1))
sampleLabels(I(1)+1:I(2))
sampleLabels(I(2)+1:I(3))
%%
% Six samples (T2, T1, T6, T4, T8 and N11) were hybridized to six arrays on
% the first chip, four samples (T3, T7, T5 and N6) on the second chip, and
% three samples (N12, N5, and N1) on the third chip.

%% Normalizing the Expression Data
% Use a boxplot to view the raw expression levels of each sample in the
% experiment.
logRawExprs = log2(rawMatrix);

%%
maboxplot(logRawExprs,'ORIENTATION', 'horizontal')
ylabel('Sample Labels')
xlabel('log2(Expression Value)')
title('Before Normalization')
%%
% The difference in intensities between samples on the same chip and
% samples on different chips does not seem too large. The first BeadChip,
% containing samples T2, T1, T6, T4, T8 and N11, seems to be slightly more
% variable than others.
%%
% Using MA and XY plots to do a pairwise comparison of the arrays on a BeadChip
% can be informative. On an MA plot, the average of the expression levels
% of two arrays (A) are plotted on the x axis, and the difference in the
% measurement (M) on the y axis. An XY plot is a scatter plot of one array
% against another. In this example, you will use the helper function
% |maxyplot| to plot MAXY plots for a pairwise comparison of the three
% arrays on the first chip hybridized with teratozoospermic samples (T2, T1
% and T6).
% 
% *Note*: You can also use the |mairplot| function to create the MA or IR
% (Intensity/Ratio) plots for comparison of specific arrays. 
inspectIdx = 1:3;
maxyplot(rawMatrix, inspectIdx)
suptitle('Before Normalization')

%%
% In an MAXY plot, the MA plots for all pairwise comparisons are in the
% upper-right corner. In the lower-left corner are the XY plots of the
% comparisons. The MAXY plot shows the two arrays, T1 and T2, to be quite
% similar, while different from the other array, T6. 

%% 
% The expression boxplots and MAXY plots reveal that, there are differences
% in expression levels within chips and between chips; hence, the data
% requires normalization. Use the |quantilenorm| function to apply quantile
% normalization to the raw data. 
% 
% *Note*: You can also try invariant set normalization using the
% |mainvarsetnorm| function.
normExprs = rawMatrix;
normExprs(:, :) = quantilenorm(rawMatrix.(':')(':'));

%%
log2NormExprs = log2(normExprs);

%%
% Display and inspect the normalized expression levels in a boxplot.
figure;
maboxplot(log2NormExprs,'ORIENTATION', 'horizontal')
ylabel('Sample Labels')
xlabel('log2(Expression Value)')
title('After Quantile Normalization')

%%
% Display and inspect the MAXY plot of the three arrays (T2, T1 and T6) on
% the first chip after the normalization.
maxyplot(normExprs, inspectIdx)
suptitle('After Quantile Normalization')

%%
% Many of the genes in this study are not expressed, or have only small
% variability across the samples. Remove these genes using non-specific
% filtering.

%%
% Use the |genelowvalfilter| function to filter out genes with very low
% absolute expression values.
[mask, log2NormExprs] = genelowvalfilter(log2NormExprs);
detectionConf = detectionConf(mask, :);
%%
% Use the |genevarfilter| function to filter out genes with a small
% variance across samples.
[mask, log2NormExprs] = genevarfilter(log2NormExprs);
detectionConf = detectionConf(mask, :);

%% Importing Feature Metadata from a BeadChip Annotation File
% Microarray manufactures usually provide annotations of a collection of
% features for each type of chip. The chip annotation files contain
% metadata such as the gene name, symbol, NCBI accession number, chromosome
% location and pathway information. Before assembling an |ExpressionSet|
% object for the experiment. Obtain the annotations about the features or
% probes on the BeadChip. You can download the |Human_WG-6.csv| annotation
% file for Sentrix Human 6 (or WG6) BeadChips from the Support page at the
% <http://www.illumina.com Illumina> web site and save the file locally.
% Read the annotation file into MATLAB as a |dataset| array. Modify this
% line to contain the path and directory to which you downloaded the
% annotation file. 
annotPath = fullfile('C:', 'Examples', 'illuminagedemo', 'Annotation'); 

%%
WG6Annot = dataset('xlsfile', fullfile(annotPath, 'Human_WG-6.csv'));

%%
% Inspect the properties of this |dataset| array.
get(WG6Annot)

%%
% Get the names of variables describing the features on the Sentrix Human
% 6 BeadChips. 
fDataVariables = get(WG6Annot, 'VarNames')

%%
% Check the number of probe target IDs in the annotation file.
numel(WG6Annot.Target)

%%
% Because the expression data in this example is only a small set of the
% full expression values, you will work with only the features in the
% |DataMatrix| object |log2NormExprs|. Find the matching features in
% |log2NormExprs| and |WG6Annot.Target|.  
[commTargets, fI, WGI] =intersect(rownames(log2NormExprs), WG6Annot.Target);
        
%% Building an ExpressionSet Object - Experiment Data Values 
% You can store the preprocessed expression values and detection limits of
% the annotated probe targets as an |ExptData| object.
fNames = rownames(log2NormExprs);
TNExptData = bioma.data.ExptData({log2NormExprs(fI, :), 'ExprsValues'},...
                                 {detectionConf(fI, :), 'DetectionConfidences'})                             

%% Building an ExpressionSet Object - Sample Information
% The sample data in the |Header.Samples| field of the |TNGEOData|
% structure can be overwhelming and difficult to navigate through. From the
% data in |Header.Samples| field, you can gather the essential information
% about the samples, like the sample titles, GEO sample accession numbers,
% etc., in the experiment, and store the sample data as a |MetaData|
% object. 

%%
% Retrieve the descriptions about sample characteristics.
sampleChars = cellfun(@(x) char(regexp(x, '\w*tile', 'match')),...
               TNGEOData.Header.Samples.characteristics_ch1, 'UniformOutput',false)
           
%%
% Create a |dataset| array to store the sample data you just extracted.
sampleDS = dataset({TNGEOData.Header.Samples.geo_accession', 'GSM'},...
                   {strtok(TNGEOData.Header.Samples.title)', 'Type'},...
                   {sampleChars', 'Characteristics'}, 'obsnames', sampleLabels') 
               
%%
% Store the sample metadata as an object of the |MetaData| class, including
% a short description for each variable.
TNSData = bioma.data.MetaData(sampleDS,...
    {'Sample GEO accession number',...
    'The spermic type of individual whoes sample was collected ',...
    'Sample individual fertility characteristics'})

%% Building an ExpressionSet Object - Feature Annotations
% The collection of feature metadata for Sentrix Human 6 BeadChips is large
% and diverse. Select information about features that are unique to the
% experiment and save the information as a |MetaData| object. Extract
% annotations of interests, for example, |Accession| and |Symbol|.
fIdx = ismember(fDataVariables, {'Accession', 'Symbol'});

%%
featureAnnot = WG6Annot(WGI, fDataVariables(fIdx));

%%
featureAnnot = set(featureAnnot, 'ObsNames', WG6Annot.Target(WGI));

%%
% Create a |MetaData| object for the feature annotation information with
% brief descriptions about the two variables of the metadata.
WG6FData = bioma.data.MetaData(featureAnnot, ...
            {'Accession number of probe target', 'Gene Symbol of probe target'})
        
%% Building an ExpressionSet Object - Experiment Information
% Most of the experiment descriptions in the |Header.Series| field of the
% |TNGEOData| structure can be reorganized and stored as a |MIAME| object,
% which you will use to assemble the |ExpressionSet| object for the
% experiment. 
TNExptInfo = bioma.data.MIAME(TNGEOData)

%% Building an ExpressionSet Object - Assembling the ExpressionSet Object
% Now that you've created all the components, you can create an object of
% the |ExpressionSet| class to store the expression values, sample
% information, chip feature annotations and description information about
% this experiment.
TNExprSet = bioma.ExpressionSet(TNExptData, 'sData', TNSData,...
                                            'fData', WG6FData,...
                                            'eInfo', TNExptInfo)
%% 
% *Note*: The |ExprsValues| element in the |ExptData| object, |TNExptData|,
% is renamed to |Expressions| in |TNGeneExprSet|.

%%
% You can save an object of |ExpressionSet| class as a |MAT| file for
% further data analysis.
save TNGeneExprSet TNExprSet

%%
clear all 

%% Profiling Gene Expression - Grouping Samples
% Load the experiment data saved from the previous section. You will find
% differentially expressed genes between the teratozoospermia and normal
% samples.
load TNGeneExprSet

%%
% Group samples into two variables: |Tz|, consisting of eight
% teratozoospermia samples, and |Ns|, consisting of five normospermic
% reproductive samples. From the expression data of all 13 samples, extract
% the data of the two different groups.
TNSamples = sampleNames(TNExprSet);
Tz = strncmp(TNSamples, 'T', 1);
Ns = strncmp(TNSamples, 'N', 1);
nTz = sum(Tz)
nNs = sum(Ns)

%% Profiling Gene Expression - Permutation T-tests
% To identify the differential changes in the levels of transcripts in
% normospermic |Ns| and teratozoospermic |Tz| samples, compare the gene
% expression values between the two groups of data: |Tz| and |Ns|.
TNExprs = expressions(TNExprSet);
TzData = TNExprs(:,Tz);
NsData = TNExprs(:,Ns);
meanTzData = mean(TzData,2); 
meanNsData = mean(NsData,2);
groupLabels = [TNSamples(Tz), TNSamples(Ns)];

%%
% Perform a permutation t-test using the |mattest| function to permute the
% columns of the gene expression data matrix of |TzData| and |NsData|.
% Note: Depending on the sample size, it may not be feasible to
% consider all possible permutations. Usually, a random subset of
% permutations are considered in the case of a large sample size. 

%%
% Use the |nchoosek| function in Statistics Toolbox(TM) to determine the
% number of all possible permutations of the samples in this example.
perms = nchoosek(1:nTz+nNs, nTz);
nPerms = size(perms,1)

%%
% Use the |PERMUTE| option of the |mattest| function to compute the
% p-values of the 1,287 permutations.
pValues = mattest(TzData, NsData, 'Permute', nPerms);

%%
% You can also compute the differential score from the p-values using this
% anonymous function [1]. 
diffscore = @(p, avgTest, avgRef) -10*sign(avgTest - avgRef).*log10(p);

%%
% The differential score of 13 corresponds to a p-value of 0.05, the
% differential score of 20 corresponds to a p-value of 0.01, and the
% differential score of 30 corresponds to a p-value of 0.001. A positive
% differential score represents up regulation, while a negative score
% represents down regulation.
diffScores = diffscore(pValues, meanTzData, meanNsData);

%%
% Determine the number of genes considered to have a differential score
% greater than 20. Note: You may get a different number of genes due to the
% permutation test outcome.
sum(diffScores > 20)

%%
sum(diffScores < -20)

%% Profiling Gene Expression - Estimating FDR
% In multiple hypothesis testing, which simultaneously tests the null
% hypothesis of thousands of genes from microarray expression data, each
% test has a specific false positive rate, or a false discovery rate (FDR)
% (Storey et al., 2003). Estimate the FDR and q-values for each test using
% the |mafdr| function.
figure;
[pFDR, qValues] = mafdr(pValues, 'showplot', true);
diffScoresFDRQ = diffscore(qValues, meanTzData, meanNsData);

%%
% Determine the number of genes with an absolute differential score greater
% than 20. Note: You may get a different number of genes due to the
% permutation test and the bootstrap outcomes.
sum(abs(diffScoresFDRQ)>=20)

%% Profiling Gene Expression - Differential Gene Expression
% Plot the _-log10_ of p-values against fold changes in a volcano plot. 
diffStruct = mavolcanoplot(TzData, NsData, qValues,...
                                   'pcutoff', 0.01, 'foldchange', 5);

%%
% Note: From the volcano plot UI, you can interactively change the p-value
% cutoff and fold-change limit, and export differentially expressed genes. 
%%
% Determine the number of differentially expressed genes.
nDiffGenes = numel(diffStruct.GeneLabels)

%%
% Get the list of up-regulated genes for the |Tz| samples compared to the
% |Ns| samples. 
up_genes = diffStruct.GeneLabels(diffStruct.FoldChanges > 0);
nUpGenes = length(up_genes)

%%
% Get the list of down-regulated genes for the |Tz| samples compared to the
% |Ns| samples.
down_genes = diffStruct.GeneLabels(diffStruct.FoldChanges < 0);
nDownGenes = length(down_genes)

%%
% Extract a list of differentially expressed genes.
diff_geneidx = zeros(nDiffGenes, 1);
for i = 1:nDiffGenes
    diff_geneidx(i) = find(strncmpi(TNExprSet.featureNames, ...
                            diffStruct.GeneLabels{i}, length(diffStruct.GeneLabels{i})), 1);
end

%%
% You can get the subset of experiment data containing only the
% differentially expressed genes.
TNDiffExprSet = TNExprSet(diff_geneidx, groupLabels);

%% PCA and Clustering Analysis of Significant Gene Profiles
% Principal component analysis (PCA) on differentially expressed genes
% shows linear separability of the |Tz| samples from the |Ns| samples.
PCAScore = pca(TNDiffExprSet.expressions);

%%
% Display the coefficients of the first and sixth principal components.
figure;
plot(PCAScore(:,1), PCAScore(:,6), 's', 'MarkerSize',10, 'MarkerFaceColor','g');
hold on
text(PCAScore(:,1)+0.02, PCAScore(:,6), TNDiffExprSet.sampleNames)
plot([0,0], [-0.5 0.5], '--r')
set(gca, 'xtick', [], 'ytick', [], 'yticklabel',[]);
title('PCA Mapping')
xlabel('Principal Component 1')
ylabel('Principal Component 6')

%%
% You can also use the interactive tool created by the |mapcaplot| function
% in the Bioinformatics Toolbox to perform principal component analysis.
mapcaplot((TNDiffExprSet.expressions)')

%%
% Perform unsupervised hierarchical clustering of the significant gene
% profiles from the |Tz| and |Ns| groups using correlation as the distance
% metric to cluster the samples.
sampleDist = pdist(TNDiffExprSet.expressions','correlation');
sampleLink = linkage(sampleDist);
%%
figure;
dendrogram(sampleLink, 'labels', TNDiffExprSet.sampleNames,'ColorThreshold',0.5)
set(gca, 'ytick', [], 'box', 'on');
title('Hierarchical Sample Clustering')

%%
% Use the |clustergram| function to create the hierarchical clustering of
% differentially expressed genes, and apply the colormap |redbluecmap| to
% the clustergram.
cmap = redbluecmap(9);
%%
cg = clustergram(TNDiffExprSet.expressions,'Colormap',cmap);        
%%    
cg.addTitle('Hierarchical Gene Clustering')
%%
% You can also add column color markers to distinguish the two groups.
colAnnot.Labels = TNDiffExprSet.sampleNames;
colAnnot.Colors = {'g','g','g','g','g','g','g','g','r','r','r','r', 'r'};
%%
set(cg, 'ColumnLabelsColor', colAnnot, 'LabelsWithMarkers', true)

%%
% Clustering of the most differentially abundant transcripts clearly
% partitions teratozoonspermic (|Tz|) and normospermic (|Ns|) spermatozoal
% RNAs.

%% References
% [1] Platts, D.A., Le, J. M. et al. (2007) Success and failure in human
%     spermatogenesis as revealed by teratozoospermic RNAs. Human Molecular
%     Genetics, 16(7), 763-773.
%%
% [2] Storey, J.D., and Tibshirani, R. (2003). Statistical significance for
%     genomewide studies. Proc.Nat.Acad.Sci., 100(16), 9440-9445.

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20ILLUMINAGENEEXPDEMO%20in%20Bioinformatics%20Toolbox%2024.2 Provide feedback for this example.>*

displayEndOfDemoMessage(mfilename)

