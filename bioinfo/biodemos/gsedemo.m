%% Working with GEO Series Data
% This example shows how to retrieve gene expression data series (GSE) from
% the NCBI Gene Expression Omnibus (GEO) and perform basic analysis on the
% expression profiles.

% Copyright 2008-2012 The MathWorks, Inc.

%% Introduction
% The NCBI Gene Expression Omnibus (<http://www.ncbi.nlm.nih.gov/geo/ GEO>)
% is the largest public repository of high-throughput microarray
% experimental data.  GEO data have four entity types including GEO
% Platform (GPL), GEO Sample (GSM), GEO Series (GSE) and curated GEO
% DataSet (GDS). 
%%
% A Platform record describes the list of elements on the array in the
% experiment (e.g., cDNAs, oligonucleotide probesets). Each Platform record
% is assigned a unique and stable GEO accession number (GPLxxx). 
%%
% A Sample record describes the conditions under which an individual Sample
% was handled, the manipulations it underwent, and the abundance
% measurement of each element derived from it. Each Sample record is
% assigned a unique and stable GEO accession number (GSMxxx).
%%
% A Series record defines a group of related Samples and provides a focal
% point and description of the whole study. Series records may also contain
% tables describing extracted data, summary conclusions, or analyses. Each
% Series record is assigned a unique GEO accession number (GSExxx).
%%
% A DataSet record (GDSxxx) represents a curated collection of
% biologically and statistically comparable GEO Samples. GEO DataSets
% (GDSxxx) are curated sets of GEO Sample data.   
%%
% More information about GEO can be found in
% <http://www.ncbi.nlm.nih.gov/projects/geo/info/overview.html GEO
% Overview>. Bioinformatics Toolbox provides functions that can retrieve
% and parse GEO format data files. GSE, GSM, GSD and GPL data can be
% retrieved by using the |getgeodata| function. The |getgeodata| function
% can also save the retrieved data in a text file. GEO Series records are
% available in SOFT format files and in tab-delimited text format files.
% The function |geoseriesread| reads the GEO Series text format file. The
% |geosoftread| function reads the usually quite large SOFT format files. 
% 
% In this example, you will retrieve the GSE5847 data set from GEO
% database, and perform statistical testing on the data. GEO Series GSE5847 
% contains experimental data from a gene expression study of tumor stroma
% and epithelium cells from 15 inflammatory breast cancer (IBC) cases and
% 35 non-inflammatory breast cancer cases (Boersma et al. 2008).

%% Retrieving GEO Series Data
% The function |getgeodata| returns a structure containing data retrieved
% from the GEO database. You can also save the returned data in its
% original format to your local file system for use in subsequent MATLAB(R)
% sessions, and use the |geoseriesread| function to parse the GSE text
% format file.
gseData = getgeodata('GSE5847', 'ToFile', 'GSE5847.txt')

%%
% The structure returned contains a |Header| field that stores the metadata
% of the Series data, and a |Data| field that stores the Series matrix
% data.

%% Exploring GSE Data
% The GSE5847 matrix data in the |Data| field are stored as a DataMatrix
% object. A DataMatrix object is a data structure like MATLAB 2-D array,
% but with additional metadata of row names and column names. The
% properties of a DataMatrix can be accessed like other MATLAB objects.
get(gseData.Data)

%%
% The row names are the identifiers of the probe sets on the array; the
% column names are the GEO Sample accession numbers. 
gseData.Data(1:5, 1:5)

%%
% The Series metadata are stored in the |Header| field.  The structure
% contains Series information in the |Header.Series| field, and sample
% information in the |Header.Sample| field.
gseData.Header

%%
% The Series field contains the title of the experiment and the microarray
% GEO Platform ID.
gseData.Header.Series

%%
gseData.Header.Samples

%%
% The |data_processing| field contains the information of the preprocessing
% methods, in this case the Robust Multi-array Average (RMA) procedure. 
gseData.Header.Samples.data_processing(1)

%%
% The |source_name_ch1| field contains the sample source:
sampleSources = unique(gseData.Header.Samples.source_name_ch1);
%%
sampleSources{:}

%%
% The field |Header.Samples.characteristics_ch1| contains the
% characteristics of the samples.
gseData.Header.Samples.characteristics_ch1(:,1)
%%
% Determine the IBC and non-IBC labels for the samples from the
% |Header.Samples.characteristics_ch1| field as group labels.
sampleGrp = gseData.Header.Samples.characteristics_ch1(1,:);

%% Retrieving GEO Platform (GPL) Data
% The Series metadata told us the array platform id: GPL96, which is an
% Affymetrix(R) GeneChip(R) Human Genome U133 array set HG-U133A. Retrieve
% the GPL96 SOFT format file from GEO using the |getgeodata| function. For
% example, assume you used the |getgeodata| function to retrieve the GPL96
% Platform file and saved it to a file, such as |GPL96.txt|. Use the
% |geosoftread| function to parse this SOFT format file.
gplData = geosoftread('GPL96.txt')

%%
% The |ColumnNames| field of the |gplData| structure contains names of the
% columns for the data:
gplData.ColumnNames
%%
% You can get the probe set ids and gene symbols for the probe sets of
% platform GPL69. 
gplProbesetIDs = gplData.Data(:, strcmp(gplData.ColumnNames, 'ID'));
geneSymbols = gplData.Data(:, strcmp(gplData.ColumnNames, 'Gene Symbol'));
%%
% Use gene symbols to label the genes in the DataMatrix object
% |gseData.Data|. Be aware that the probe set |IDs| from the platform file
% may not be in the same order as in |gseData.Data|. In this example they
% are in the same order.
%%
% Change the row names of the expression data to gene symbols. 
gseData.Data = rownames(gseData.Data, ':', geneSymbols);
%%
% Display the first five rows and five columns of the expression data
% with row names as gene symbols. 
gseData.Data(1:5, 1:5)

%% Analyzing the Data
% Bioinformatics Toolbox and Statistics Toolbox(TM) offer a wide spectrum
% of analysis and visualization tools for microarray data analysis.
% However, because it is not our main goal to explain the analysis methods
% in this example, you will apply only a few of the functions to the gene
% expression profile from stromal cells. For more elaborate examples
% about the gene expression analysis and feature selection tools available,
% see <cnsgeneexpdemo.html Exploring Gene Expression Data> and
% <http://www.mathworks.com/products/statistics/demos.html?file=/products/demos/shipping/stats/cvsequentialfsdemo.html
% Selecting Features for Classifying High-dimensional Data>.  

%%
% The experiment was performed on IBC and non-IBC samples derived from
% stromal cells and epithelial cells. In this example, you will work with
% the gene expression profile from stromal cells. Determine the sample
% indices for the stromal cell type from the
% |gseData.Header.Samples.source_name_ch1| field. 
stromaIdx = strcmpi(sampleSources{1}, gseData.Header.Samples.source_name_ch1);

%%
% Determine number of samples from stromal cells.
nStroma = sum(stromaIdx)

%%
stromaData = gseData.Data(:, stromaIdx);
stromaGrp = sampleGrp(stromaIdx);

%%
% Determine the number of IBC and non-IBC stromal cell samples.
nStromaIBC = sum(strcmp('IBC', stromaGrp))
%%
nStromaNonIBC = sum(strcmp('non-IBC', stromaGrp))

%%
% You can also label the columns in |stromaData| with the group labels:
stromaData = colnames(stromaData, ':', stromaGrp);

%%
% Display the histogram of the normalized gene expression measurements of a
% specified gene. The x-axes represent the normalized expression level.
% For example, inspect the distribution of the gene expression values of
% these genes.
fID = 331:339;
%%
zValues = zscore(stromaData.(':')(':'), 0, 2);
bw = 0.25;
edges = -10:bw:10;
bins = edges(1:end-1) + diff(edges)/2;
%%
histStroma = histc(zValues(fID, :)', edges) ./ (stromaData.NCols*bw);
%%
figure;
for i = 1:length(fID)
    subplot(3,3,i); 
    bar(edges, histStroma(:,i), 'histc')
    xlim([-3 3])
    if i <= length(fID)-3
        set(gca, 'XtickLabel', [])
    end
    title(sprintf('gene%d - %s', fID(i), stromaData.RowNames{fID(i)}))
end
suptitle('Gene Expression Value Distributions')
%% 
% The gene expression profile was accessed using the Affymetrix GeneChip
% more then 22,000 features on a small number of samples (~100). Among the
% 47 tumor stromal samples, there are 13 IBC and 34 non-IBC. But not all
% the genes are differentially expressed between IBC and non-IBC tumors.
% Statistical tests are needed to identify a gene expression signature that
% distinguish IBC from non-IBC stromal samples.
%%
% Use |genevarfilter| to filter out genes with a small variance across
% samples.
[mask, stromaData] = genevarfilter(stromaData);

%%
stromaData.NRows
%%
% Apply a t-statistic on each gene and compare _p-values_ for each gene to
% find significantly differentially expressed genes between IBC and non-IBC
% groups by permuting the samples (1,000 times for this example).
%%
randn('state', 0)
[pvalues, tscores]=mattest(stromaData(:, 'IBC'), stromaData(:, 'non-IBC'),...
                           'Showhist', true', 'showplot', true, 'permute', 1000);
%%
% Select the genes at a specified p-value.
sum(pvalues < 0.001)
%%
% There are about 50 genes selected directly at _p-values_ < 0.001. 

%%
% Sort and list the top 20 genes:
testResults = [pvalues, tscores];
testResults = sortrows(testResults);
testResults(1:20, :)

%% References
% [1] Boersma, B.J., Reimers, M., Yi, M., Ludwig, J.A., et al. "A stromal
%     gene signature associated with inflammatory breast cancer", Int J
%     Cancer 15;122(6), pp 1324-1332, 2008.

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20GSEDEMO%20in%20Bioinformatics%20Toolbox%204.2 Suggest an enhancement for this example.>*

displayEndOfDemoMessage(mfilename)
