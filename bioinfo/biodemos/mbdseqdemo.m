%% Exploring Genome-wide Differences in DNA Methylation Profiles
% This example shows how to perform a genome-wide analysis of DNA
% methylation in the human by using genome sequencing.

%  Copyright 2012 The MathWorks, Inc.

%%
% Note: For enhanced performance, MathWorks recommends that you run this
% example on a 64-bit platform, because the memory footprint is close to 2
% GB. On a 32-bit platform, if you receive |"Out of memory"| errors when
% running this example, try increasing the virtual memory (or swap space)
% of your operating system or try setting the |3GB| switch (32-bit
% Windows(R) XP only). These techniques are described in this
% <http://www.mathworks.com/support/tech-notes/1100/1107.html document>.

%% Introduction
% DNA methylation is an epigenetic modification that modulates gene
% expression and the maintenance of genomic organization in normal and
% disease processes. DNA methylation can define different states of the
% cell, and it is inheritable during cell replication. Aberrant DNA
% methylation patterns have been associated with cancer and tumor suppressor
% genes.
%
% In this example you will explore the DNA methylation profiles of two
% human cancer cells: parental HCT116 colon cancer cells and DICERex5
% cells. DICERex5 cells are derived from HCT116 cells after the truncation
% of the DICER1 alleles. Serre et al. in [1] proposed to study DNA
% methylation profiles by using the MBD2 protein as a methyl CpG binding
% domain and subsequently used high-throughput sequencing (HTseq). This
% technique is commonly know as MBD-Seq. Short reads for two replicates of 
% the two samples have been submitted to NCBI's SRA
% <http://www.ncbi.nlm.nih.gov/sra/ archive> by the authors of [1]. There
% are other technologies available to interrogate DNA methylation status of
% CpG sites in combination with HTseq, for example MeDIP-seq or the use of
% restriction enzymes. You can also analyze this type of data sets
% following the approach presented in this example.

%% Data Sets
% You can obtain the unmapped single-end reads for four sequencing
% experiments from the
% <ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP001/SRP001414/
% NCBI FTP site>. Short reads were produced using Illumina(R)'s Genome
% Analyzer II. Average insert size is 120 bp, and the length of short reads
% is 36 bp.

%%
% This example assumes that you:
% 
% (1) downloaded the files |SRR030222.sra|, |SRR030223.sra|,
% |SRR030224.sra| and |SRR030225.sra| containing the unmapped short 
% reads for two replicates of from the DICERex5 sample and two replicates
% from the HCT116 sample respectively. Converted them to FASTQ-formatted
% files using the <http://www.ncbi.nlm.nih.gov/books/NBK47540/ NCBI SRA
% Toolkit>. 
% 
% (2) produced SAM-formatted files by mapping the short reads to the
% reference human genome (NCBI Build 37.5) using the Bowtie [2] algorithm.
% Only uniquely mapped reads are reported.
%
% (3) compressed the SAM formatted files to BAM and ordered them by
% reference name first, then by genomic position by using SAMtools [3]. 


%%
% This example also assumes that you downloaded the reference human genome 
% (GRCh37.p5). You can use the |bowtie-inspect| command to reconstruct the
% human reference directly from the bowtie indices. Or you may download the
% reference from the NCBI repository by uncommenting the following line:

% getgenbank('NC_000009','FileFormat','fasta','tofile','hsch9.fasta');

%% Creating a MATLAB(R) Interface to the BAM-Formatted Files
% To explore the signal coverage of the HCT116 samples you need to
% construct a |BioMap|. |BioMap| has an interface that provides direct
% access to the mapped short reads stored in the BAM-formatted file, thus
% minimizing the amount of data that is actually loaded into memory.
% Use the function |baminfo| to obtain a list of the existing references
% and the actual number of short reads mapped to each one.

info = baminfo('SRR030224.bam','ScanDictionary',true);
fprintf('%-35s%s\n','Reference','Number of Reads');
for i = 1:numel(info.ScannedDictionary)
    fprintf('%-35s%d\n',info.ScannedDictionary{i},...
            info.ScannedDictionaryCount(i));
end

%%
% In this example you will focus on the analysis of chromosome 9. Create a
% |BioMap| for the two HCT116 sample replicates.

bm_hct116_1 = BioMap('SRR030224.bam','SelectRef','gi|224589821|ref|NC_000009.11|')
bm_hct116_2 = BioMap('SRR030225.bam','SelectRef','gi|224589821|ref|NC_000009.11|')

%%
% Using a binning algorithm provided by the |getBaseCoverage| method, you
% can plot the coverage of both replicates for an initial inspection. For
% reference, you can also add the ideogram for the human chromosome 9 to
% the plot using the |chromosomeplot| function. 

figure
ha = gca;
hold on
n = 141213431;               % length of chromosome 9
[cov,bin] = getBaseCoverage(bm_hct116_1,1,n,'binWidth',100);
h1 = plot(bin,cov,'b');      % plots the binned coverage of bm_hct116_1
[cov,bin] = getBaseCoverage(bm_hct116_2,1,n,'binWidth',100);
h2 = plot(bin,cov,'g');      % plots the binned coverage of bm_hct116_2
chromosomeplot('hs_cytoBand.txt', 9, 'AddToPlot', ha) % plots an ideogram along the x-axis
axis(ha,[1 n 0 100])         % zooms-in the y-axis
fixGenomicPositionLabels(ha) % formats tick labels and adds datacursors
legend([h1 h2],'HCT116-1','HCT116-2','Location','NorthEast')
ylabel('Coverage')
title('Coverage for two replicates of the HCT116 sample')
set(gcf,'Position',max(get(gcf,'Position'),[0 0 900 0])) % resize window

%%
% Because short reads represent the methylated regions of the DNA, there is
% a correlation between aligned coverage and DNA methylation. Observe the
% increased DNA methylation close to the chromosome telomeres; it is known
% that there is an association between DNA methylation and the role of
% telomeres for maintaining the integrity of the chromosomes. In the
% coverage plot you can also see a long gap over the chromosome centromere. 
% This is due to the repetitive sequences  present in the centromere, which
% prevent us from aligning short reads to a unique position in this region. 
% For the data sets used in this example, only about 30% of the short reads
% were uniquely mapped to the reference genome.

%% Correlating CpG Islands and DNA Methylation
% DNA methylation normally occurs in CpG dinucleotides. Alteration of the
% DNA methylation patterns can lead to transcriptional silencing,
% especially in the gene promoter CpG islands. But, it is also known that
% DNA methylation can block CTCF binding and can silence miRNA
% transcription among other relevant functions. In general, it is expected
% that mapped reads should preferably align to CpG rich regions. 
%
% Load the human chromosome 9 from the reference file |hs37.fasta|. For
% this example, it is assumed that you recovered the reference from the 
% Bowtie indices using the |bowtie-inspect| command; therefore |hs37.fasta|
% contains all the human chromosomes. To load only the chromosome 9 you can
% use the option nave-value pair |BLOCKREAD| with the |fastaread| function. 

chr9 = fastaread('hs37.fasta','blockread',9)

%% 
% Use the |cpgisland| function to find the CpG clusters. Using the standard  
% definition for CpG islands [4], 200 or more bp islands with 60% or
% greater CpGobserved/CpGexpected ratio, leads to 1682 GpG islands found in
% chromosome 9.

cpgi = cpgisland(chr9.Sequence)

%%
% Use the |getCounts| method to calculate the ratio of aligned bases that
% are inside CpG islands. For the first replicate of the sample HCT116, 
% the ratio is close to 45%. 

aligned_bases_in_CpG_islands = getCounts(bm_hct116_1,cpgi.Starts,cpgi.Stops,'method','sum')
aligned_bases_total = getCounts(bm_hct116_1,1,n,'method','sum')
ratio = aligned_bases_in_CpG_islands ./ aligned_bases_total

%%
% You can explore high resolution coverage plots of the two sample
% replicates and observe how the signal correlates with the CpG islands.
% For example, explore the region between 23,820,000 and 23,830,000 bp.
% This is the 5' region of the human gene ELAVL2.

r1 = 23820001; % set the region limits
r2 = 23830000;
fhELAVL2 = figure; % keep the figure handle to use it later
hold on
% plot high-resolution coverage of bm_hct116_1
h1 = plot(r1:r2,getBaseCoverage(bm_hct116_1,r1,r2,'binWidth',1),'b');
% plot high-resolution coverage of bm_hct116_2
h2 = plot(r1:r2,getBaseCoverage(bm_hct116_2,r1,r2,'binWidth',1),'g');

% mark the CpG islands within the [r1 r2] region 
for i=1:numel(cpgi.Starts)
   if cpgi.Starts(i)>r1 && cpgi.Stops(i)<r2 % is CpG island inside [r1 r2]?
      px = [cpgi.Starts([i i]) cpgi.Stops([i i])]; % x-coordinates for patch
      py = [0 max(ylim) max(ylim) 0];              % y-coordinates for patch
      hp = patch(px,py,'r','FaceAlpha',.1,'EdgeColor','r','Tag','cpgi');
   end
end

axis([r1 r2 0 20])            % zooms-in the y-axis
fixGenomicPositionLabels(gca) % formats tick labels and adds datacursors
legend([h1 h2 hp],'HCT116-1','HCT116-2','CpG Islands')
ylabel('Coverage')
xlabel('Chromosome 9 position')
title('Coverage for two replicates of the HCT116 sample')

%% Statistical Modelling of Count Data
% To find regions that contain more mapped reads than would be expected
% by chance, you can follow a similar approach to the one described by
% Serre et al. [1]. The number of counts for non-overlapping contiguous 100
% bp windows is statistically modeled.
%
% First, use the |getCounts| method to count the number of mapped reads
% that start at each window. In this example you use a binning approach
% that considers only the start position of every mapped read, following
% the approach of Serre et al. However, you may also use the |OVERLAP| and
% |METHOD| name-value pairs in |getCounts| to compute more accurate
% statistics. For instance, to obtain the maximum coverage for each window
% considering base pair resolution, set |OVERLAP| to 1 and |METHOD| to
% |MAX|. 

n = numel(chr9.Sequence); % length of chromosome
w = 1:100:n; % windows of 100 bp

counts_1 = getCounts(bm_hct116_1,w,w+99,'independent',true,'overlap','start');
counts_2 = getCounts(bm_hct116_2,w,w+99,'independent',true,'overlap','start');

%%
% First, try to model the counts assuming that all the windows with counts
% are biologically significant and therefore from the same distribution. 
% Use the negative bionomial distribution to fit a model the count data. 

nbp = nbinfit(counts_1);

%%
% Plot the fitted model over a histogram of the empirical data.

figure
hold on
emphist = histc(counts_1,0:100); % calculate the empirical distribution
bar(0:100,emphist./sum(emphist),'c','grouped') % plot histogram
plot(0:100,nbinpdf(0:100,nbp(1),nbp(2)),'b','linewidth',2); % plot fitted model
axis([0 50 0 .001])
legend('Empirical Distribution','Negative Binomial Fit')
ylabel('Frequency')
xlabel('Counts')
title('Frequency of counts for 100 bp windows (HCT116-1)')

%%
% The poor fitting indicates that the observed distribution may be due to
% the mixture of two models, one that represents the background and one
% that represents the count data in methylated DNA windows.
%
% A more realistic scenario would be to assume that windows with a small
% number of mapped reads are mainly the background (or null model). Serre
% et al. assumed that 100-bp windows contaning four or more reads are
% unlikely to be generated by chance. To estimate a good approximation to
% the null model, you can fit the left body of the emprirical distribution
% to a truncated negative binomial distribution. To fit a truncated
% distribution use the |mle| function. First you need to define an 
% anonymous function that defines the right-truncated version of |nbinpdf|.

rtnbinpdf = @(x,p1,p2,t) nbinpdf(x,p1,p2) ./ nbincdf(t-1,p1,p2);

%%
% Define the fitting function using another anonymous function.
rtnbinfit = @(x,t) mle(x,'pdf',@(x,p1,p2) rtnbinpdf(x,p1,p2,t),'start',nbinfit(x),'lower',[0 0]);

%%
% Before fitting the real data, let us assess the fiting procedure with
% some sampled data from a known distribution.

nbp = [0.5 0.2];              % Known coefficients
x = nbinrnd(nbp(1),nbp(2),10000,1); % Random sample
trun = 6;                     % Set a truncation threshold

nbphat1 = nbinfit(x);         % Fit non-truncated model to all data
nbphat2 = nbinfit(x(x<trun)); % Fit non-truncated model to truncated data (wrong)
nbphat3 = rtnbinfit(x(x<trun),trun); % Fit truncated model to truncated data

figure
hold on
emphist = histc(x,0:100);     % Calculate the empirical distribution
bar(0:100,emphist./sum(emphist),'c','grouped') % plot histogram
h1 = plot(0:100,nbinpdf(0:100,nbphat1(1),nbphat1(2)),'b-o','linewidth',2);
h2 = plot(0:100,nbinpdf(0:100,nbphat2(1),nbphat2(2)),'r','linewidth',2);
h3 = plot(0:100,nbinpdf(0:100,nbphat3(1),nbphat3(2)),'g','linewidth',2);
axis([0 25 0 .2])
legend([h1 h2 h3],'Neg-binomial fitted to all data',...
                  'Neg-binomial fitted to truncated data',...
                  'Truncated neg-binomial fitted to truncated data')
ylabel('Frequency')
xlabel('Counts')

%% Identifying Significant Methylated Regions
% For the two replicates of the HCT116 sample, fit a right-truncated
% negative binomial distribution to the observed null model using the
% |rtnbinfit| anonymous function previously defined.

trun = 4;  % Set a truncation threshold (as in [1])
pn1 = rtnbinfit(counts_1(counts_1<trun),trun); % Fit to HCT116-1 counts
pn2 = rtnbinfit(counts_2(counts_2<trun),trun); % Fit to HCT116-2 counts

%%
% Calculate the p-value for each window to the null distribution.
pval1 = 1 - nbincdf(counts_1,pn1(1),pn1(2));
pval2 = 1 - nbincdf(counts_2,pn2(1),pn2(2));

%%
% Calculate the false discovery rate using the |mafdr| function. Use the
% name-value pair |BHFDR| to use the linear-step up (LSU) procedure ([6])
% to calculate the FDR adjusted p-values. Setting the FDR < 0.01 permits
% you to identify the 100-bp windows that are significantly methylated.

fdr1 = mafdr(pval1,'bhfdr',true);
fdr2 = mafdr(pval2,'bhfdr',true);

w1 = fdr1<.01; % logical vector indicating significant windows in HCT116-1
w2 = fdr2<.01; % logical vector indicating significant windows in HCT116-2
w12 = w1 & w2; % logical vector indicating significant windows in both replicates

Number_of_sig_windows_HCT116_1 = sum(w1)
Number_of_sig_windows_HCT116_2 = sum(w2)
Number_of_sig_windows_HCT116 = sum(w12)

%%
% Overall, you identified 1662 and 1674 non-overlapping 100-bp windows in
% the two replicates of the HCT116 samples, which indicates there is
% significant evidence of DNA methylation. There are 1346 windows that are
% significant in both replicates. 
%
% For example, looking again in the promoter region of the ELAVL2 human
% gene you can observe that in both sample replicates, multiple 100-bp
% windows have been marked significant.

figure(fhELAVL2) % bring back to focus the previously plotted figure
plot(w(w1)+50,counts_1(w1),'bs') % plot significant windows in HCT116-1
plot(w(w2)+50,counts_2(w2),'gs') % plot significant windows in HCT116-2
axis([r1 r2 0 100])
title('Significant 100-bp windows in both replicates of the HCT116 sample')

%% Finding Genes With Significant Methylated Promoter Regions
% DNA methylation is involved in the modulation of gene expression. For
% instance, it is well known that hypermethylation is associated with the
% inactivation of several tumor suppresor genes. You can study in this
% data set the methylation of gene promoter regions. 
%
% First, download from Ensembl a tab-separated-value (TSV) table with all
% protein encoding genes to a text file, |ensemblmart_genes_hum37.txt|. For
% this example, we are using Ensamble release 64. Using Ensembl's
% <http://www.ensembl.org/biomart/martview/ BioMart> service, you can
% select a table with the following attributes: chromosome name, gene
% biotype, gene name, gene start/end, and strand direction.
%
% Use the provided helper function |ensemblmart2gff| to convert the
% downloaded TSV file to a GFF formatted file. Then use |GFFAnnotation| to
% load the file into MATLAB and create a subset with the genes present
% in chromosome 9 only. This results 800 annotated protein-coding genes in
% the Ensembl database.

GFFfilename = ensemblmart2gff('ensemblmart_genes_hum37.txt');
a = GFFAnnotation(GFFfilename)
a9 = getSubset(a,'reference','9')
numGenes = a9.NumEntries

%%
% Find the promoter regions for each gene. In this example we consider the
% proximal promoter as the -500/100 upstream region.

downstream = 500;
upstream   = 100;

geneDir = strcmp(a9.Strand,'+');  % logical vector indicating strands in the forward direction

% calculate promoter's start position for genes in the forward direction
promoterStart(geneDir) = a9.Start(geneDir) - downstream;
% calculate promoter's end position for genes in the forward direction
promoterStop(geneDir) = a9.Start(geneDir) + upstream;
% calculate promoter's start position for genes in the reverse direction
promoterStart(~geneDir) = a9.Stop(~geneDir) - upstream;
% calculate promoter's end position for genes in the reverse direction
promoterStop(~geneDir) = a9.Stop(~geneDir) + downstream;

%%
% Use a |dataset| as a container for the promoter information, as we can
% later add new columns to store gene counts and p-values.

promoters = dataset({a9.Feature,'Gene'});
promoters.Strand = char(a9.Strand);
promoters.Start = promoterStart';
promoters.Stop = promoterStop';

%%
%
% Find genes with significant DNA methylation in the promoter region by
% looking at the number of mapped short reads that overlap at least one
% base pair in the defined promoter region. 

promoters.Counts_1 = getCounts(bm_hct116_1,promoters.Start,promoters.Stop,...
                               'overlap',1,'independent',true);
promoters.Counts_2 = getCounts(bm_hct116_2,promoters.Start,promoters.Stop,...
                               'overlap',1,'independent',true);

%%
% Fit a null distribution for each sample replicate and compute the
% p-values: 
trun = 5;  % Set a truncation threshold
pn1 = rtnbinfit(promoters.Counts_1(promoters.Counts_1<trun),trun); % Fit to HCT116-1 promoter counts
pn2 = rtnbinfit(promoters.Counts_2(promoters.Counts_2<trun),trun); % Fit to HCT116-2 promoter counts
promoters.pval_1 = 1 - nbincdf(promoters.Counts_1,pn1(1),pn1(2)); % p-value for every promoter in HCT116-1
promoters.pval_2 = 1 - nbincdf(promoters.Counts_2,pn2(1),pn2(2)); % p-value for every promoter in HCT116-2

Number_of_sig_promoters =  sum(promoters.pval_1<.01 & promoters.pval_2<.01)

Ratio_of_sig_methylated_promoters = Number_of_sig_promoters./numGenes

%%
% Observe that only 74 (out of 800) genes in chromosome 9 have
% significantly DNA methylated regions (pval<0.01 in both replicates).
% Display a report of the 30 genes with the most significant methylated
% promoter regions. 

[~,order] = sort(promoters.pval_1.*promoters.pval_2);
promoters(order(1:30),[1 2 3 4 5 7 6 8])

%% Finding Intergenic Regions that are Significantly Methylated
% Serre et al. [1] reported that, in these data sets, approximately 90% of
% the uniquely mapped reads fall outside the 5' gene promoter regions.
% Using a similar approach as before, you can find genes that have
% intergenic methylated regions. To compensate for the varying lengths of
% the genes, you can use the maximum coverage, computed base-by-base,
% instead of the raw number of mapped short reads. Another alternative
% approach to normalize the counts by the gene length is to set the
% |METHOD| name-value pair to |rpkm| in the |getCounts| function.

intergenic = dataset({a9.Feature,'Gene'});
intergenic.Strand = char(a9.Strand);
intergenic.Start = a9.Start;
intergenic.Stop = a9.Stop;

intergenic.Counts_1 = getCounts(bm_hct116_1,intergenic.Start,intergenic.Stop,...
                      'overlap','full','method','max','independent',true);
intergenic.Counts_2 = getCounts(bm_hct116_2,intergenic.Start,intergenic.Stop,...
                      'overlap','full','method','max','independent',true);
trun = 10; % Set a truncation threshold
pn1 = rtnbinfit(intergenic.Counts_1(intergenic.Counts_1<trun),trun); % Fit to HCT116-1 intergenic counts
pn2 = rtnbinfit(intergenic.Counts_2(intergenic.Counts_2<trun),trun); % Fit to HCT116-2 intergenic counts
intergenic.pval_1 = 1 - nbincdf(intergenic.Counts_1,pn1(1),pn1(2)); % p-value for every intergenic region in HCT116-1
intergenic.pval_2 = 1 - nbincdf(intergenic.Counts_2,pn2(1),pn2(2)); % p-value for every intergenic region in HCT116-2

Number_of_sig_genes =  sum(intergenic.pval_1<.01 & intergenic.pval_2<.01)

Ratio_of_sig_methylated_genes = Number_of_sig_genes./numGenes

[~,order] = sort(intergenic.pval_1.*intergenic.pval_2);

intergenic(order(1:30),[1 2 3 4 5 7 6 8])

%%
% For instance, explore the methylation profile of the BARX1 gene, the
% sixth significant gene with intergenic methylation in the previous list.
% The GTF formatted file |ensemblmart_barx1.gtf| contains structural
% information for this gene obtained from Ensembl using the
% <http://www.ensembl.org/biomart/martview/ BioMart> service.
%
% Use |GTFAnnotation| to load the structural information into MATLAB. There
% are two annotated transcripts for this gene.

barx1 = GTFAnnotation('ensemblmart_barx1.gtf')
transcripts = getTranscriptNames(barx1)

%%
% Plot the DNA methylation profile for both HCT116 sample replicates with
% base-pair resolution. Overlay the CpG islands and plot the exons for each
% of the two transcripts along the bottom of the plot.

range = barx1.getRange;
r1 = range(1)-1000; % set the region limits
r2 = range(2)+1000;
figure
hold on
% plot high-resolution coverage of bm_hct116_1
h1 = plot(r1:r2,getBaseCoverage(bm_hct116_1,r1,r2,'binWidth',1),'b');
% plot high-resolution coverage of bm_hct116_2
h2 = plot(r1:r2,getBaseCoverage(bm_hct116_2,r1,r2,'binWidth',1),'g');

% mark the CpG islands within the [r1 r2] region 
for i=1:numel(cpgi.Starts)
    if cpgi.Starts(i)>r1 && cpgi.Stops(i)<r2 % is CpG island inside [r1 r2]?
       px = [cpgi.Starts([i i]) cpgi.Stops([i i])];  % x-coordinates for patch
       py = [0 max(ylim) max(ylim) 0];               % y-coordinates for patch
       hp = patch(px,py,'r','FaceAlpha',.1,'EdgeColor','r','Tag','cpgi');
    end
end

% mark the exons at the bottom of the axes
for i = 1:numel(transcripts)
    exons = getSubset(barx1,'Transcript',transcripts{i},'Feature','exon');
    for j = 1:exons.NumEntries
       px = [exons.Start([j j]);exons.Stop([j j])]'; % x-coordinates for patch
       py = [0 1 1 0]-i*2-1;                         % y-coordinates for patch
       hq = patch(px,py,'b','FaceAlpha',.1,'EdgeColor','b','Tag','exon'); 
    end
end

axis([r1 r2 -numel(transcripts)*2-2 80])  % zooms-in the y-axis
fixGenomicPositionLabels(gca) % formats tick labels and adds datacursors
ylabel('Coverage')
xlabel('Chromosome 9 position')
title('High resolution coverage in the BARX1 gene')
legend([h1 h2 hp hq],'HCT116-1','HCT116-2','CpG Islands','Exons','Location','NorthWest')

%%
% Observe the highly methylated region in the 5' promoter region
% (right-most CpG island). Recall that for this gene trasciption occurs in
% the reverse strand. More interesting, observe the highly methylated
% regions that overlap the initiation of each of the two annotated
% transcripts (two middle CpG islands).

%% Differential Analysis of Methylation Patterns
% In the study by Serre et al. another cell line is also analyzed. New
% cells (DICERex5) are derived from the same HCT116 colon cancer cells
% after truncating the DICER1 alleles. It has been reported in literature
% [5] that there is a localized change of DNA methylation at small number
% of gene promoters. In this example, you be find significant 100-bp
% windows in two sample replicates of the DICERex5 cells following the same 
% approach as the parental HCT116 cells, and then you will search
% statistically significant differences between the two cell lines.
%
% The helper function |getWindowCounts| captures the similar steps to find
% windows with significant coverage as before. |getWindowCounts| returns
% vectors with counts, p-values, and false discovery rates for each new
% replicate.

bm_dicer_1 = BioMap('SRR030222.bam','SelectRef','gi|224589821|ref|NC_000009.11|');
bm_dicer_2 = BioMap('SRR030223.bam','SelectRef','gi|224589821|ref|NC_000009.11|');
[counts_3,pval3,fdr3] = getWindowCounts(bm_dicer_1,4,w,100);
[counts_4,pval4,fdr4] = getWindowCounts(bm_dicer_2,4,w,100);
w3 = fdr3<.01; % logical vector indicating significant windows in DICERex5_1
w4 = fdr4<.01; % logical vector indicating significant windows in DICERex5-2
w34 = w3 & w4; % logical vector indicating significant windows in both replicates
Number_of_sig_windows_DICERex5_1 = sum(w3) 
Number_of_sig_windows_DICERex5_2 = sum(w4)
Number_of_sig_windows_DICERex5 = sum(w34)

%%
% To perform a differential analysis you use the 100-bp windows that are
% significant in at least one of the samples (either HCT116 or DICERex5).

wd = w34 | w12; % logical vector indicating windows included in the diff. analysis

counts = [counts_1(wd) counts_2(wd) counts_3(wd) counts_4(wd)];
ws = w(wd); % window start for each row in counts

%%
% Use the function |manorm| to normalize the data. The |PERCENTILE|
% name-value pair lets you filter out windows with very large number of
% counts while normalizing, since these windows are mainly due to
% artifacts, such as repetitive regions in the reference chromosome. 

counts_norm = round(manorm(counts,'percentile',90).*100);

%% 
% Use the function |mattest| to perform a two-sample t-test to identify
% differentially covered windows from the two different cell lines.
pval = mattest(counts_norm(:,[1 2]),counts_norm(:,[3 4]),'bootstrap',true,...
               'showhist',true,'showplot',true);

%%
% Create a report with the 25 most significant differentially covered
% windows. While creating the report use the helper function
% |findClosestGene| to determine if the window is intergenic, intragenic,
% or if it is in a proximal promoter region.
[~,ord] = sort(pval);
fprintf('Window Pos       Type                  p-value   HCT116     DICERex5\n\n');
for i = 1:25
    j = ord(i);
    [~,msg] = findClosestGene(a9,[ws(j) ws(j)+99]);
    fprintf('%10d  %-25s %7.6f%5d%5d %5d%5d\n',  ...
         ws(j),msg,pval(j),counts_norm(j,:));
end

%%
% Plot the DNA methylation profile for the promoter region of gene
% FAM189A2, the most signicant differentially covered promoter region from
% the previous list. Overlay the CpG islands and the FAM189A2 gene.

range = getRange(getSubset(a9,'Feature','FAM189A2'));
r1 = range(1)-1000;
r2 = range(2)+1000;
figure
hold on
% plot high-resolution coverage of all replicates
h1 = plot(r1:r2,getBaseCoverage(bm_hct116_1,r1,r2,'binWidth',1),'b');
h2 = plot(r1:r2,getBaseCoverage(bm_hct116_2,r1,r2,'binWidth',1),'g');
h3 = plot(r1:r2,getBaseCoverage(bm_dicer_1,r1,r2,'binWidth',1),'r');
h4 = plot(r1:r2,getBaseCoverage(bm_dicer_2,r1,r2,'binWidth',1),'m');

% mark the CpG islands within the [r1 r2] region 
for i=1:numel(cpgi.Starts)
    if cpgi.Starts(i)>r1 && cpgi.Stops(i)<r2 % is CpG island inside [r1 r2]?
       px = [cpgi.Starts([i i]) cpgi.Stops([i i])]; % x-coordinates for patch
       py = [0 max(ylim) max(ylim) 0];              % y-coordinates for patch
       hp = patch(px,py,'r','FaceAlpha',.1,'EdgeColor','r','Tag','cpgi');
    end
end

% mark the gene at the bottom of the axes
px = range([1 1 2 2]);
py = [0 1 1 0]-2;
hq = patch(px,py,'b','FaceAlpha',.1,'EdgeColor','b','Tag','gene');

axis([r1 r1+4000 -4 30]) % zooms-in 
fixGenomicPositionLabels(gca) % formats tick labels and adds datacursors
ylabel('Coverage')
xlabel('Chromosome 9 position')
title('DNA Methylation profiles along the promoter region of the FAM189A2 gene.')
legend([h1 h2 h3 h4 hp hq],'HCT116-1','HCT116-2','DICERex5-1','DICERex5-2','CpG Islands','FAM189A2 Gene','Location','NorthEast')

%%
% Observe that the CpG islands are clearly unmethylated for both of the
% DICERex5 replicates. 

%% References
% [1] Serre, D., Lee, B.H., and Ting A.H. "MBD-isolated Genome Sequencing
%     provides a high-throughput and comprehensive survey of DNA
%     methylation in the human genome", Nucleic Acids Research, 
%     38(2), pp 391-399, 2010.
%%
% [2] Langmead, B., Trapnell, C., Pop, M., and Salzberg, S.L. "Ultrafast
%     and Memory-efficient Alignment of Short DNA Sequences to the Human
%     Genome", Genome Biology, 10:R25, pp 1-10, 2009.
%%
% [3] Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N.,
%     Marth, G., Abecasis, G., Durbin, R. and 1000 Genome Project Data
%     Processing Subgroup "The Sequence Alignment/map (SAM) Format and
%     SAMtools", Bioinformatics, 25, pp 2078-2079, 2009.
%%
% [4] Gardiner-Garden, M and Frommer, M. "CpG islands in vertebrate
%     genomes", J.Mol.Biol. 196, pp 261-282, 1987. 
%%
% [5] Ting, A.H., Suzuki, H., Cope, L., Schuebel, K.E., Lee, B.H.,
%     Toyota, M., Imai, K., Shinomura, Y., Tokino, T. and Baylin, S.B. "A
%     Requirement for DICER to Maintain Full Promoter CpG Island %
%     Hypermethylation in Human Cancer Cells", Cancer Research, 68, 2570,
%     April 15, 2008.
%
%%
% [6] Benjamini, Y., Hochberg, Y., "Controlling the false discovery rate: a
%     practical and powerful approach to multiple testing", Journal of the
%     Royal Statistical Society, 57, pp 289-300, 1995.
%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20MBDSEQDEMO%20in%20Bioinformatics%20Toolbox%204.2 Provide feedback for this example.>*

displayEndOfDemoMessage(mfilename)
