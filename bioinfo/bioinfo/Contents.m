% Bioinformatics Toolbox
% Version 4.3.1 (R2013b) 08-Aug-2013
%
% File I/O
%   BioIndexedFile   - Read access to large text files using an index file.
%   affyread         - Read Affymetrix GeneChip files. 
%   affyprobeseqread - Read Affymetrix GeneChip probe sequence file.
%   agferead         - Read Agilent Feature Extraction format data.
%   bamindexread     - Read the index of a BAM formatted file.
%   baminfo          - Return a summary of the content of a BAM file.
%   bamread          - Read a BAM formatted file.
%   blastread        - Read an NCBI BLAST format report file.
%   blastreadlocal   - Read local BLAST format report file.
%   celintensityread - Read probe intensities from Affymetrix CEL files.
%   cytobandread     - Read cytogenetic banding information.
%   emblread         - Read an EMBL format file.
%   fastainfo        - Return a summary of the contents of a FASTA file.
%   fastaread        - Read a sequence from a FASTA format file or URL.
%   fastawrite       - Write a sequence to a FASTA format file.
%   fastqinfo        - Return a summary of the contents of a FASTQ file.
%   fastqread        - Read a FASTQ format file.
%   fastqwrite       - Write sequences to a FASTQ format file.
%   galread          - Read GenePix GAL file.
%   genbankread      - Read a GenBank format file.
%   genpeptread      - Read a GenPept format file.
%   geoseriesread    - Read Gene Expression Omnibus (GEO) GSE format data.
%   geosoftread      - Read Gene Expression Omnibus (GEO) SOFT format data.
%   getblast         - Get a BLAST report from NCBI.
%   getembl          - Get sequence data from EMBL.
%   getgenbank       - Get sequence data from GenBank.
%   getgenpept       - Get sequence data from GenPept.
%   getgeodata       - Get Gene Expression Omnibus (GEO) data.
%   gethmmalignment  - Get a multiple alignment from the PFAM database.
%   gethmmprof       - Get a HMM from the PFAM database.
%   gethmmtree       - Get a phylogenetic tree from the PFAM database.
%   getpdb           - Get sequence data from PDB.
%   gprread          - Read GenePix GPR file.
%   ilmnbsread       - Read data exported from Illumina BeadStudio. 
%   imageneread      - Read ImaGene format results file.
%   jcampread        - Read JCAMP-DX file.
%   multialignread   - Read a multiple sequence alignment file.
%   multialignwrite  - Write a multiple sequence alignment file.
%   mzxmlinfo        - Return information about mzXML file 
%   mzxmlread        - Read mzXML file.
%   pdbread          - Read a PDB format file.
%   pdbwrite         - Write a PDB format file.
%   pfamhmmread      - Read a PFAM format HMM profile.
%   phytreeread      - Read NEWICK tree formatted file.
%   saminfo          - Return a summary of the content of a SAM file.
%   samread          - Read a SAM formatted file.
%   soapread         - Read a SOAP aligner formatted file.
%   scfread          - Read an SCF format trace file.
%   sffinfo          - Return a summary of the content of an SFF file.
%   sffread          - Read an SFF format file.
%   sptread          - Read SPOT format file.
%   tgspcinfo        - Return information about an SPC format file.
%   tgspcread        - Reads a Thermo-Galactic SPC format file.
%
% Sequence Conversion
%   aa2int          - Convert from amino acid to integer representation.
%   aa2nt           - Convert a sequence of amino acids to nucleotides.
%   dna2rna         - Convert a sequence of DNA nucleotides to RNA.
%   int2aa          - Convert from integer to amino acid representation.
%   int2nt          - Convert from integer to nucleotide representation.
%   nt2aa           - Convert a sequence of nucleotides to amino acids.
%   nt2int          - Convert from nucleotide to integer representation.
%   rnaconvert      - Convert RNA structure between bracket and matrix representation.
%   rna2dna         - Convert a sequence of RNA nucleotides to DNA.
%   seq2regexp      - Convert a sequence that contains wildcards to a regular expression.
%   seqcomplement   - Calculate the complementary strand of a DNA sequence.
%   seqrcomplement  - Calculate the reverse complement of a DNA sequence.
%   seqreverse      - Reverse a sequence.
%
% Sequence Statistics
%   aacount         - Report amino acid counts in a sequence.
%   atomiccomp      - Calculate atomic composition of a protein.
%   basecount       - Report nucleotide base counts in a sequence.
%   BioReadQualityStatistics - quality statistics from short-read sequences.
%   codonbias       - Report codon usage per amino acid for a DNA sequence.
%   codoncount      - Report codon counts in a sequence.
%   cpgisland       - Locate CpG islands in a DNA sequence.
%   dimercount      - Report dimer counts in a sequence.
%   isoelectric     - Estimate the isoelectric point of a protein sequence.
%   molweight       - Calculate molecular weight of a peptide sequence.
%   nmercount       - Report n-mer counts in a sequence.
%   ntdensity       - Plot nucleotide density along the sequence.
%   seqwordcount    - Report word counts for a sequence.
%
% Sequence Utilities
%   aminolookup      - Lookup table for peptide symbols.
%   baselookup       - Lookup table for nucleotide symbols.
%   cleave           - Cleave a protein with an enzyme.
%   evalrasmolscript - Send Rasmol script to a molecule viewer.
%   featuresmap      - Graphical map showing the features of a GenBank structure.
%   featuresparse    - Parse features from GenBank, GenPept, or, EMBL data.
%   geneticcode      - Mapping for the genetic code.
%   joinseq          - Join two sequences.
%   ngsbrowser       - Browse short-read sequence alignment.
%   molviewer        - Visualize molecules.
%   oligoprop        - DNA oligonucleotide sequence properties.
%   palindromes      - Find palindromes in a sequence.
%   pdbdistplot      - Visualization of inter-molecular distances in a PDB file.
%   proteinplot      - GUI for protein analysis.
%   ramachandran     - Ramachandran plot for PDB data.
%   randseq          - Generate a random sequence from a finite alphabet.
%   rebasecuts       - Find restriction enzymes that cut a sequence.
%   restrict         - Split a sequence at a restriction site.
%   revgeneticcode   - Reverse mapping for the genetic code.
%   rnaplot          - Plot RNA secondary structure.
%   rnafold          - Predict secondary structure of a RNA sequence.
%   seqconsensus     - Compute the consensus sequence for a set of sequences.      
%   seqdisp          - Format long sequences for easy viewing.
%   seqlogo          - Display sequence logos for DNA and protein sequences.
%   seqmatch         - Find matches for every string in a library.
%   seqprofile       - Compute the sequence profile of a multiple alignment.
%   seqshoworfs      - Graphical display of Open Reading Frames in a sequence.
%   seqshowwords     - Graphical display of words in a sequence.
%   seqviewer        - Visualize biological sequences.
%
% Sequence Alignment
%   align2cigar      - Calculate the CIGAR of aligned sequences.
%   blastformat      - Run local version of BLAST formatdb.
%   blastlocal       - Run local version of BLAST.
%   blastncbi        - Generate a remote NCBI BLAST request.
%   bowtie           - Map short reads using the Burrows-Wheeler transform.
%   bowtiebuild      - Generate an index using Burrows-Wheeler transform.
%   cigar2align      - Align sequences using a CIGAR string.
%   localalign       - Find optimal and suboptimal local alignments of two sequences.
%   multialign       - Progressive multiple sequence alignment.
%   ngsbrowser       - Interactive browser for exploring short read sequence alignment data.
%   nwalign          - Needleman-Wunsch global alignment.
%   profalign        - Needleman-Wunsch global alignment of two profiles.
%   seqalignviewer   - Visualize and edit a sequence alignment.
%   seqdotplot       - Create a dotplot of two sequences.
%   showalignment    - Visualization of pairwise sequence alignment.
%   swalign          - Smith-Waterman local alignment.
%
% Sequence Data Containers
%   BioRead          - Class representing a collection of sequences with quality scores.
%   BioMap           - Class representing a collection of sequences with alignment information.
%   GFFAnnotation    – Class representing a collection of GFF annotations.
%   GTFAnnotation    – Class representing a collection of GTF annotations.
%
% Statistical Learning
%   classify         - Discriminant analysis. (Statistics toolbox)
%   classperf        - Evaluate the performance of a classifier.
%   crossvalind      - Cross-validation index generation.
%   kmeans           - K-means clustering. (Statistics toolbox)
%   knnclassify      - K-Nearest neighbor classifier.
%   knnimpute        - Impute missing data using the nearest neighbor method.
%   optimalleaforder - Reorders a hierarchical binary cluster tree. (Statistics toolbox)
%   randfeatures     - Randomized subset feature selection.
%   rankfeatures     - Ranks key features by class separability criteria.
%   samplealign      - Aligns two data sets containing sequential observations. 
%   svmclassify      - Classify using a support vector machine classifier (Statistics toolbox). 
%   svmtrain         - Train a support vector machine classifier (Statistics toolbox). 
%   treefit          - Classification tree fitting (Statistics toolbox).
%
% Protein Analysis
%   aacount          - Show the amino acid composition of a protein sequence.
%   aminolookup      - Lookup table for peptide symbols.
%   atomiccomp       - Calculate atomic composition of a protein.
%   cleave           - Cleave a protein with an enzyme.
%   cleavelookup     - Display cleavage rules of enzymes or compounds.
%   evalrasmolscript - Send Rasmol script to a molecule viewer.
%   isoelectric      - Estimate the isoelectric point of a protein sequence.
%   molviewer        - Visualize molecules.
%   molweight        - Calculate molecular weight of a peptide sequence.
%   pdbdistplot      - Visualization of inter-molecular distances in a PDB file.
%   pdbsuperpose     - Superpose the 3-D structures of two proteins.
%   pdbtransform     - Apply a linear transformation to the 3D structure of a molecule.
%   proteinpropplot  - Plot hydrophobicity and other properties of a sequence.
%   proteinplot      - GUI for protein analysis.
%   ramachandran     - Ramachandran plot for PDB data.
%
% Trace tools
%   scfread         - Read SCF format trace data.
%   traceplot       - View nucleotide trace plots.
%
% Profile Hidden Markov Models
%   gethmmalignment - Get a multiple alignment from the PFAM database.
%   gethmmprof      - Get a HMM from the PFAM database.
%   gethmmtree      - Get a phylogenetic tree from the PFAM database.
%   hmmprofalign    - Sequence alignment to a profile HMM.
%   hmmprofestimate - Estimate the parameters of a profile HMM.
%   hmmprofgenerate - Generate a random sequence from a profile HMM.
%   hmmprofmerge    - Align the output strings of several profile alignments.
%   hmmprofstruct   - Create a profile HMM structure.
%   pfamhmmread     - Read a PFAM format HMM profile.
%   showhmmprof     - Plot an HMM profile.
%
% Phylogenetic Tree Tools
%   dnds            - Estimate synonymous and nonsynonymous substitution rates.
%   dndsml          - DNDS using maximum likelihood.
%   phytree         - Class representing a phylogenetic tree object.
%   phytreeread     - Read NEWICK tree formatted file.
%   phytreeviewer   - Visualize and edit phylogenetic trees.
%   phytreewrite    - Save a phylogenetic tree object as a NEWICK format file.
%   seqlinkage      - Construct a phylogenetic tree from pairwise distances.
%   seqneighjoin    - Neighbor-joining for phylogenetic tree reconstruction.
%   seqpdist        - Pairwise distance between sequences.
%
% Phylogenetic Tree Methods
%   phytree/cluster      - Construct clusters from a phylogenetic tree.
%   phytree/get          - Get information about a phylogenetic tree object.
%   phytree/getbyname    - Select branches and leaves by name.
%   phytree/getcanonical - Calculates the canonical form of a phylogenetic tree.
%   phytree/getmatrix    - Converts a tree into a relationship matrix.
%   phytree/getnewickstr - Creates a NEWICK formatted string.
%   phytree/pdist        - Compute the pairwise patristic distance.
%   phytree/plot         - Render a phylogenetic tree.
%   phytree/prune        - Reduce a phylogenetic tree by removing branch nodes.
%   phytree/reroot       - Changes the root of a phylogenetic tree.
%   phytree/reorder      - Changes the leaf order of a phylogenetic tree.
%   phytree/select       - Select tree leaves and branches.
%   phytree/subtree      - Extracts a subtree.
%   phytree/view         - View a phylogenetic tree in phytreeviewer.
%   phytree/weights      - Tree-based sequence weights.
%
% Microarray Data Containers
%   bioma.ExpressionSet     - Class to contain microarray gene expression experiment data.
%   bioma.data.DataMatrix   - Two dimensional data array with row and column names.
%   bioma.data.ExptData     - Data container class to store experiment data values.
%   bioma.data.MetaData     - Metadata collection of variables and their description.
%   bioma.data.MIAME        - Class for storing a description about a microarray experiment.
%
% Microarray Data Analysis and Visualization
%   cghcbs          - Compute circular binary segmentation on array CGH data.
%   cghfreqplot     - Frequency plot of copy number alterations.
%   chromosomeplot  - Plot chromosome ideograms.
%   clustergram     - Clustergram plot.
%   HeatMap         - False color image of numeric array.
%   maboxplot       - Box plot of microarray data.
%   mafdr           - Compute false discovery rates of gene expression data.
%   maimage         - Pseudocolor plot of microarray spatial data.
%   mairplot        - Intensity plot of microarray signals.
%   maloglog        - Log-log plot of microarray data.
%   mapcaplot       - Principal Component plot of expression profile data.
%   mattest         - Unpaired student's t-test of microarray data.
%   mavolcanoplot   - Volcano plot of expression profile data.
%   microplateplot  - Creates a visualization of a microtiter plate.
%   redbluecmap     - Generate red and blue colormap.
%   redgreencmap    - Generate red and green colormap.
%
% Microarray Normalization and Filtering
%   affygcrma         - Perform GCRMA procedure for multiple Affymetrix GeneChips.
%   affyinvarsetnorm  - Invariant set normalization of Affymetrix probe-level data.
%   affyprobeaffinities - Compute probe affinities of Affymetrix GeneChip.
%   affyrma           - Perform RMA procedure for multiple Affymetrix GeneChips.
%   exprprofrange     - Calculate range of expression profiles.
%   exprprofvar       - Calculate variance of expression profiles.
%   gcrma             - GCRMA measure gene expression of Affymetrix microarray data.
%   gcrmabackadj      - GCRMA background adjustment of Affymetrix probe-level data.
%   geneentropyfilter - Remove genes with entropy expression values.
%   genelowvalfilter  - Remove genes with low expression values.
%   generangefilter   - Remove genes with small expression ranges.
%   genevarfilter     - Remove genes with small expression variance.
%   mainvarsetnorm    - Rank invariant set normalization.
%   malowess          - Lowess normalization.
%   manorm            - Normalization by scaling and centering.
%   quantilenorm      - Quantile normalization.
%   rmabackadj        - RMA background adjustment of Affymetrix probe-level data.
%   rmasummary        - RMA summarization of multiple Affymetrix microarray data.
%   zonebackadj       - Zone based background adjustment of Affymetrix probe-level data.
%
% Microarray Utility Functions
%   affysnpintensitysplit - Splits Affymetrix SNP probe intensity matrix for allele A and B.
%   affysnpquartets       - Create a table of SNP quartets for a probe set.
%   ilmnbslookup          - Look up Illumina probe sequence and annotation information.
%   magetfield            - Extract data from microarray structure.
%   probelibraryinfo      - Get library information for a probe.
%   probesetlink          - Show probe set information from NetAffx.
%   probesetlookup        - Get gene information for a probe set.
%   probesetplot          - Plot probe set values.
%   probesetvalues        - Get probe set values from CEL and CDF information.
%
% Gene Ontology Functions and Methods
%   geneont                        - Creates a Gene Ontology (GO) object
%   geneont.geneont.getancestors   - Finds the ancestors of a GO term
%   geneont.geneont.getdescendants - Finds the descendents of a GO term
%   geneont.geneont.getmatrix      - Converts a GO Object into a relationship matrix
%   geneont.geneont.getrelatives   - Finds the related terms for a GO term
%   goannotread                    - Extract data from microarray structure.
%   num2goid                       - Converts numeric values to GO IDs
%
% Bioanalytics and Mass-Spectrometry Preprocessing and Visualization
%   isotopicdist      - Calculate isotope mass distribution and density function.
%   jcampread         - Read JCAMP-DX file.
%   msalign           - Signal calibration and alignment by reference peaks.
%   msbackadj         - Background estimation and correction.
%   msdotplot         - Create a dot plot of an LCMS or GCMS dataset.
%   msheatmap         - Heat map image of a set of spectra.
%   mslowess          - Non-parametric smoothing using Lowess method.
%   msnorm            - Normalization of a set of spectra.
%   mspalign          - Peak binning and dynamic programming peak alignment.
%   mspeaks           - Peak detection with wavelet denoising.
%   msppresample      - Signal resampling from peak information.
%   msresample        - Resample with antialias filtering.
%   mssgolay          - Least-squares polynomial smoothing.
%   msviewer          - Plot a spectrum or a set of spectra.
%   mzcdf2peaks       - Convert an mzCDF structure to a list of peaks.
%   mzcdfinfo         - Return information about netCDF file.
%   mzcdfread         - Read netCDF file with mass-spectrometric data.
%   mzxml2peaks       - Convert an mzXML structure to a list of peaks.
%   mzxmlinfo         - Return information about mzXML file.
%   mzxmlread         - Read mzXML file.
%   samplealign       - Constrained dynamic programming alignment and warping.
%
% Graph Theory Algorithms
%   graphallshortestpaths - Find distance of all shortest paths.
%   graphconncomp         - Strong and weak connected components.
%   graphisdag            - Check if graph is DAG.
%   graphisomorphism      - Map between two isomorphic graphs.
%   graphisspantree       - Check if graph is a spanning tree.
%   graphmaxflow          - Max-flow (and min-cut) algorithm.
%   graphminspantree      - Find the minimal spanning tree.
%   graphpred2path        - Covert from a predecessor list to a path.
%   graphshortestpath     - Find the shortest path between two nodes.
%   graphtopoorder        - Topological order of a DAG.
%   graphtraverse         - Depth first search and breadth first search.
% 
% Graph Visualization Methods
%   biograph                           - Create a bioinformatics graph object.
%   biograph.biograph.dolayout         - Calculate node and edge positions.
%   biograph.biograph.getmatrix        - Get the relationship matrix.
%   biograph.biograph.getnodesbyid     - Get handles to nodes.
%   biograph.biograph.getedgesbynodeid - Get handles to edges.
%   biograph.biograph.view             - Render a graph in its viewer.
%   biograph.node.getancestors         - Find ancestors.
%   biograph.node.getdescendants       - Find descendants.
%   biograph.node.getrelatives         - Find neighbors.
%  
% Scoring Matrices
%   blosum            - BLOSUM family of matrices.
%   dayhoff           - Dayhoff matrix.
%   gonnet            - Gonnet variation on PAM250.
%   nuc44             - Nuc44 nucleotide matrix.
%   pam               - PAM family of matrices.
%
% Tutorials, demos and examples.
%   acghhmmdemo        - Bayesian hidden Markov modeling of array CGH data.
%   affydemo           - Example of working with Affymetrix GeneChip data.
%   affypreprocessdemo - Preprocessing Affymetrix microarray data at the probe level.
%   affysnpcnvdemo     - Analyzing Affymetrix SNP arrays for DNA copy number variants.
%   aligndemo          - Basic sequence alignment tutorial demo. 
%   alignscoringdemo   - Tutorial showing the use of scoring matrices. 
%   alignsigdemo       - Demo of how to estimate the significance of alignments.
%   bacacghdemo        - Detecting DNA copy number alteration in array-based CGH data.
%   birdfludemo        - Investigating the Bird Flu virus (H5N1).
%   biodbdemo          - Example of connecting to local databases.
%   biodistcompdemo    - Batch processing through sequential and parallel computing.
%   biographdemo       - Working with BIOGRAPH objects.
%   biomemorymapdemo   - Using memory mapping to work with whole genome data.
%   bioperldemo        - Example of calling Bioperl functions.
%   cancerdetectdemo   - Data mining analysis for mass spectrometry profiles.
%   chipseqpedemo      - Exploring protein-DNA binding sites from paired-end ChIP-seq data.
%   clustergramdemo    - Clustergram functionality examples.
%   cnsgeneexpdemo     - Exploring gene expression data.
%   diffprotdemo       - Differential proteomics and metabolomics analysis of LCMS.
%   evemotifdemo       - Identifying over-represented regulatory motifs. 
%   dndsdemo           - Analyzing synonymous and non-synonymous substitution rates.
%   geneontologydemo   - Example of working with Gene Ontology data.
%   graphtheorydemo    - Example of working with graph theoretic functions.
%   gsedemo            - Example of working with Gene Expression Omnibus Series data.
%   gutmicrobiomedemo  - Analyzing the human distal gut microbiome.
%   hmmprofdemo        - HMM profile alignment tutorial example.
%   hivdemo            - Analyzing the origin of the HIV with phylogenetic trees.
%   illuminagedemo     - Example of working with Illumina BeadChip data.
%   ilmnsolexademo     - Working with Illumina/Solexa next-generation sequencing data.
%   lcmsdemo           - Visualizing and preprocessing LCMS data.
%   maexptdemo         - Examples of working with microarray experiment data structures.
%   mbdseqdemo         - Exploring genome-wide differences in DNA methylation profiles.
%   metagenomicdemo    - Metagenomic analysis of a Sargasso Sea Sample.
%   molviewerdemo      - Visualizing the three-dimensional structure of a molecule.
%   mousedemo          - Microarray normalization and visualization example.
%   msgademo           - Mass spectra data analysis with Genetic Algorithms.
%   mspreprodemo       - Preprocessing of raw mass spectrometry data.
%   ncbieutilsdemo     - Accessing NCBI Entrez Databases with E-Utilities.
%   phybootdistdemo    - Confidence estimation of trees by bootstrapping.
%   primerdemo         - Primer design tutorial example.
%   primatesdemo       - Building a phylogenetic tree for the hominidae species.
%   rnademo            - Predicting and visualizing the secondary structure of RNA.
%   rnaseqdedemo       - Identifying differentially expressed genes from RNA-seq data.
%   sarsdemo           - Reconstructing the origin and the diffusion of the SARS epidemic.
%   secstructnnetdemo  - Predicting protein secondary structure using a neural network.
%   seqstatsdemo       - Sequence statistics tutorial example.
%   sffdemo            - Working with SFF files from 454 Genome Sequencer.
%   wholegenomedemo    - Comparing whole genomes.
%   yeastdemo          - Microarray data analysis example.
%

% Utility Functions Without Reference Pages
%   bioinfochecknargin - Check the number of input arguments to a function.

%   GenePix is a registered trademark of Axon Instruments, Inc. 
%   GeneChip and Affymetrix are registered trademark of Affymetrix, Inc. 
%   Illumina is a registered trademark of Illumina, Inc.

%   Copyright 2003-2013 The MathWorks, Inc. 


