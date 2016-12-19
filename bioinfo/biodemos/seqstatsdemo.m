%% Calculating and Visualizing Sequence Statistics
% This example shows basic sequence manipulation techniques and computes 
% some useful sequence statistics. It also illustrates how to look for 
% coding regions (such as proteins) and pursue further analysis of them. 

%   Copyright 2002-2013 The MathWorks, Inc.


%% The Human Mitochondrial Genome
% In this example you will explore the DNA sequence of the human 
% mitochondria. Mitochondria are structures, called organelles, that are
% found in the cytoplasm of the cell in hundreds to thousands for each
% cell. Mitochondria are generally the major energy production center in
% eukaryotes, they help to degrade fats and sugars. 

%%
% The consensus sequence of the human mitochondria genome has accession
% number NC_012920. You can |getgenbank| function to get the latest
% annotated sequence from GenBank(R) into the MATLAB(R) workspace.
%
%   mitochondria_gbk = getgenbank('NC_012920');

%%
% For your convenience, previously downloaded sequence is included in a
% MAT-file. Note that data in public repositories is frequently curated and 
% updated; therefore the results of this example might be slightly different
% when you use up-to-date datasets.

load mitochondria 

%%
% Copy just the DNA sequence to a new variable |mitochondria|. You can
% access parts of the DNA sequence by using regular MATLAB indexing
% commands.
mitochondria = mitochondria_gbk.Sequence;
mitochondria_length = length(mitochondria)
first_300_bases = seqdisp(mitochondria(1:300))

%%
% You can look at the composition of the nucleotides with the |ntdensity|
% function. 

figure
ntdensity(mitochondria)

%%
% This shows that the mitochondria genome is A-T rich. The GC-content is
% sometimes used to classify organisms in taxonomy, it may vary between
% different species from ~30% up to ~70%. Measuring GC content is also
% useful for identifying genes and for estimating the annealing temperature
% of DNA sequence.

%% Calculating Sequence Statistics
% Now, you will use some of the sequence statistics functions in the
% Bioinformatics Toolbox(TM) to look at various properties of the human
% mitochondrial genome. You can count the number of bases of the whole
% sequence using the |basecount| function. 

bases = basecount(mitochondria)

%%
% These are on the 5'-3' strand. You can look at the reverse complement
% case using the |seqrcomplement| function.

compBases = basecount(seqrcomplement(mitochondria))

%%
% As expected, the base counts on the reverse complement strand are
% complementary to the counts on the 5'-3' strand.  
%%
% You can use the chart option to |basecount| to display a pie chart of the
% distribution of the bases.

figure
basecount(mitochondria,'chart','pie');
title('Distribution of Nucleotide Bases for Human Mitochondrial Genome');

%%
% Now look at the dimers in the sequence and display the information in a
% bar chart using |dimercount|.

figure
dimers = dimercount(mitochondria,'chart','bar')
title('Mitochondrial Genome Dimer Histogram');

%% Exploring the Open Reading Frames (ORFs)
% In a nucleotide sequence an obvious thing to look for is if there are any
% open reading frames. An ORF is any sequence of DNA or RNA that can be
% potentially translated into a protein. The function |seqshoworfs| can be
% used to visualize ORFs in a sequence. 

%%
% Note: In the HTML tutorial only the first page of the output is shown,
% however when running the example you will be able to inspect the complete
% mitochondrial genome using the scrollbar on the figure.

seqshoworfs(mitochondria);

%%
% If you compare this output to the genes shown on the NCBI page there seem
% to be slightly fewer ORFs, and hence fewer genes, than expected.

%%
% Vertebrate mitochondria do not use the Standard genetic code so some
% codons have different meaning in mitochondrial genomes. For more
% information about using different genetic codes in MATLAB see the help
% for the function |geneticcode|. The |GeneticCode| option to the
% |seqshoworfs| function allows you to look at the ORFs again but this time
% with the vertebrate mitochondrial genetic code.

%%
% In the human mitochondrial DNA sequence some genes are also started by
% alternative start codons [1]. Use the |AlternativeStartCodons| option to
% the |seqshoworfs| function to search also for these ORFs.

%%
% Notice that there are now two much larger ORFs on the third reading
% frame: One starting at position 4470 and the other starting at 5904.
% These correspond to the ND2 (NADH dehydrogenase subunit 2) and COX1
% (cytochrome c oxidase subunit I) genes. 

orfs = seqshoworfs(mitochondria,'GeneticCode','Vertebrate Mitochondrial',...
        'AlternativeStartCodons',true)


%% Inspecting Annotated Features
% You can also look at all the features that have been annotated to the
% human mitochondrial genome. Explore the complete GenBank entry
% |mitochondria_gbk| with the |featuresparse| function. Particularly, you
% can explore the annotated coding sequences (CDS) and compare them with
% the ORFs previously found. Use the |Sequence| option to the
% |featuresparse| function to extract, when possible, the DNA sequences
% respective to each feature. The |featuresparse| function will complement
% the pieces of the source sequence when appropriate.

features = featuresparse(mitochondria_gbk,'Sequence',true)
coding_sequences = features.CDS;
coding_sequences_id = sprintf('%s ',coding_sequences.gene)

%%
ND2CDS = coding_sequences(2) % ND2 is in the 2nd position 
COX1CDS = coding_sequences(3) % COX1 is in the 3rd position 

%% 
% Create a map indicating all the features found in this GenBank entry
% using the |featuresmap| function.
[h,l] = featuresmap(mitochondria_gbk,{'CDS','tRNA','rRNA','D_loop'},...
                                      [2 1 2 2 2],'Fontsize',9);
legend(h,l,'interpreter','none');
title('Homo sapiens mitochondrion, complete genome')
    
%% Extracting and Analyzing the ND2 and COX1 Proteins
% You can translate the DNA sequences that code for the ND2 and COX1 proteins
% by using the |nt2aa| function. Again the |GeneticCode| option must be
% used to specify the vertebrate mitochondrial genetic code.

ND2 = nt2aa(ND2CDS,'GeneticCode','Vertebrate Mitochondrial');
disp(seqdisp(ND2))

%%
COX1 = nt2aa(COX1CDS,'GeneticCode','Vertebrate Mitochondrial');
disp(seqdisp(COX1))

%%
% You can get a more complete picture of the amino acid content with
% |aacount|. 

figure
subplot(2,1,1)
ND2aaCount = aacount(ND2,'chart','bar');
title('Histogram of Amino Acid Count for the ND2 Protein');

subplot(2,1,2)
COX1aaCount = aacount(COX1,'chart','bar');
title('Histogram of Amino Acid Count for the COX1 Protein');

%%
% Notice the high leucine, threonine and isoleucine content and also the
% lack of cysteine or aspartic acid.

%%
% You can use the |atomiccomp| and |molweight| functions to calculate more
% properties about the ND2 protein.

ND2AtomicComp = atomiccomp(ND2)
ND2MolWeight = molweight(ND2)

%%
% For further investigation of the properties of the ND2 protein, use
% |proteinplot|. This is a graphical user interface (GUI) that allows you
% to easily create plots of various properties, such as hydrophobicity, of
% a protein sequence. Click on the "Edit" menu to create new properties,
% to modify existing property values, or, to adjust the smoothing
% parameters. Click on the "Help" menu in the GUI for more information on
% how to use the tool.   

proteinplot(ND2)

%%
% You can also programmatically create plots of various properties of the
% sequence using |proteinpropplot|.

figure
proteinpropplot(ND2,'PropertyTitle','Parallel beta strand')

%% Calculating the Codon Frequency using all the Genes in the Human Mitochondrial Genome

%%
% The |codoncount| function counts the number of occurrences of each codon
% in the sequence and displays a formatted table of the result. 

codoncount(ND2CDS) 

%%
% Notice that in the ND2 gene there are more CTA, ATC and ACC codons than
% others. You can check what amino acids these codons get translated into
% using the |nt2aa| and |aminolookup| functions.

CTA_aa = aminolookup('code',nt2aa('CTA')) 
ATC_aa = aminolookup('code',nt2aa('ATC'))
ACC_aa = aminolookup('code',nt2aa('ACC'))

%% 
% To calculate the codon frequency for all the genes you can concatenate
% them into a single sequence before using the function |codoncount|. You
% need to ensure that the codons are complete (three nucleotides each) so
% the read frame of the sequence is not lost at the concatenation.

numCDS = numel(coding_sequences);
CDS = cell(numCDS,1);
for i = 1:numCDS
     seq = coding_sequences(i).Sequence;
     CDS{i} = seq(1:3*floor(length(seq)/3));
end
allCDS = [CDS{:}];
codoncount(allCDS)

%%
% Use the |figure| option to the |codoncount| function to show a heat map
% with the codon frequency. Use the |geneticcode| option to overlay a grid
% on the figure that groups the  synonymous codons according with the
% Vertebrate Mitochondrial genetic code. Observe the particular bias of
% Leucine (codons 'CTN').

figure
count = codoncount(allCDS,'figure',true,'geneticcode','Vertebrate Mitochondrial');
title('Human Mitochondrial Genome Codon Frequency')

%% References
% [1] Barrell, B.G., Bankier A.T., and Drouin J. (1979). A different
% genetic code in human mitochondria. Nature 282, 189-194.

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20SEQSTATSDEMO%20in%20Bioinformatics%20Toolbox%204.2 Suggest an enhancement for this example.>*

displayEndOfDemoMessage(mfilename)
