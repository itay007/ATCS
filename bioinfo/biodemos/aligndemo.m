%% Aligning Pairs of Sequences
% This example shows how to extract some sequences from GenBank(R), find 
% open reading frames (ORFs), and then aligns the sequences using global 
% and local alignment algorithms. 

%  Copyright 2002-2013 The MathWorks, Inc.


%% Accessing NCBI Data from the MATLAB(R) Workspace
% One of the many fascinating sections of the NCBI web site is the
% <http://www.ncbi.nlm.nih.gov/bookshelf/br.fcgi?book=gnd Genes and 
% diseases section>. This section provides a comprehensive introduction
% to medical genetics.

%%
% In this example you will be looking at genes associated with
% <http://www.ncbi.nlm.nih.gov/bookshelf/br.fcgi?book=gnd&part=taysachsdisease
% Tay-Sachs Disease>. Tay-Sachs is an autosomal recessive disease caused by
% mutations in both alleles of a gene (HEXA, which codes for the alpha
% subunit of hexosaminidase A) on chromosome 15. 

%%
% The NCBI reference sequence for HEXA has accession number
% <http://www.ncbi.nlm.nih.gov/nuccore/189181665 NM_000520>.
% You can use the |getgenbank| function to retreive the sequence
% information from the NCBI data repository and load it into MATLAB(R). 
%
%   humanHEXA = getgenbank('NM_000520');

%%
% By doing a BLAST search or by searching in the mouse genome you can find
% an orthogonal gene, <http://www.ncbi.nlm.nih.gov/nuccore/26348756
% AK080777> is the accession number for a mouse hexosaminidase A gene.
%
%   mouseHEXA = getgenbank('AK080777');

%%
% For your convenience, previously downloaded sequences are included in a
% MAT-file. Note that data in public repositories is frequently curated and 
% updated; therefore the results of this example might be slightly different
% when you use up-to-date datasets.

load('hexosaminidase.mat','humanHEXA','mouseHEXA')

%% Exploring the Open Reading Frames (ORFs)
% You can use the function |seqshoworfs| to look for ORFs in the sequence
% for the human HEXA gene. Notice that the longest ORF is on the first
% reading frame. The output value in the variable |humanORFs| is a
% structure giving the position of the start and stop codons for all the
% ORFs on each reading frame.

humanORFs = seqshoworfs(humanHEXA.Sequence) 

%%
% Now look at the ORFs in the mouse HEXA gene. In this case the ORF is also
% on the first frame.

mouseORFs = seqshoworfs(mouseHEXA.Sequence)

%% Aligning the Sequences
% The first step is to use global sequence alignment to look for
% similarities between these sequences. 
% You could look at the alignment between the nucleotide sequences, but it is
% generally more instructive to look at the alignment between the protein
% sequences, in this example we know that the sequences are coding
% sequences. Use the |nt2aa| function to convert the nucleotide sequences
% into the corresponding amino acid sequences. Observe that the HEXA gene
% occurs in the first frame for both sequences, otherwise you should use
% the input argument |Frame| to specify an alternative coding frame.

humanProtein = nt2aa(humanHEXA.Sequence);
mouseProtein = nt2aa(mouseHEXA.Sequence);

%% 
% One of the easiest ways to look for similarity between sequences is with
% a dot plot. 

seqdotplot(mouseProtein,humanProtein)
xlabel('Human hexosaminidase A');ylabel('Mouse hexosaminidase A');

%%
% With the default settings, the dot plot is a little difficult to
% interpret, so you can try a slightly more stringent dot plot.

seqdotplot(mouseProtein,humanProtein,4,3)
xlabel('Human hexosaminidase A');ylabel('Mouse hexosaminidase A');

%% 
% The diagonal line indicates that there is probably a good alignment so 
% you can now take a look at the global alignment using the function
% |nwalign| which uses the Needleman-Wunsch algorithm.

[score, globalAlignment] = nwalign(humanProtein,mouseProtein);

%%
% The function |showalignment| displays the alignment in a MATLAB figure
% window with matching and similar residues highlighted in different colors. 

showalignment(globalAlignment);

%% Refining With Semi-global Alignment 
% The alignment is very good except for the terminal segments. For
% instance, notice the sparse matched pairs in the first positions. This
% occurs because a global alignment attempts to force the matching all the
% way to the ends and there is point where the penalty for opening new gaps
% is comparable to the score of matching residues. In some cases it is
% desirable to remove the gap penalty added at the ends of a global
% alignment; this allows you to better match this pair of sequences. This
% technique is commonly known as 'semi-global' alignment or 'glocal'
% alignment.

[score, globalAlignment] = nwalign(humanProtein,mouseProtein,'glocal',true);
showalignment(globalAlignment);

%% Refining the Alignment by Extracting the Protein Sequence
% Another way to refine your alignment is by using only the protein
% sequences. Notice that the aligned region is delimited by start
% ( M-methionine ) and stop ( * ) amino acids in the sequences. If the
% sequence is shortened so that only the translated regions are considered,
% then it seems likely that you will get a better alignment. Use the |find|
% command to look for the index of the start amino acid in each sequence:

humanStart = find(humanProtein == 'M',1)
mouseStart = find(mouseProtein == 'M',1)

%%
% Similarly, use the |find| command to look for the index of the first
% stop occurring after the start of the translation. Special care needs to
% be taken because there is also a stop at the very beginning of the
% |humanProtein| sequence.

humanStop = find(humanProtein(humanStart:end)=='*',1) + humanStart - 1
mouseStop = find(mouseProtein(mouseStart:end)=='*',1) + mouseStart - 1

%%
% Use these indices to truncate the sequences.

humanSeq = humanProtein(humanStart:humanStop);
humanSeqFormatted = seqdisp(humanSeq)
mouseSeq = mouseProtein(mouseStart:mouseStop);
mouseSeqFormatted = seqdisp(mouseSeq)


%%
% If you align these two sequences and then view them you will see a very
% good global alignment.

[score, alignment] = nwalign(humanSeq,mouseSeq);
showalignment(alignment);

%%
% Open reading frame information is also available from the output of the
% |seqshoworfs| command, but the indices are based on the nucleotide
% sequences. Use these indices to trim the original nucleotide sequences
% and then translate them to amino acids.

humanPORF = nt2aa(humanHEXA.Sequence(humanORFs(1).Start(1):humanORFs(1).Stop(1)));
mousePORF = nt2aa(mouseHEXA.Sequence(mouseORFs(1).Start(1):mouseORFs(1).Stop(1)));
[score, ORFAlignment] = nwalign(humanPORF,mousePORF);
showalignment(ORFAlignment);

%% 
% Alternatively, you can use the coding region information (CDS) from the
% GenBank data structure to find the coding region of the genes. 

idx = humanHEXA.CDS.indices;
humanCodingRegion = humanHEXA.Sequence(idx(1):idx(2));
idx = mouseHEXA.CDS.indices;        
mouseCodingRegion = mouseHEXA.Sequence(idx(1):idx(2));

%% 
% You can also get the translation of the coding regions from this structure.

humanTranslatedRegion = humanHEXA.CDS.translation;
mouseTranslatedRegion = mouseHEXA.CDS.translation;

%% Local Alignment
% Instead of truncating the sequences to look for better alignment, an
% alternative approach is to use a local alignment. The function |swalign|
% performs local alignment using the Smith-Waterman algorithm. This shows a
% very good alignment for the whole coding region and reasonable similarity
% for a few residues beyond at both the ends of the gene. 

[score, localAlignment] = swalign(humanProtein,mouseProtein);
showalignment(localAlignment);

%% Alignment of Complementary DNA Sequences
% All the sequence alignment functions provided in MATLAB can be
% customized. For example, by modifying the rows and columns of a scoring
% matrix you can align sequences by complement and not by identity.
% In this case you can reorder the |NUC44| scoring matrix; a positive score
% is given for complements while a negative score is given otherwise. The
% first 30 nucleotides from the mouse HEXA gene will be aligned to its
% complement.

[M, info] = nuc44;
map = nt2int(seqcomplement(info.Order))
Mc = M(:,map)
[score, compAlignment] = nwalign(mouseHEXA.Sequence(1:30), ...
    seqcomplement(mouseHEXA.Sequence(1:30)), 'SCORINGMATRIX', ...
        Mc, 'ALPHABET', 'NT')
    
%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20ALIGNDEMO%20in%20Bioinformatics%20Toolbox%204.2 Suggest an enhancement for this example.>*

displayEndOfDemoMessage(mfilename)
    
