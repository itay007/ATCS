%% Calling Bioperl Functions from MATLAB(R)
% This example shows the interoperability between MATLAB(R) and Bioperl -
% passing arguments from MATLAB to Perl scripts and pulling BLAST search
% data back to MATLAB.
%
% NOTE: Perl and the Bioperl modules must be installed to run the Perl
% scripts in this example. Since version 1.4, Bioperl modules have a
% warnings.pm dependency requiring at least version 5.6 of Perl. If you
% have difficulty running the Perl scripts, make sure your PERL5LIB
% environment variable includes the path to your Bioperl installation or
% try running from the Bioperl installation directory. See the links at
% http://www.perl.com and http://bioperl.org/ for current release files and
% complete installation instructions.

%   Copyright 2003-2012 The MathWorks, Inc.


%% Introduction
% Gleevec(TM) (STI571 or imatinib mesylate) was the first approved drug to
% specifically turn off the signal of a known cancer-causing protein.
% Initially approved to treat chronic myelogenous leukemia (CML), it is
% also effective for treatment of gastrointestinal stromal tumors (GIST).
%
% If you have access to the Internet, run this command to learn more:
% web('http://www.cancer.gov/clinicaltrials/digestpage/gleevec')

%%
% Research has identified several gene targets for Gleevec including:
% Proto-oncogene tyrosine-protein kinase ABL1 (NP_009297), Proto-oncogene
% tyrosine-protein kinase Kit (NP_000213), and Platelet-derived growth
% factor receptor alpha precursor (NP_006197).

target_ABL1 = 'NP_009297';
target_Kit = 'NP_000213';
target_PDGFRA = 'NP_006197';

%% Accessing Sequence Information
% You can load the sequence information for these proteins from local
% GenPept text files using *genpeptread*.

ABL1_seq = getfield(genpeptread('ABL1_gp.txt'), 'Sequence');
Kit_seq = getfield(genpeptread('Kit_gp.txt'), 'Sequence');
PDGFRA_seq = getfield(genpeptread('PDGFRA_gp.txt'), 'Sequence');

%%
% Alternatively, you can obtain protein information directly from the
% online GenPept database maintained by the National Center for
% Biotechnology Information (NCBI).
%
% Run these commands to download data from NCBI:

% ABL1_seq = getgenpept(target_ABL1, 'SequenceOnly', true);
% Kit_seq = getgenpept(target_Kit, 'SequenceOnly', true);
% PDGFRA_seq = getgenpept(target_PDGFRA, 'SequenceOnly', true);

%%
% The MATLAB *whos* command gives information about the size of these
% sequences.

whos ABL1_seq
whos Kit_seq
whos PDGFRA_seq

%% Calling Perl Programs from MATLAB
% From MATLAB, you can harness existing Bioperl modules to run a BLAST
% search on these sequences. MW_BLAST.pl is a Perl program based on the
% RemoteBlast Bioperl module. It reads sequences from FASTA files, so
% start by creating a FASTA file for each sequence.

fastawrite('ABL1.fa', 'ABL1 Proto-oncogene tyrosine-protein kinase (NP_009297)', ABL1_seq);
fastawrite('Kit.fa', 'Kit Proto-oncogene tyrosine-protein kinase (NP_000213)', Kit_seq);
fastawrite('PDGFRA.fa', 'PDGFRA alpha precursor (NP_006197)', PDGFRA_seq);

%%
% BLAST searches can take a long time to return results, and the Perl
% program MW_BLAST includes a repeating sleep state to await the report.
% Sample results have been included with this example, but if you have an
% Internet connection and want to try running the BLAST search with the
% three sequences, uncomment the following commands. MW_BLAST.pl will save
% the BLAST results in three files on your disk, ABL1.out, Kit.out and
% PDGFRA.out. The process can take 15 minutes or more.

% try
%     perl('MW_BLAST.pl','blastp','pdb','1e-10','ABL1.fa','Kit.fa','PDGFRA.fa');
% catch
%     error(message('bioinfo:bioperldemo:PerlError'))
% end

%%
% Here is the Perl code for MW_BLAST:

type MW_BLAST.pl

%%
% The next step is to parse the output reports and find scores >= 100. You
% can then identify hits found by more than one protein for further
% research, possibly identifying new targets for drug therapy.

try
    protein_list = perl('MW_parse.pl', 'ABL1.out', 'Kit.out', 'PDGFRA.out')
catch
    error(message('bioinfo:bioperldemo:PerlError'))
end

%%
% This is the code for MW_parse:

type MW_parse.pl

%% Calling MATLAB Functions within Perl Programs
% If you are running on Windows(R), it is also possible to call MATLAB
% functions from Perl. You can launch MATLAB in an Automation Server mode
% by using the /Automation switch in the MATLAB startup command
% (e.g. D:\applications\matlab7x\bin\matlab.exe /Automation).
%
% Here's a script to illustrate the process of launching an automation
% server, calling MATLAB functions and passing variables between Perl and
% MATLAB.

type MATLAB_from_Perl.pl

%% Protein Analysis Tools in the Bioinformatics Toolbox(TM)
% MATLAB offers additional tools for protein analysis and further research
% with these proteins. For example, to access the sequences and run a full
% Smith-Waterman alignment on the tyrosine kinase domain of the human
% insulin receptor (pdb 1IRK) and the kinase domain of the human lymphocyte
% kinase (pdb 3LCK), load the sequence data:

IRK = pdbread('pdb1irk.ent');
LCK = pdbread('pdb3lck.ent');

% Run these commands to bring the data from the Internet:
% IRK = getpdb('1IRK');
% LCK = getpdb('3LCK');

%%
% Now perform a local alignment with the Smith-Waterman algorithm. MATLAB
% uses BLOSUM 50 as the default scoring matrix for AA strings with a gap
% penalty of 8. Of course, you can change any of these parameters.

[Score, Alignment] = swalign(IRK, LCK, 'showscore', true);
%%
showalignment(Alignment);

%%
% MATLAB and the Bioinformatics Toolbox(TM) offer additional tools for
% investigating nucleotide and amino acid sequences. For example,
% *pdbdistplot* displays the distances between atoms and amino acids in a
% PDB structure, while *ramachandran* generates a plot of the torsion angle
% PHI and the torsion angle PSI of the protein sequence. The toolbox
% function *proteinplot* provides a graphical user interface (GUI) to
% easily import sequences and plot various properties such as
% hydrophobicity.

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20BIOPERLDEMO%20in%20Bioinformatics%20Toolbox%204.2 Suggest an enhancement for this example.>*

displayEndOfDemoMessage(mfilename)
