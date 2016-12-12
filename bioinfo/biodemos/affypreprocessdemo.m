%% Preprocessing Affymetrix(R) Microarray Data at the Probe Level
% This example shows how to use MATLAB(R) and Bioinformatics Toolbox(TM)
% for preprocessing Affymetrix(R) oligonucleotide microarray probe-level
% data with two preprocessing techniques, Robust Multi-array Average (RMA)
% and GC Robust Multi-array Average (GCRMA).
 
%   Copyright 2006-2012 The MathWorks, Inc.

%% Introduction
% With Affymetrix oligonucleotide microarray platforms, gene expression is
% measured using probe sets consisting of 11 to 20 perfect match (PM)
% probes (25 nucleotides in length) complementary to target mRNA sequences.
% Each probe set also has the same number of mismatch (MM) probes, in which
% the 13th nucleotide has been changed to its complement. The PM probes are
% designed for gene specific hybridization. The control MM probe
% measurements are thought to comprise most of the background non-specific
% binding, such as cross-hybridization. A PM probe and its corresponding MM
% probe are referred to as a probe pair.

%%
% The measured probe intensities and locations from a hybridized microarray
% are stored in a CEL file. For each Affymetrix microarray platform, the
% information relating probe pairs to probe set IDs, and to locations on
% the array is stored in a CDF library file. The probe sequence information
% is stored in a sequence file (FASTA or tab-separated format).

%%
% In general, preprocessing Affymetrix probe-level expression data consists
% of three steps: background adjustment, normalization, and summarization
% at the probe set level as a measure of the expression level of
% corresponding mRNA. Many methods exist for the statistical procedures of
% these three steps. Two popular techniques, RMA (Irizarry et al., 2003)
% and GCRMA (Wu et al., 2004), are used in this example.

%%
% *Note:* This example shows the RMA and GCRMA preprocessing procedures to
% compute expression values from input CEL files in step-by-step detail,
% using several functions. You can also complete the same RMA or GCRMA
% techniques in one function call by using the Bioinformatics Toolbox
% |affyrma| or |affygcrma| functions, respectively.

%%
% A publicly available dataset containing Affymetrix microarray
% measurements of 42 tumor tissues of the embryonal central nervous system
% (CNS, Pomeroy et al., 2002) is used for this example. You will import and
% access the probe level data of multiple arrays, and then perform
% expression level measurements with RMA and GCRMA preprocessing methods.

%% Importing Data
% The CNS experiment was conducted using the Affymetrix HuGeneFL GeneChip(R)
% array, and the data were stored in CEL files. Information related to each
% probe is contained in the Affymetrix Hu6800 CDF library file. 

%%
% If you dont already have the Hu6800 CDF library file, download the
% <http://www.affymetrix.com/Auth/support/downloads/library_files/hugenefl_libraryfile.zip
% HuGeneFL Genome Array library zip file>. Extract the Hu6800.CDF file into
% a directory, such as |C:\Examples\affypreprocessdemo\libfiles|. Note: You
% will have to register in order to access the library files, but you do
% not have to run the setup.exe file in the archive.

%%
% The CNS dataset (CEL files) is available at the
% <http://www.broad.mit.edu/mpr/CNS/ CNS experiment> Website. To complete
% this example, download the CEL files of the CNS dataset into a directory,
% such as |C:\Examples\affypreprocessdemo\data|. Unzip the CEL file
% archives. Note: This dataset contains more CEL files than are needed for
% this example.

%%
% CNS_DataA_Sample_CEL.txt, a file provided with Bioinformatics Toolbox,
% contains a list of the 42 CEL filenames used for this example, and the
% samples (10 medulloblastomas, 10 rhabdoid, 10 malignant glioma, 8
% supratentorial PNETS, and 4 normal human cerebella) to which they belong.
% Load this data into two MATLAB variables.

fid = fopen('CNS_DataA_Sample_CEL.txt','r');
ftext = textscan(fid,'%q%q');
fclose(fid);
samples = ftext{1};
cels = ftext{2};

%%
% The function |celintensityread| can read multiple CEL files and access a
% CDF library file. It returns a MATLAB structure containing the probe
% information and probe intensities. The matrices of PM and MM intensities
% from multiple CEL files are stored in the |PMIntensities| and
% |MMIntensities| fields. In each probe intensity matrix, the column
% indices correspond to the order in which the CEL files were read, and
% each row corresponds to a probe. Create a MATLAB structure of PM and MM
% probe intensities by loading data from the CEL files from the directory
% where the CEL files are stored, and pass in the path to where you stored
% the CDF library file. (Note: |celintensityread| will report the progress
% to the MATLAB command window. You can turn the progress report off by
% setting the input parameter |VERBOSE| to false.) 
celPath = 'C:\Examples\affypreprocessdemo\data';
libPath = 'C:\Examples\affypreprocessdemo\libfiles'; 
probeData = celintensityread(cels, 'Hu6800.CDF',...
                 'celpath', celPath, 'cdfpath', libPath, 'pmonly', false)


%%
% Determine the number of CEL files loaded.
nSamples = probeData.NumChips

%%
% Determine the number of probe sets on a HuGeneFL array.
nProbeSets = probeData.NumProbeSets

%% 
% Determine the number of probes on a HuGeneFL array.
nProbes = probeData.NumProbes


%%
% To perform GCRMA preprocessing, the probe sequence information of the
% HuGeneFL array is also required. The Affymetrix
% <http://www.affymetrix.com/support/index.affx support site> provides
% probe sequence information for most of the available arrays, either as
% FASTA formatted or tab-delimited files. This example assumes you have the
% HuGeneFL_probe_tab file in the |C:\Examples\affypreprocessdemo\libfiles|
% directory. Use the function |affyprobeseqread| to parse the sequence file
% and return the probe sequences in an |nProbes| x 25 matrix of integers
% that represents the PM probe sequence bases, with rows corresponding to
% the probes on the chip and columns corresponding to the base positions of
% the 25-mer.
seqPath = 'C:\Examples\affypreprocessdemo\libfiles'; 
S = affyprobeseqread('HuGeneFL_probe_tab', 'Hu6800.CDF',...
                'seqpath', seqPath, 'cdfpath', libPath, 'seqonly', true)
            
            
%% Preprocessing Probe-Level Expression Data 
% The RMA procedure uses only PM probe intensities for background
% adjustment (Irizarry et al., 2003), while GCRMA adjusts background using
% probe sequence information and MM control probe intensities to estimate
% non-specific binding (Wu et al., 2004).  Both RMA and GCRMA are preceded
% by quantile normalization (Bolstad et al., 2003) and median polish
% summarization (Irizarry et al., 2003) of PM intensities. 

%% Using the RMA Procedure
% The RMA background adjustment method corrects PM probe intensities chip
% by chip. The PM probe intensities are modeled as the sum of a normal
% noise component and an exponential signal component. Use |rmabackadj| to
% background adjust the PM intensities in the CNS data. You can inspect the
% intensity distribution histogram and the estimated background adjustment
% of a specific chip by setting the input parameter |SHOWPLOT| to the
% column index of the chip.
pms_bg = rmabackadj(probeData.PMIntensities, 'showplot', 1);

%%
% Several nonlinear normalization methods have been successfully applied to
% Affymetrix microarray data. The RMA procedure normalizes the probe-level
% data with a quantile normalization method. Use |quantilenorm| to
% normalize the background adjusted PM intensities in the CNS data. Note:
% If you are interested in a rank-invariant set normalization method, use
% the |affyinvarsetnorm| function instead. 
pms_bgnorm = quantilenorm(pms_bg);

%%
% A median polish procedure is applied to the PM intensities in
% summarization. To calculate the expression values, use |rmasummary| to
% summarize probe intensities of each probe set across multiple chips. The
% expression values are the probe set intensity summaries on a log-2 scale.
cns_rma_exp = rmasummary(probeData.ProbeIndices, pms_bgnorm);

%% Using the GCRMA Procedure
% The GCRMA procedure adjusts for optical noise and non-specific binding
% (NSB) taking into account the effect of the stronger bonding of G/C pairs
% (Naef et al., 2003, Wu et al., 2004). GCRMA uses probe sequence
% information to estimate probe affinities for computing non-specific
% binding. The probe affinity is modeled as a sum of the position-dependent
% base effects. Usually, the probe affinities are estimated from the MM
% intensities of an NSB experiment. If NSB data is not available, the probe
% affinities can still be estimated from sequence information and MM probe
% intensities normalized by the probe set median intensity (Naef et al.,
% 2003). 

%%
% For the CNS dataset, use the data from the microarray hybridized with the
% normal cerebella sample (Brain_Ncer_1) to compute the probe affinities
% for the HuGeneFL array. Use |affyprobeaffinities| to estimate the probe
% affinities of an Affymetrix microarray. Use the |SHOWPLOT| input
% parameter to inspect a plot showing the effects of base A, C, G, and T at
% the 25 positions. 
figure
idx = find(strcmpi('Brain_Ncer_1', samples));
[pmAlpha, mmAlpha] = affyprobeaffinities(S.SequenceMatrix,...
                       probeData.MMIntensities(:, idx), 'showplot', true);

%%
% Note: There are 496 probes on a HuGeneFL array that do not have sequence
% information; the affinities for these probes were NaN.

%% 
% With the probe affinities available, the amount of NSB can be estimated
% by fitting a LOWESS curve through MM probe intensities vs. MM probe
% affinities. The function |gcrmabackadj| performs optical and NSB
% corrections. The input parameter |SHOWPLOT| shows a plot of the optical
% noise adjusted MM intensities against its affinities, and the smooth fit 
% of a specified chip. You can compute the background intensities with one
% of two estimation methods, Maximum Likelihood Estimate (MLE) and
% Empirical-Bayes (EB), which computes the posterior mean of specific
% binding given prior observed intensities. Here you will background adjust
% four arrays using both estimation methods. (Note: |gcrmabackadj| will
% report the progress to the MATLAB command window. You can turn the
% progress report off by setting the input parameter |VERBOSE| to false.) 

%%
% Background adjust the first four chips using GCRMA-MLE method, and
% inspect the plot of intensity vs. affinity for data from the third array.
pms_MLE_bg = gcrmabackadj(probeData.PMIntensities(:,1:4),...
                            probeData.MMIntensities(:, 1:4),...
                            pmAlpha, mmAlpha, 'showplot', 3);
                      
% Adjust YLIM for better view
ylim([4 16]);

%%
% Background adjust the first four chips using the GCRMA-EB method.
% Processing with this method is more computationally intensive and will
% take longer.
pms_EB_bg = gcrmabackadj(probeData.PMIntensities(:,1:4),...
                            probeData.MMIntensities(:, 1:4),...
                            pmAlpha, mmAlpha, 'method', 'EB');
                      
%%
% You can continue the preprocessing with the |quatilenorm| and
% |rmasummary| functions, or use the |gcrma| function to do everything. The
% |gcrma| function performs background adjustment and returns expression
% measures of background adjusted PM probe intensities using the same
% normalization and summarization methods as RMA. You can also pass in the
% sequence matrix instead of affinities. The function will automatically
% compute the affinities in this case. (Note: |gcrma| will report the
% progress to the MATLAB command window. You can turn the progress report
% off by setting the input parameter |VERBOSE| to false.) 
cns_mle_exp = gcrma(probeData.PMIntensities, probeData.MMIntensities,...
                    probeData.ProbeIndices, pmAlpha, mmAlpha);

%% Inspecting the Background Adjustment Results
% Use boxplots to inspect the PM intensity distributions of the first four
% chips  with three background adjustment procedures. 
figure
subplot(4,1,1)
maboxplot(log2(probeData.PMIntensities(:, 1:4)), samples(1:4),...
          'title','Raw Intensities', 'orientation', 'horizontal')
subplot(4,1,2)
maboxplot(log2(pms_bg(:,1:4)), samples(1:4),...
          'title','After RMA Background Adjustment','orient','horizontal')  
subplot(4,1,3)
maboxplot(log2(pms_MLE_bg), samples(1:4),...
          'title','After GCRMA-MLE Background Adjustment','orient','horizontal')  
subplot(4,1,4)
maboxplot(log2(pms_EB_bg), samples(1:4),...
          'title','After GCRMA-EB Background Adjustment','orient','horizontal')  

%%
% Use boxplots to inspect the background corrected and normalized PM
% intensity distributions of the first four chips  with three background
% adjustment procedures. 
pms_MLE_bgnorm = quantilenorm(pms_MLE_bg);
pms_EB_bgnorm  = quantilenorm(pms_EB_bg);

figure
subplot(3,1,1)
maboxplot(log2(pms_bgnorm(:, 1:4)), samples(1:4),...
          'title','Normalized RMA Background Adjusted Intensity',...
          'orientation', 'horizontal')
subplot(3,1,2)
maboxplot(log2(pms_MLE_bgnorm), samples(1:4),...
          'title','Normalized GCRMA-MLE Background Adjusted Intensity',...
          'orientation', 'horizontal')
subplot(3,1,3)
maboxplot(log2(pms_EB_bgnorm), samples(1:4),...
          'title','Normalized GCRMA-EB Background Adjusted Intensity',...
          'orientation', 'horizontal')

%% Final Remarks
% You can perform importing of data from CEL files and all three
% preprocessing steps of the RMA and GCRMA techniques shown in this example
% by using the |affyrma| and |affygcrma| functions respectively.

%% 
% For more information on gene expression analysis with Bioinformatics
% Toolbox, see <cnsgeneexpdemo.html Exploring Gene Expression Data>.

%%
% Affymetrix and GeneChip are registered trademarks of Affymetrix, Inc.

%% References
% [1] Pomeroy, S.L., Tamayo, P., Gaasenbeek, M., Sturla, L.M., Angelo, M.,
%     McLaughlin, M.E., Kim, J.Y., Goumnerova, L.C., Black, P.M., Lau, C.,
%     Allen, J.C., Zagzag, D., Olson, J.M., Curran, T., Wetmore, C.,
%     Biegel, J.A., Poggio, T., Mukherjee, S., Rifkin, R., Califano, A.,
%     Stolovitzky, G., Louis, DN, Mesirov, J.P., Lander, E.S., and Golub,
%     T.R. (2002). Prediction of central nervous system embryonal tumour
%     outcome based on gene expression. Nature, 415(6870) 436-442.
%%
% [2] Irizarry, R.A., Hobbs, B., Collin, F., Beazer-Barclay, Y.D., 
%     Antonellis, K.J., Scherf, U., and Speed. T.P. (2003). Exploration,
%     normalization, and summaries of high density oligonucleotide array
%     probe level data. Biostatistics, 4, 249-264.
%%
% [3] Wu, Z., Irizarry, R.A., Gentleman, R., Murillo, F.M., and Spencer, F.
%     (2004). A model based background adjustment for oligonucleotide
%     expression arrays. J. Amer. Stat. Assoc. 99, 909-917.
%%
% [4] Bolstad, B.M., Irizarry R.A., Astrand, M., and Speed, T.P. (2003). A
%     comparison of normalization methods for high density oligonucleotide
%     array data based on variance and bias. Bioinformatics, 19, 185-193.
%%
% [5] Naef, F., and Magnasco, M.O. (2003). Solving the riddle of the bright
%     mismatches: labeling and effective binding in oligonucleotide arrays.
%     Phys. Rev. E 68, 011906. 

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20AFFYPREPROCESSDEMO%20in%20Bioinformatics%20Toolbox%204.2 Suggest an enhancement for this example.>*

displayEndOfDemoMessage(mfilename)

