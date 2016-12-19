%% Genetic Algorithm Search for Features in Mass Spectrometry Data
% This example shows how to use the *Global Optimization Toolbox* with the
% *Bioinformatics Toolbox(TM)* to optimize the search for features to
% classify mass spectrometry (SELDI) data.

%   Copyright 2005-2012 The MathWorks, Inc.


%% Introduction
% Genetic algorithms optimize search results for problems with large data
% sets. You can use the MATLAB(R) genetic algorithm function to solve these 
% problems in Bioinformatics. Genetic algorithms have been applied to
% phylogenetic tree building, gene expression and mass spectrometry data
% analysis, and many other areas of Bioinformatics that have large and
% computationally expensive problems. This example searches for optimal
% features (peaks) in mass spectrometry data. We will look for specific
% peaks in the data that distinguish cancer patients from control patients.

%% Global Optimization Toolbox
% First familiarize yourself with the Global Optimization Toolbox. The
% documentation describes how a genetic algorithm works and how to use it
% in MATLAB. To access the documentation, use the *doc* command.

doc ga

%% Load Mass Spectrometry Data into MATLAB(R)
% The data in this example is also used in the example
% <cancerdetectdemo.html Identifying Significant Features and Classifying
% Protein Profiles> and is available from the FDA-NCI Clinical Proteomics
% Program Databank at
% http://home.ccr.cancer.gov/ncifdaproteomics/ppatterns.asp. It is a
% collection of samples from 121 ovarian cancer patients and 95 control  
% patients. This example assumes that you downloaded, uncompressed, and
% preprocessed the raw mass-spectrometry data from the FDA-NCI Web site.
% You can recreate the preprocessed data file needed for this example, 
% OvarianCancerQAQCdataset.mat, by running the script *msseqprocessing*, or
% by following the steps in the example <biodistcompdemo.html Batch
% Processing of Spectra Using Sequential and Parallel Computing>. There
% should be two variables  that are loaded into MATLAB (*MZ* and *Y*). *MZ*
% is the mass/charge vector while *Y* is the intensity values for all 216
% patients (control and cancer). To visualize this data see the example
% <cancerdetectdemo.html Identifying Significant Features and Classifying
% Protein Profiles>.  

load OvarianCancerQAQCdataset
whos

%%
% Initialize the variables used in the example.

[numPoints numSamples] = size(Y); % total number of samples and data points
id = grp2idx(grp);  % ground truth: Cancer=1, Control=2

%% Create a Fitness Function for the Genetic Algorithm
% A genetic algorithm requires an objective function, also known as the
% fitness function, which describes the phenomenon that we want to
% optimize. In this example, the genetic algorithm machinery tests small
% subsets of  M/Z values using the fitness function and then determines
% which M/Z values get passed on to or removed from each subsequent
% generation. The fitness function *biogafit* is passed to the genetic
% algorithm solver using a function handle. In this example, *biogafit*
% maximizes the separability of two classes by using a linear combination
% of 1) the a-posteriori probability and 2) the empirical error rate of a
% linear classifier (*classify*). You can create your own fitness function
% to try different classifiers or alternative methods for assessing the
% performance of the classifiers.    

type biogafit.m

%% Create an Initial Population
% Users can change how the optimization is performed by the genetic
% algorithm by creating custom functions for crossover, fitness scaling,
% mutation, selection, and population creation. In this example you will
% use the *biogacreate* function written for this example to create initial
% random data points from the mass spectrometry data. The function header
% requires specific input parameters as specified by the GA documentation.
% There is a default creation function in the toolbox for creating initial
% populations of data points.

type biogacreate.m

%% Set Genetic Algorithm Options
% The GA function uses an options structure to hold the algorithm
% parameters that it uses when performing a minimization with a genetic
% algorithm. The *gaoptimset* function will create this options structure.
% For the purposes of this example, the genetic algorithm will run only for
% 50 generations. However, you may set 'Generations' to a larger value.  

options = gaoptimset('CreationFcn', {@biogacreate,Y,id},...
                     'PopulationSize',100,...
                     'Generations',50,...
                     'Display', 'iter')
        
%% Run GA to Find 20 Discriminative Features
% Use *ga* to start the genetic algorithm function. 100 groups of 20
% datapoints each will evolve over 50 generations. Selection, crossover,
% and mutation events generate a new population in every generation. 

% initialize the random generators to the same state used to generate the
% published example 
rand('seed',1)      
randn('seed',1)     
nVars = 20;                          % set the number of desired features   
FitnessFcn = {@biogafit,Y,id};       % set the fitness function   
feat = ga(FitnessFcn,nVars,options); % call the Genetic Algorithm   
    
feat = round(feat);   
Significant_Masses = MZ(feat)   
    
cp = classperf(classify(Y(feat,:)',Y(feat,:)',id),id);   
cp.CorrectRate   
    
%% Display the Features that are Discriminatory   
% To visualize which features have been selected by the genetic algorithm,   
% the data is plotted with peak positions marked with red vertical lines.   
    
xAxisLabel = 'Mass/Charge (M/Z)';       % x label for plots   
yAxisLabel = 'Relative Ion Intensity';  % y label for plots   
figure; hold on;   
hC = plot(MZ,Y(:,1:15) ,'b');   
hN = plot(MZ,Y(:,141:155),'g');   
hG = plot(MZ(feat(ceil((1:60 )/3))), repmat([0 100 NaN],1,20),'r');   
xlabel(xAxisLabel); ylabel(yAxisLabel);   
axis([1900 12000 -1 40]);   
legend([hN(1),hC(1),hG(1)],{'Control','Ovarian Cancer', 'Discriminative Features'},2);   
title('Mass Spectrometry Data with Discriminative Features found by Genetic Algorithm');   
    
%%   
% Observe the interesting peak around 8100 Da., which seems to be shifted   
% to the right on healthy samples.   
axis([8082 8128 -.5 4]) 

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20MSGADEMO%20in%20Bioinformatics%20Toolbox%204.2 Provide feedback for this example.>*

displayEndOfDemoMessage(mfilename)
