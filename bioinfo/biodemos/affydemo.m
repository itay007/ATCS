%% Working with Affymetrix(R) Data
% This example shows how to use the functions in the Bioinformatics
% Toolbox(TM) for working with Affymetrix(R) GeneChip(R) data.

% Copyright 2004-2012 The MathWorks, Inc.

%% About Affymetrix Data Files
% The function |affyread| can read four types of Affymetrix data files.
% These are DAT files, which contain raw image data, CEL files which
% contain information about the intensity values of the individual probes,
% CHP files which contain information about probe sets, and EXP files,
% which contain information about experimental conditions and protocols.
% |affyread| can also read CDF and GIN library files. The CDF file contains
% information about which probes belong to which probe set and the GIN file
% contains information about the probe sets such as the gene name with
% which the probe set is associated. To learn more about the actual files,
% you can download sample data files from the
% <http://www.affymetrix.com/Auth/support/technical/sample_data/demo_data.affx
% Affymetrix Support Site>. Most of the data sets are stored in DTT
% archives. To extract the DAT, CEL and CHP files you will need to install
% the <http://www.affymetrix.com/support/technical/software_downloads.affx
% Data Transfer Tool>.

%% Downloading the E. coli Antisense Data Set
% For this example, you will need some sample data files (DAT, CEL, CHP)
% from the _E. coli_ Antisense Genome Array. Download these from
% <http://www.affymetrix.com/Auth/support/downloads/demo_data/Demo_Data_E-coli-antisense.zip
% Demo_Data_E-coli-antisense.zip> Extract the data files from the DTT
% archive using the
% <http://www.affymetrix.com/support/technical/software_downloads.affx Data
% Transfer Tool>. Modify this line to the name of the path and directory to
% which you extracted the sample data files.

exampleDataDir = 'C:\Examples\affydemo\data';

%% Downloading E. coli Antisense Library Files
% In addition to the data files, you will also need Ecoli_ASv2.CDF and
% Ecoli_ASv2.GIN, the library files for the _E. coli_ Antisense Genome
% Array. You may already have these files if you have any Affymetrix
% GeneChip software installed on your machine. If not, get the library
% files by downloading and unzipping the
% <http://www.affymetrix.com/Auth/support/downloads/library_files/ecoli_antisense_libraryfile.zip 
% _E. coli_ Antisense Genome Array zip file>
%
% Note that you will have to register in order to access the library files.
%
% You only have to unzip the files, you do not have to run the Setup.exe
% file in the archive.  
%
% Modify this line to the name of the path and directory to which you
% extracted the library files.

libDir = 'C:\Examples\affydemo\libfiles';

%% Image Files (DAT Files)
% The raw image data from the chip scanner is saved in the DAT file. If you
% use |affyread| to read a DAT file you will see that it creates a MATLAB(R)
% structure.

datStruct = affyread(fullfile(exampleDataDir,'Ecoli-antisense-121502.dat'))


%%
% You can access fields of the structure using the dot notation.

datStruct.NumRows

%% Displaying an Image File
% You can use the |imagesc| command to display the image.
datFigure = figure;
imagesc(datStruct.Image);

%%
% You can change the |colormap| from the default |jet| to another using the
% |colormap| command.

colormap pink

%%
% You can zoom in on a particular area by using the Zoom In tool with the
% mouse, or by using the |axis| command. Notice that this stretches the
% y-axis.

axis([1900 2800 160 650])

%% 
% You can use the |axis image| command to set the correct aspect ratio.

axis image
axis([1900 2800 160 650])


%% Probe Results Files (CEL Files)
% The information about each probe on the chip is extracted from the image
% data by the Affymetrix image analysis software. The information is stored
% in the CEL file. |affyread| reads a CEL file into a structure. Notice
% that many of the fields are the same as those in the DAT structure.

celStruct = affyread(fullfile(exampleDataDir,'Ecoli-antisense-121502.CEL'))

%%
% The CEL file contains information about where each probe is on the chip
% and also the intensity values for the probe. You can use the |maimage|
% function to display the chip.

celFigure = figure;
maimage(celStruct)

%%
% Again, you can zoom in on a specific region.  

axis([200 340 0 70])

%%
% If you compare the image created from the CEL file and the image created
% from the DAT file, you will notice that the CEL image is lower
% resolution. This is because there is only one pixel per probe in this
% image, whereas the DAT file image has many pixels per probe.

%%
% The structures created by |affyread| can be very large. It is a good idea
% to clear them from memory once they are no longer needed.

clear datStruct
close(datFigure); close(celFigure);

%%
% The |Probes| field of the CEL structure contains information about the
% individual probes. There are eight values per probe. These are stored in
% the |ProbeColumnNames| field of the structure.

celStruct.ProbeColumnNames

%%
% So if you look at one row of the |Probes| field of the CEL structure you
% will see eight values corresponding to the X position, Y position,
% intensity, and so forth.

celStruct.Probes(1:10,:)

%% Results Files (CHP Files)
% The CHP file contains the results of the experiment. These include the
% average signal measures for each probe set as determined by the
% Affymetrix software and information about which probe sets are called as
% present, absent or marginal and the p-values for these calls.

chpStruct = affyread(fullfile(exampleDataDir,'Ecoli-antisense-121502.CHP'),libDir)

%%
% The |ProbeSets| field contains information about the probe sets. This
% includes some library information, such as the ID and the type of probe
% set, and also results information such as the calculated signal value and
% the |Present/Absent/Marginal| call information. The call is given in the
% |Detection| field of the ProbeSets structure. The 'argG_b3172_at' probe
% set is called as being |'Present'|.

chpStruct.ProbeSets(5213)
            
%%
% However, the 'IG_2069_3319273_3319712_rev_at' probe set is called |'Absent'|.

chpStruct.ProbeSets(5216)
            
%%
% And the 'yhbX_b3173_at' probe set is called 'Marginal'.

chpStruct.ProbeSets(5215)
%%
% You can calculate how many probe sets are called as being |'Present'|,

numPresent = sum(strcmp('Present',{chpStruct.ProbeSets.Detection}))

%%
% |'Absent'|,

numAbsent = sum(strcmp('Absent',{chpStruct.ProbeSets.Detection}))
%%
% and |'Marginal'|.

numMarginal = sum(strcmp('Marginal',{chpStruct.ProbeSets.Detection}))

%%
% |maboxplot| will display a box plot of the log2 signal values for all
% probe sets.

maboxplot(chpStruct,'Signal','title',chpStruct.Name)

%% Library Files (CDF Files)
% The CHP file gives summary information about probe sets but if you want
% more detailed information about how the individual probes in a probe set
% behave you need to connect the probe information in the CEL file to the
% corresponding probe sets. This information is stored in the CDF library
% file associated with a chip type. The CDF files are typically stored in a
% central library directory. 

cdfStruct = affyread('Ecoli_ASv2.CDF',libDir)

%%
% Most of the information in the file is about the probe sets. In this
% example there are 7312 regular probe sets and 13 QC probe sets. The
% |ProbeSets| field of the structure is a 7325x1 array of structures.

cdfStruct.ProbeSets

%%
% A probe set record contains information about the name, type and number
% of probe pairs in the probe set.
probeSetIndex = 5213;
cdfStruct.ProbeSets(probeSetIndex)

%%
% The information about where the probes for a probe set are on the chip is
% stored in the |ProbePairs| field. This is a matrix with one row for each
% probe pair and six columns. The information in the columns corresponds to
% the |ProbeSetColumnNames| of the CDF structure.

cdfStruct.ProbeSetColumnNames
cdfStruct.ProbeSets(probeSetIndex).ProbePairs

%%
% The first column shows the probe group number. The second column shows
% the probe direction. The group number is always 1 for expression arrays.
% Direction 1 corresponds to 'sense' and 2 corresponds to 'anti-sense'.
% The remaining columns give the X and Y coordinates of the PM and MM
% probes on the chip. You can use these coordinates to find the index of a 
% probe in the celStruct.

PMX = cdfStruct.ProbeSets(probeSetIndex).ProbePairs(1,3);
PMY = cdfStruct.ProbeSets(probeSetIndex).ProbePairs(1,4);
theProbe = find((celStruct.Probes(:,1) == PMX) & ...
                       (celStruct.Probes(:,2) == PMY))
    
%%
% You can then extract all the information about this probe from the CEL
% structure.

celStruct.Probes(theProbe,:)

%%
% If you want to do this lookup for all probes, you can use the function
% |probelibraryinfo|. This creates a matrix with one row per probe and
% three columns. The first column is the index of the probe set to which
% the probe belongs. The second column contains the probe pair index and
% the third column indicates if the probe is a perfect match (1) or
% mismatch (-1) probe. Notice that index of the probe pair index is 1
% based.

probeinfo = probelibraryinfo(celStruct,cdfStruct);

probeinfo(theProbe,:)

%%
% The function |probesetvalues| does the reverse of this lookup and creates
% a matrix of information from the CEL and CDF structures containing all
% the information about a given probe set. This matrix has 20 columns
% corresponding to |ProbeSetNumber|, |ProbePairNumber|, |UseProbePair|,
% |Background|, |PMPosX|, |PMPosY|, |PMIntensity|, |PMStdDev|, |PMPixels|,
% |PMOutlier|, |PMMasked|, |MMPosX|, |MMPosY|, |MMIntensity|, |MMStdDev|,
% |MMPixels|, |MMOutlier|, |MMMasked|, |Group|, and |Direction|.
probeName = cdfStruct.ProbeSets(probeSetIndex).Name;
psvals = probesetvalues(celStruct,cdfStruct,probeName);
sprintf( ['%4d %2d %d %d  PM: %3d %3d %5.1f %5.1f %2d %d %d',...
          '  MM: %3d %3d %5.1f %5.1f %2d %d %d %d %d\n'],psvals')

%%
% You can extract the intensity values from the matrix and look at some of
% the statistics of the data.

pmIntensity = psvals(:,7);
mmIntensity = psvals(:,14);
boxplot([pmIntensity,mmIntensity],'labels',{'Perfect Match','Mismatch'})
title(sprintf('Boxplot of raw intensity values for probe set %s',...
    probeName),'interpreter','none')
% Use interpreter none to prevent the TeX interpreter treating the _ as
% subscript.

%% Plotting the Probe Set Values
% Now that you have the intensity values for the probes, you can plot the
% values for the perfect match and mismatch probes.

figure
plot(pmIntensity,'b'); hold on
plot(mmIntensity,'r'); hold off
title(sprintf('Probe intensity values for probe set %s',...
    probeName),'interpreter','none')
%% 
% Alternatively, you can use the function |probesetplot| to create this
% plot directly from the CEL and CDF structures. The showstats option adds
% the mean, and lines for +/- one standard deviation for both the perfect
% match and the mismatch probes to the plot.

probesetplot(celStruct,cdfStruct,probeName,'showstats',true);

%% Gene Names and Probe Set IDs
% The Affymetrix probe set IDs are not particularly descriptive. The
% mapping between the probe set IDs and the gene IDs is stored in the GIN library
% file. This is a text file so you can open it in an editor and browse
% through the file, or you can use |affyread| to read the information into
% a structure.

ginStruct = affyread('Ecoli_ASv2.GIN',libDir)

%%
% You can search through the structure for a particular probe set.
% Alternatively, you can use the function |probesetlookup| to find
% information about the gene for a probe set. 

info = probesetlookup(cdfStruct,probeName)

%% Getting Sequence Information About a Probe Set
% The function |probesetlink| will link out to the NetAffx(TM) Web site to show
% the actual sequences used for the probes. Note that you will need to be a
% registered user of NetAffx to access this information.

probesetlink(cdfStruct,probeName);

%%
% Affymetrix, GeneChip, and NetAffx are registered trademarks of
% Affymetrix, Inc.

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20AFFYDEMO%20in%20Bioinformatics%20Toolbox%204.2 Suggest an enhancement for this example.>*

displayEndOfDemoMessage(mfilename)

