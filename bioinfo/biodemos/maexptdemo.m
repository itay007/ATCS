%% Working with Objects for Microarray Experiment Data
% This example shows how to create and manipulate MATLAB(R)'s containers
% designed for storing data from a microarray experiment.

%   Copyright 2009-2012 The MathWorks, Inc.

%% Containers for Gene Expression Experiment Data
% Microarray experimental data are very complex, usually consisting of data
% and information from a number of different sources. Storing and managing
% the large and complex data sets in a coherent manner has always being a
% challenge. Bioinformatics Toolbox(TM) provides a set of objects to
% represent the different pieces of data from a microarray experiment. 

%%
% It is easier to manage all the data from a microarray experiment if the
% different pieces can be organized and stored into a single data
% structure. The |ExpressionSet| class is a single convenient data
% structure for storing and coordinating the different data objects from a
% microarray gene expression experiment. 

%%
% An |ExpressionSet| object consists of these four components that are
% common to all microarray gene expression experiments: 
% 
% _Experiment data:_ Expression values from microarray experiments. These
% data are stored as an instance of the |ExptData| class.
% 
% _Sample information:_ The metadata describing the samples in the
% experiment. The sample metadata are stored as an instance of the
% |MetaData| class.
%
% _Array feature annotations:_ The annotations about the features or probes
% on the array used in the experiment. The annotations can be stored as an
% instance of the |MetaData| class.
%
% _Experiment descriptions:_ Information to describe the experiment methods
% and conditions. The information can be stored as an instance of the
% |MIAME| class. 

%%
% The |ExpressionSet| class coordinates and validates these data
% components. The class provides methods for retrieving and setting the
% data stored in an |ExpressionSet| object. An |ExpressionSet| object also
% behaves like many other MATLAB data structures that can be subsetted
% and copied. 

%% Experiment Data
% In a microarray gene expression experiment, the measured expression
% values for each feature per sample can be represented as a
% two-dimensional (2D) matrix. The matrix has _F_ rows and _S_ columns,
% where _F_ is the number of features on the array, and _S_ is the number
% of samples on which the expression values were measured. A |DataMatrix|
% object is designed to contain this type of data. It is a 2D matrix with
% row and column names. A |DataMatrix| object can be indexed not only by
% its row and column numbers, logical vectors, but also by its row and
% column names. But linear indexing is not supported. 
% 
% For example, create a Datamatrix with row and column names:
dm = bioma.data.DataMatrix(rand(5,4), 'RowNames','Feature', 'ColNames', 'Sample')

%% 
% The function |size| returns the number of rows and columns in a
% |DataMatrix| object. 
size(dm)

%%
% You can index into a |DataMatrix| object like other MATLAB numeric arrays
% by using row and column numbers. For example, access the elements at rows
% 1 and 2, column 3 of |dm|:
dm(1:2, 3)

%%
% You can also index into a |DataMatrix| object by using its row and column
% names. Reassign the elements in row 2 and 3, column 1 and 4 to different
% values:
dm({'Feature2', 'Feature3'}, {'Sample1', 'Sample4'}) = [2, 3; 4, 5]

%%
% The example gene expression data used in this example is a small set of
% data from a microarray experiment profiling adult mouse gene expression
% patterns in common strains on the Affymetrix(R) MG-U74Av2 array [1]. The
% file |mouseExprsData.txt| contains the small set of expression values in
% a table format.

%%
% Read the expression values from the file |mouseExprsData.txt| into MATLAB
% Workspace as a |DataMatrix| object: 
exprsData = bioma.data.DataMatrix('file', 'mouseExprsData.txt');

%%
class(exprsData)
%%
% Get the properties of the |DataMatrix| object, |exprsData|.
get(exprsData)

%%
% Check the sample names:
colnames(exprsData)

%%
% View the first 10 rows and 5 columns: 
exprsData(1:10, 1:5)

%%
% Many of the basic MATLAB array operations also work with a |DataMatrix|
% object. For example, you can log2 transform the expression values:
exprsData_log2 = log2(exprsData);
%%
% View the first 10 rows and 5 columns
exprsData_log2(1:10, 1:5)

%%
% Change the |Name| property to be descriptive about |exprsData_log2|:
exprsData_log2 = set(exprsData_log2, 'Name', 'Log2 Based mouseExprsData');
%%
get(exprsData_log2)

%%
% In a microarray experiment, the data set often contains one or more
% matrices that have the same number of rows and columns and identical row
% names and column names, like in two-color microarray experiments.
% |ExptData| class is designed to contain and coordinate one or more data
% matrices with same dimensional properties, i.e. same dimension size,
% identical row names and column names. The data values are stored in an
% |ExptData| object as |DataMatrix| objects.  Each |DataMatrix| object is
% considered an element in an |ExptData| object. The |ExptData| class is
% responsible for data validation and coordination between these
% |DataMatrix| objects. For the purposes of this example, you will store
% the gene expression data of natural scale and log2 base expression values
% separately in an instance of |ExptData| class. 
mouseExptData = bioma.data.ExptData(exprsData, exprsData_log2,...
                    'ElementNames', {'natualExprs', 'log2Exprs'})
                
%%
% Access a |DataMatrix| element in |mouseExptData| using the element name. 
exprsData2 = mouseExptData('log2Exprs');

%%
get(exprsData2)

%%
% |ExptData| does not allow input matrices of different size or
% |DataMatrix| objecs with different row or column names. It would error in
% following case.
try
    mouseExptData = bioma.data.ExptData(exprsData, dm,...
        'ElementNames', {'naturalExprs', 'log2Exprs'})
catch ME
    disp(ME.message)
end

%% Sample Metadata
% The metadata about the samples in a microarray experiment can be
% represented as a table with _S_ rows and _V_ columns, where _S_ is the
% number of samples, and _V_ is the number of variables. The contents of
% the table are the values for each sample per variable. For example, the
% file |mouseSampleData.txt| contains such a table. Alternately, this table
% of variable values can be stored in a |dataset| array.
%% 
% Users often find that simple column names do not provide enough
% information about the variables. What is the name supposed to represent?
% What units are the variables measured in? Another table can contain such
% description metadata about variables. In this table of metadata, rows
% represent variables and at least one column contains a description of
% each variables. For example, the file |mouseSampleData.txt| contains
% descriptions about the sample variables (The lines are each prefaced with
% a _#_ symbol. The metadata about the variables can also be stored in a
% |dataset| array.

%%
% The |MetaData| class is designed for storing and manipulating variable
% values and their metadata in a coordinated fashion. You can read the
% |mouseSampleData.txt| file into MATLAB as a |MetaData| object.
sData = bioma.data.MetaData('file', 'mouseSampleData.txt', 'vardescchar', '#')

%%
% The properties of |MetaData| class provide information about the size and
% dimension labels. There are 26 rows of samples and 5 columns of variables
% in the example sample data file.
sData.NSamples

%%
sData.NVariables

%%
% The variable values and the variable descriptions for the samples are
% stored as two |dataset| arrays in a |MetaData| class. The |MetaData|
% class provides access methods to the variable values and the meta
% information describing the variables. Access the sample metadata using
% the |variableValues| method.
sData.variableValues

%%
% View a summary of the sample metadata and the variables it contains.
summary(sData.variableValues)

%%
% The |sampleNames| and |variableNames| methods are convenient ways to
% access the names of samples and variables. Retrieve the variable names of
% the |sData| object.
variableNames(sData)

%%
% You can retrieve the meta information about the variables describing the
% samples using the |variableDesc| method. In this example, it contains
% only the descriptions about the variables.
variableDesc(sData)

%%
% You can subset the sample data |sData| object the same way as a |dataset|
% array. 
sData(3:6, :)

%%
% To see the mouse strain of the 2nd and 14th samples.
sData.Strain([2 14])

%%
% Note that the row names in |sData| and the column names in |exprsData|
% are the same. It is important feature of relationship between the
% expression data and the sample data in the same experiment. 
all(ismember(sampleNames(sData), colnames(exprsData)))

%% Feature Annotation Metadata
% The gene expression data in the example is obtained using an Affymetrix
% MG-U74Av2 array. The metadata about the features or probe set on an array
% can be very large and diverse, and important for the experiment. The chip
% makers usually provide a specific annotation file about the features of
% each type of array. The metadata can be stored as a |MetaData| object for
% a specific experiment. In this example, the annotation file for the
% MG-U74Av2 array can be downloaded from the Affymetrix web site. Download
% the file and read it into MATLAB as a |dataset| array. Header lines in
% the annotation file have been manually removed. Alternatively you can use
% the |Range| option in the |dataset| constructor. MATLAB removes blank
% spaces in the variable names (columns) so they may be used as valid
% MATLAB identifiers; this is indicated in the command prompt by displaying
% a warning.
mgU74Av2 = dataset('xlsfile', 'MG_U74Av2_annot.csv');

%%
% Inspect the properties of this |dataset| array.
get(mgU74Av2)

%%
% Check and see the number of probe set IDs in the annotation file.
numel(mgU74Av2.ProbeSetID)

%%
% Retrieve the names of variables describing the features on the MG-U74Av2
% array and view the first 20 variable names. 
fDataVariables = get(mgU74Av2, 'VarNames');
fDataVariables(1:20)'

%%
% Set the |ObsNames| property to the probe set IDs, this allows you to
% access individual gene annotations by indexing the |dataset| with probe
% set IDS.
mgU74Av2 = set(mgU74Av2,'ObsNames',mgU74Av2.ProbeSetID);
mgU74Av2('100709_at',{'GeneSymbol','ChromosomalLocation'})

%%
% In many cases, not all the information read from the array annotation
% file is useful, it is better to store only the annotation information
% applicable to the experiment. In this example, extract annotations 
% |GeneTitle|, |GeneSymbol|, |ChromosomalLocation|, and |Pathway| for the
% features unique to the data in |exprsData|. 
mgU74Av2 = mgU74Av2(:,{'GeneTitle',...
                       'GeneSymbol',...
                       'ChromosomalLocation',...
                       'Pathway'});

%%
% Because the expression data in this example is only a small set of the
% full expression values, you will work with only the features in the
% |exprsData| |DataMatrix| object. Find the matching features in
% |exprsData|.  
mgU74Av2 = mgU74Av2(rownames(exprsData),:);
get(mgU74Av2)

%%
% You can store the feature annotation |dataset| array as an instance of
% the |MetaData| class.
fData = bioma.data.MetaData(mgU74Av2)

%%
% Notice that there are not descriptions for the feature variables in the
% |fData| |MetaData| object. You can add descriptions about the variables
% in |fData| using the |variableDesc| method.
fData = variableDesc(fData, {'Gene title of a probe set',...
                             'Probe set gene symbol',...
                             'Probe set chromosomal locations',...
                             'The pathway the genes involved in'})

%% Experiment Information
% The |MIAME| class is a flexible data container designed for a collection
% of basic descriptions about a microarray experiment, for instance,
% investigators or laboratory where the experiment was done, and
% description about the array designs. The |MIAME| class is designed to be
% light-weight and loosely follows the Minimum Information About a
% Microarray Experiment (MIAME) specification [2]. The information can be
% accessed through the 14 properties of the |MIAME| class. 

%%
% Create a |MIAME| object by providing some basic information.
expDesc = bioma.data.MIAME('investigator', 'Jane OneName',...
                           'lab',          'Bioinformatics Laboratory',...
                           'title',        'Example Gene Expression Experiment',...
                           'abstract',     'An example of using microarray objects.',...
                           'other',        {'Notes: Created from a text files.'})
%%
% Another way to create an |MIAME| object is from GEO series data. The
% |MIAME| class will populate the corresponding properties from the data
% structure returned by the |getgeodata| function. Create an |MIAME| object
% for the experiment information about the mouse gene profile experiment in
% the example. The dataset is available in the GEO database with a series
% accession number of |GSE3327| [1]. Note: The |GSE3327| dataset is quite
% large it takes some time to download.
geoSeries = getgeodata('GSE3327')

%%
exptGSE3327 = bioma.data.MIAME(geoSeries)
%%
% View the abstract of the experiment and its PubMed IDs.
exptGSE3327.Abstract

%%
exptGSE3327.PubMedID

%% Assembling an ExpressionSet Object
% The |ExpressionSet| class is designed specifically for microarray gene
% expression experiment data. Assemble an |ExpressionSet| object for the
% example mouse gene expression experiment from the different data objects
% you just created. 
exptSet = bioma.ExpressionSet(exprsData, 'SData', sData,...
                                         'FData', fData,...
                                         'Einfo', exptGSE3327)

%%
% You can also create an |ExpressionSet| object with only the expression
% values in a |DataMatrix| or a numeric matrix.
miniExprSet = bioma.ExpressionSet(exprsData)

%% Saving and Loading an ExpressionSet Object
% The data objects for a microarray experiment can be saved as _MAT_ files.
% Save the |ExpressionSet| object |exptSet| to a _MAT_ file named
% |mouseExpressionSet.mat|.
save mouseExpressionSet exptSet

%%
% Clear all the variables from the MATLAB Workspace.
clear all

%%
% Load the _MAT_ file |mouseExpressionSet| into the MATLAB Workspace.
load mouseExpressionSet

%%
% Inspect the loaded |ExpressionSet| object.
exptSet.elementNames

%%
exptSet.NSamples

%%
exptSet.NFeatures

%% Accessing Data Components of an ExpressionSet Object
% A number of methods are available to access and update data stored in an
% |ExpressionSet| object. In this example, you will explore some of the
% data access methods and basic operations of the |ExpressionSet| class. 
%%
% You can also access the columns of the sample data using dot notation.
exptSet.Strain(1:5)

%%
% Retrieve the feature names using the |featureNames| method. In this
% example, the feature names are the probe set identifiers on the array.
featureNames(exptSet, 1:5)

%%
% The unique identifier of the samples can be accessed via the
% |sampleNames| method.
exptSet.sampleNames(1:5)

%%
% The |sampleVarNames| method lists the variable names in the sample data.
exptSet.sampleVarNames
%%
% Extract the |dataset| array containing sample information.
sDataset = sampleVarValues(exptSet)

%%
% Retrieve the |ExptData| object containing expression values. There are
% can be more than one |DataMatrix| object with identical dimensions in an
% |ExptData| object. While in an |ExpressionSet| object, there is always a
% element |DataMatrix| object named |Expressions| containing the expression
% matrix.
exptDS = exptData(exptSet)

%%
% Extract only the expression |DataMatrix| instance.
dMatrix = expressions(exptSet);

%%
% The returned expression |DataMatrix| should be identical to the
% |exprsData| |DataMatrix| object that you created earlier. 
get(dMatrix)

%%
% Get PubMed IDs for the experiment stored in |exptSet|.
exptSet.pubMedID

%% Subsetting an ExpressionSet Object
% Subsetting is a very useful operation. Subsetting an |ExpressionSet|
% object is very similar to subsetting a |DataMatrix| object or a |dataset|
% array. The first indexing argument subsets the features and the second
% argument subsets the samples. For example, you can create a new
% |ExpressionSet| object consisting of the first five features and the
% samples named |A|, |B|, and |C|.
mySet = exptSet(1:5, {'A', 'B', 'C'})

%%
% Inspect the subset |mySet|
size(mySet)

%%
featureNames(mySet)
%%
sampleNames(mySet)

%%
% Another example is to create a subset consisting of only the samples from
% hippocampus tissues.
hippocampusSet = exptSet(:, nominal(exptSet.Source)== 'hippocampus')

%%
hippocampusSet.Source

%%
hippocampusExprs = expressions(hippocampusSet);

%%
get(hippocampusExprs)

%%
% You can find more details about the basic operations and available
% methods for the microarray experiment data objects in the help and
% reference pages.

%% References
% [1] Hovatta, I., Tennant, R. S., Helton, R., et al. "Glyoxalase 1 and
%     glutathione reductase 1 regulate anxiety in mice", Nature, 438, 
%     pp 662-666, 2005.
%%
% [2] http://www.mged.org/Workgroups/MIAME/miame_1.1.html.

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20MAEXPTDEMO%20in%20Bioinformatics%20Toolbox%204.2 Suggest an enhancement for this example.>*

displayEndOfDemoMessage(mfilename)
