% MSSEQPROCESSING Script to create OvarianCancerQAQCdataset.mat (used in
% CANCERDETECTDEMO). Before running this file initialize the variable
% "repository" to the full path where you placed you mass-spectrometry
% files. For Example:
%
%   repository = 'F:/MassSpecRepository/OvarianCD_PostQAQC/';
% 
% or
%
%   repository = '/home/username/MassSpecRepository/OvarianCD_PostQAQC/';
%
% The approximate time of execution is 18 minutes (Pentium 4, 4GHz). If you
% have the Parallel Computing Toolbox refer to BIODISTCOMPDEMO to see
% how you can speed this analysis up.

%   Copyright 2003-2008 The MathWorks, Inc.


repositoryC = [repository 'Cancer/'];
repositoryN = [repository 'Normal/'];

filesCancer = dir([repositoryC '*.txt']);
NumberCancerDatasets = numel(filesCancer);
fprintf('Found %d Cancer mass-spectrograms.\n',NumberCancerDatasets)
filesNormal = dir([repositoryN '*.txt']);
NumberNormalDatasets = numel(filesNormal);
fprintf('Found %d Control mass-spectrograms.\n',NumberNormalDatasets)

files = [ strcat('Cancer/',{filesCancer.name}) ...
          strcat('Normal/',{filesNormal.name})];
N = numel(files);   % total number of files

fprintf('Total %d mass-spectrograms to process...\n',N)

[MZ,Y] = msbatchprocessing(repository,files);

disp('Finished; normalizing and saving to OvarianCancerQAQCdataset.mat.')
Y = msnorm(MZ,Y,'QUANTILE',0.5,'LIMITS',[3500 11000],'MAX',50);

grp = [repmat({'Cancer'},size(filesCancer));...
       repmat({'Normal'},size(filesNormal))];
   
save OvarianCancerQAQCdataset.mat Y MZ grp
