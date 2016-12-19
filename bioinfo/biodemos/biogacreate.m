function pop = biogacreate(GenomeLength,FitnessFcn,options,Y,id)
%BIOGACREATE Population creation function for MSGADEMO
%   
%   This function creates a population matrix with dimensions of 
%   options.PopulationSize rows by the number of independent variables 
%   (GenomeLength) columns. These values are integers that correspond to
%   randomly selected rows of the mass spectrometry data Y. Each row of the
%   population matrix is a random sample of row indices of the mass spec
%   data.
%   Note: This function's input arguments are required by the GA function.
%   See GA documentation for further detail.

%   Copyright 2005 The MathWorks, Inc.


ranked_features = rankfeatures(Y,id);
pop = zeros(options.PopulationSize,GenomeLength);
pop(:) = randsample(ranked_features(1:numel(pop)),numel(pop));
