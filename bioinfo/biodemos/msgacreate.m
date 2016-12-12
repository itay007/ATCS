function pop = msgacreate(GenomeLength,FitnessFcn,options,numPoints)
%MSGACREATE Population creation function for MSGADEMO
%   
%   This function creates a population matrix with dimensions of 
%   options.PopulationSize rows by the number of independent variables 
%   (GenomeLength) columns. These values are integers that correspond to
%   randomly selected rows of the mass spectrometry data Y. Each row of the
%   population matrix is a random sample of row indices of the mass spec
%   data
%   Note: This function's input arguments are required by the genetic
%   algorithm solver. See the Global Optimization Toolbox documentation for
%   further detail. 

%   Copyright 2003-2009 The MathWorks, Inc.


NumberOfIntervals = ceil(numPoints/GenomeLength);
totalPopulationSize = sum(options.PopulationSize);
range = 1:NumberOfIntervals-1:numPoints; 
indiv = zeros(1,GenomeLength); 
pop = zeros(totalPopulationSize,length(range)-1);
for i = 1:totalPopulationSize  %Every individual
    for j = 1: length(range)-1
        limit = range(j+1) - range(j);
        indiv(j) = range(j) + floor(rand*limit); %choose a random
    end
    pop(i,:) = indiv; 
end

