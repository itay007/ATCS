function classPerformance = msgafit(thePopulation,Y,id,numSamples)
%MSGAFIT The fitness function for MSGADEMO
%
%   This function uses the knnclassify function from the Bioinformatics
%   Toolbox to measure how well mass spectrometry data is grouped. The
%   input argument thePopulation is a vector of row indices from the mass
%   spectrometry data Y. Data points from these rows are grouped by k
%   nearest-neighbor classification and compared against the known
%   classification (id). The number of correctly classified rows is
%   returned as a negative percentage of correctly classified points.

%   Copyright 2005 The MathWorks, Inc.


thePopulation = uint16(thePopulation);
c = knnclassify(Y(thePopulation,:)',Y(thePopulation,:)', ...
    double(id),5,'corr','consensus');
classPerformance = sum(c==id)/numSamples; % get percent of total properly classified
classPerformance = -classPerformance;     % reverse the sign in order to minimize
