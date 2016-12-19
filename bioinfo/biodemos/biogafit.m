function classPerformance = biogafit(thePopulation,Y,id)
%BIOGAFIT The fitness function for BIOGAMSDEMO
%
%   This function uses the classify function to measure how well mass
%   spectrometry data is grouped using certain masses. The input argument
%   thePopulation is a vector of row indices from the mass spectrometry
%   data Y. Classification performance is a linear combination of the error
%   rate and the posteriori probability of the classifier. 

%   Copyright 2003-2005 The MathWorks, Inc.


thePopulation = round(thePopulation);
try
   [c,e,p]  = classify(Y(thePopulation,:)',Y(thePopulation,:)',double(id));
   cp = classperf(id,c);
   classPerformance = 100*cp.ErrorRate + 1 - mean(max(p,[],2));
catch
   classPerformance = Inf;
end


