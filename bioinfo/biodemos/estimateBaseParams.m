function varargout = estimateBaseParams(countDM, sFactors, outOpt)
% Return estimate of the mean, the variance of counts. And the variance
% functions for predicting how much variance one should expect for counts
% at a certain level.
% countDM  - a DataMatrix containing the count table
% sFactors - DataMatrix containing the size factor for the experimental
%            condition.

% Transfer to common scale
baseCounts = dmbsxfun(@rdivide, countDM, sFactors);
%  Compute base mean - the average of the counts from a sample
baseMean = mean(baseCounts, 2);
% Compute variance - sample variance on the common scale
baseVar = var(baseCounts, 0, 2);

outFlag = 1;
switch outOpt
    case 'MeanAndVar'
       outFlag = 1;
    case 'SmoothFunc'
       outFlag = 2;
    case 'Diagnostic'
       outFlag = 3; 
end

if outFlag == 1
   varargout{1} = baseMean;
   varargout{2} = baseVar;
   return;
end

zidx = baseMean > 0;


% Remove the genes with zero means
baseVar_nz = baseVar(zidx);
baseMean_nz = baseMean(zidx);

% Try LOWESS to get the smooth function
varSmoothF_X = log(baseMean_nz);
varSmoothF_Y = malowess(varSmoothF_X, baseVar_nz);

% Uniqueness required for interp1
[varSmoothF_X_U, uidx] = unique(varSmoothF_X);
varSmoothF_Y_U = varSmoothF_Y(uidx);

if outFlag == 2    
    varargout{1} = varSmoothF_X_U;
    varargout{2} = varSmoothF_Y_U;
    return;
end

if outFlag == 3
    df = size(countDM, 2) - 1;
    raw_var_fit = interp1(varSmoothF_X_U, varSmoothF_Y_U, log(baseMean),...
                          'linear', 0);
    varRatio = baseVar ./ raw_var_fit;
    varargout{1} = chi2cdf(df * varRatio, df);
    return;
end
end