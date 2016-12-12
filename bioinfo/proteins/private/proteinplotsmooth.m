function [data, dataInd] = proteinplotsmooth(data,smoothMethod,ws,edgeval,nterm,cterm)
%PROTEINPLOTSMOOTH helper function for proteinplot and proteinpropplot.
%
% PROTEINPLOTSMOOTH smooths the input data using the supplied method and
% options.


%   Copyright 2003-2006 The MathWorks, Inc.

mid = (ws + 1)/2;
switch smoothMethod
    case 'linear'
        fb = ones(1,ws);
        fs = linspace(edgeval,1,mid);
        fb(1:mid) = fs;
        fb(mid:ws) = fliplr(fs);
        fb = fb / sum(fb);
        data = filter(fb,1,data);
        data = data(mid:(end - mid + 1),:);
        dataInd = (nterm + (ws-1)/2):(cterm - (ws-1)/2);
    case  'exponential'
        fb = ones(1,ws);
        fb(1) = edgeval; fb(ws) = edgeval;
        if edgeval < 1e-6
            edgeval = 1e-6;
        end
        k = log(edgeval) /  (mid - 1);
        fb(2:mid-1) = exp(k) .^((mid-2):-1:1);
        fb((mid+1): (ws-1)) = fliplr(fb(2:mid-1));
        fb = fb / sum(fb);
        data = filter(fb,1,data);
        data = data(mid:(end - mid + 1),:);
        dataInd = (nterm + (ws-1)/2):(cterm - (ws-1)/2);
    case 'lowess'
        for n = 1:size(data,2),
            data(:,n) = bioinfoprivate.masmooth([],data(:,n),ws,'lowess');
        end
        dataInd = nterm:cterm;
end
