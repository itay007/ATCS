function absvar = exprprofvar(data,varargin)
%EXPRPROFVAR calculates the variance of expression profiles.
%
%   EXPRPROFVAR(DATA) calculates the variance of each expression profile in
%   dataset DATA. DATA can be a MATLAB numeric array or a DataMatrix
%   object.
%
% 	If no output arguments are specified, a histogram bar plot of the
% 	variance is displayed.
%     
%   EXPRPROFVAR(...,'SHOWHIST',TF) displays a histogram of the variance of
%   data if TF is true. 
%
%   Example:
%
%       load yeastdata
%       datavar = exprprofvar(yeastvalues,'showhist',true);
%
%   See also EXPRPROFRANGE, GENERANGEFILTER, GENEVARFILTER.

% Copyright 2003-2008 The MathWorks, Inc.


showhist = false;
if nargout == 0
    showhist = true;
end

if nargin > 1
    if rem(nargin,2)== 0
        error(message('bioinfo:exprprofvar:IncorrectNumberOfArguments', mfilename));
    end
    
    okargs = {'showhist',''};
    for j=1:2:nargin-2
        [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
        
        switch(k)
            case 1  % showhist
                showhist = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        end   
    end
end 

absvar = var(data,[],2);

%== Output DataMatrix
if nargout > 0 && isa(data, 'bioma.data.DataMatrix') 
    absvar = bioma.data.DataMatrix(absvar, data.RowNames, []);
end

if showhist
    if isa(data, 'bioma.data.DataMatrix')
        numbuckets = max(10,numel(absvar, ':', ':')/100);
    else
        numbuckets = max(10,numel(absvar)/100);
    end
    
    numbuckets = ceil(min(numbuckets,100));
    hist(absvar,numbuckets);
    title('Profile Variances');
end
