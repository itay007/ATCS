function expr = affyrma(probeData, CDFFile, varargin)
%AFFYRMA Performs Robust Multi-array Average (RMA) procedures.
%
%   EXPRESSION = AFFYRMA(CELFILES, CDFFILE) reads Affymetrix CELFILES and
%   the associated library CDFFILE, processes the probe intensities using
%   RMA background adjustment, normalization and summarization procedures,
%   and then outputs the log2 based expression value for each probe set to
%   a DataMatrix object, EXPRESSION, where each row corresponds to a probe
%   set, and each column corresponds to a CEL file. The probe set IDs are
%   the row names of EXPRESSION, and the CEL file names are the column
%   names of EXPRESSION. CELFILES can be a string specifying a CEL filename
%   or a cell array of CEL file names. CDFFILE is a string of the CDF file
%   name. If CELFILES is set to wildcard '*', then all the CEL files in the
%   current directory will be read. If CELFILES or CDFFILE is set to '', an
%   Open File dialog opens for you to select CEL files or a CDF file
%   interactively. 
% 
%   EXPRESSION = AFFYRMA(PROBEINTENSTRUCT) uses the RMA procedures to
%   process the probe intensity structure returned by the CELINTENSITYREAD
%   function. 
%
%   EXPRESSION = AFFYRMA(..., 'CELPATH', CELPATH) allows you to specify the
%   directory where the CELFILES are stored.
%
%   EXPRESSION = AFFYRMA(..., 'CDFPATH', CDFPATH) allows you to specify the
%   directory where the CDFFILE is stored.
% 
%   EXPRESSION = AFFYRMA(...,'METHOD',METHOD) sets the estimation method
%   for the background adjustment model parameters. METHOD can be 'RMA'
%   (default) or 'MLE'. For more details of the estimation methods See the
%   help for RMABACKADJ.
% 
%   EXPRESSION = AFFYRMA(...,'TRUNCATE',TF) uses truncated Gaussian as
%   background noise model by default. For more details of the background
%   noise model, see the help for RMABACKADJ.
% 
%   EXPRESSION = AFFYRMA(...,'MEDIAN', true) takes the median of the ranked
%   values instead of the mean for normalization. For more information, see
%   the help for QUANTILENORM.
% 
%   EXPRESSION = AFFYRMA(..., 'OUTPUT', TYPE), returns output expression
%   values in the form defined by TYPE. TYPE can be 'log', 'log2'
%   (default), 'log10', 'linear' scale, or a function handle. If a
%   function handle is passed, the output is transformed by the
%   transformation defined by the function.
% 
%   AFFYRMA(...,'SHOWPLOT',',SP) plots the histogram of the background
%   adjusted log2(PM) values, and the convoluted probability distribution
%   function with estimated mu, sigma and alpha for the chip specified by
%   index SP. SP may also be a vector contains the indices in CELFILES or
%   PROBEINTENSTRUCT. Use 'all' for SP to show plots for all the chips. For
%   more information, see the help for RMABACKADJ. 
%
%   AFFYRMA(...,'VERBOSE',false) turns off the display of the progress of
%   the preprocessing process. Default is true.
% 
%   Example:
%
%       % Process expression values from the CEL files in current directory
%       % and a CDF file in % a specified directory
%       expr = affyrma('*','HG_U95Av2.CDF',...
%                      'cdfpath', 'D:\Affymetrix\LibFiles\HGGenome');
% 
%   See also  AFFYGCRMA, CELINTENSITYREAD, GCRMA, MAFDR, MATTEST,
%   QUANTILENORM, RMABACKADJ, RMASUMMARY.

%   Copyright 2008 The MathWorks, Inc.


%==Check input
bioinfochecknargin(nargin,1,mfilename);

%==Check if the first input is a structure
if isstruct(probeData) 
    if ~checkProbeStruct(probeData)
        error(message('bioinfo:affyrma:InvalidProbeStructure'));
    end
    
    if nargin > 1
        varargin = {CDFFile, varargin{:}};
    end
else
    bioinfochecknargin(nargin,2,mfilename);
end

%== Parse input options
inStruct = parse_inputs(varargin{:});

%== If not output and not plot, just return;
if nargout == 0 && isempty(inStruct.showplot)
   return; 
end

%== Read CEL files
if ~isstruct(probeData)
    try
        probeData = celintensityread(probeData, CDFFile,...
                                     'celpath', inStruct.celpath,...
                                     'cdfpath', inStruct.cdfpath,...
                                     'verbose', inStruct.verbose);
    catch ME
        bioinfoprivate.bioerrorrethrow(mfilename, ME);
    end 
end

%== RMABACKADJ
try
    if inStruct.verbose
        disp('Perform background adjustment.');
    end
    
    if isempty(inStruct.showplot)
        probeData.PMIntensities = rmabackadj(probeData.PMIntensities,...
            'method', inStruct.method,...
            'truncate', inStruct.truncate);
    else  
        probeData.PMIntensities = rmabackadj(probeData.PMIntensities,...
            'method', inStruct.method,...
            'truncate', inStruct.truncate,...
            'showplot', inStruct.showplot);
    end
catch ME
    bioinfoprivate.bioerrorrethrow(mfilename, ME);
end

% If not output assignment skipt the rest.
if nargout == 0
   return; 
end

%== Quantilenorm
try
    if inStruct.verbose
        disp('Normalizing.');
    end
    
    probeData.PMIntensities = quantilenorm(probeData.PMIntensities,...
                                           'median', inStruct.median);
catch ME
    bioinfoprivate.bioerrorrethrow(mfilename, ME);
end

%==RMASUMMARY
try 
    if inStruct.verbose
        disp('Calculating expression.');
    end
    expr = rmasummary(probeData.ProbeIndices, probeData.PMIntensities,...
                      'output', inStruct.output);
catch ME
    bioinfoprivate.bioerrorrethrow(mfilename, ME);
end

%== Output
if isfield(probeData, 'CELNames')
    expr = bioma.data.DataMatrix(expr, 'RowNames', probeData.ProbeSetIDs,...
                                       'ColNames', probeData.CELNames);
else
    expr = bioma.data.DataMatrix(expr, 'RowNames', probeData.ProbeSetIDs);
end

end % affyrma

%------------------Hlper functions----------------
function inStruct = parse_inputs(varargin)
% Parse input PV pairs.

% Check for the right number of inputs
if rem(nargin,2)== 1
     error(message('bioinfo:affyrma:IncorrectNumberOfArguments', mfilename))
end

% Allowed inputs
okargs = {'celpath', 'cdfpath',...
          'method', 'truncate',...
          'median', 'output',...
          'showplot', 'verbose'};

% Defaults
inStruct.celpath = pwd;         % celintensityread
inStruct.cdfpath = '';          % celintensityread
inStruct.method = 'rma';        % rmabackadj
inStruct.truncate = true;       % rmabackadj
inStruct.median = false;        % quantilenorm
inStruct.output = @log2;       % rmasummary
inStruct.showplot = [];          % rmabackadj
inStruct.verbose = true;        % all

for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1 % cel path
            inStruct.celpath = pval;
        case 2 % cdfpath
            inStruct.cdfpath = pval;
        case 3 % method
            if ischar(pval)
                okmethods = {'rma','mle'};
                nm = strmatch(lower(pval), okmethods);
                if isempty(nm)
                    error(message('bioinfo:affyrma:MethodNameNotValid'));
                elseif length(nm)>1
                    error(message('bioinfo:affyrma:AmbiguousMethodName', pval));
                else
                    inStruct.method = pval;
                end
            else
                error(message('bioinfo:affyrma:MethodNameNotValid'));
            end
        case 4 % truncate
            inStruct.truncate =  bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 5 % median
            inStruct.median =  bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 6 % output
            outputtype = pval;
            if ischar(outputtype)
                okmethods = {'log','log2','log10','natural','linear'};
                output = bioinfoprivate.pvpair(pval, lower(outputtype), okmethods, mfilename);
                switch output
                    case 1
                        inStruct.output = @log;
                    case 2
                        inStruct.output = @log2;
                    case 3
                        inStruct.output = @log10;
                    case 4 %'natural' is deprecated (R2012b)
                        error(message('bioinfo:affyrma:incompatibleOutputScale'));
                    case 5
                        inStruct.output = @(x)x; %No op
                end
            elseif isa(outputtype, 'function_handle')
                inStruct.output = outputtype;
            else
                error(message('bioinfo:affyrma:OutputTypeNotFunctionHandle'));
            end
        case 7 % showplot
            inStruct.showplot = pval;
        case 8 % verbose
            inStruct.verbose = bioinfoprivate. opttf(pval,okargs{k},mfilename);
    end
end
end %parse_input

function valid = checkProbeStruct(probeStruct)
% Checks for probe intensity structure field names

valid =  isfield(probeStruct,'NumChips') && ...
         isfield(probeStruct,'NumProbeSets') && ...
         isfield(probeStruct,'NumProbes') && ...
         isfield(probeStruct,'ProbeSetIDs') && ...
         isfield(probeStruct,'ProbeIndices') && ...
         isfield(probeStruct,'GroupNumbers') && ...
         isfield(probeStruct,'PMIntensities');
end % checkProbeStruct function
