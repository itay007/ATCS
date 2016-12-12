function expr = affygcrma(probeData, CDFFile, Seqs, varargin)
%AFFYGCRMA Performs GC Robust Multi-array Average (GCRMA) procedure.
%
%   EXPRESSION = AFFYGCRMA(CELFILES, CDFFILE, SEQS) reads Affymetrix
%   CELFILES, the associated library CDFFILE, and the associated sequence
%   file, processes the probe intensities using GCRMA procedures, and then
%   outputs the log2 based expression value for each probe set to a
%   DataMatrix object, EXPRESSION, where each row corresponds to a probe
%   set, and each column corresponds to a CEL file. The probe set IDs are
%   the row names of EXPRESSION, and the CEL file names are the column
%   names of EXPRESSION. CELFILES can be a string specifying a CEL file
%   name or a cell array of CEL file names. CDFFILE is a string specifying
%   a CDF file name. If CELFILES is set to wildcard '*', then all the CEL
%   files in the current directory will be read. If CELFILES or CDFFILE is
%   set to ' ', an Open File dialog opens for you to select CEL files or a
%   CDF file interactively. SEQS can be a probe sequence data file name or
%   a sequence matrix returned by the AFFYPROBESEQREAD function.
% 
%   EXPRESSION = AFFYGCRMA(PROBEINTENSTRUCT, SEQS) uses the GCRMA procedure
%   to processes the probe intensity structure returned by the
%   CELINTENSITYREAD function. SEQS is a sequence matrix returned by the
%   AFFYPROBESEQREAD function.
%
%   EXPRESSION = AFFYGCRMA(..., 'CELPATH', CELPATH) allows you to specify
%   the directory where the CELFILES are stored.
%
%   EXPRESSION = AFFYGCRMA(..., 'CDFPATH', CDFPATH) allows you to specify
%   the directory where the CDFFILE is stored.
% 
%   EXPRESSION = AFFYGCRMA(..., 'SEQPATH', SEQPATH) allows you to specify
%   the directory where the probe sequence data file SEQS is stored.
% 
%   EXPRESSION = AFFYGCRMA(...,'CHIPINDEX', CID) computes affinities with
%   sequence information and mismatch probe intensities from the chip
%   specified by index CID. Default CID = 1 for using the first chip.
% 
%   EXPRESSION = AFFYGCRMA(...,'OPTICALCORR', TF) performs optical
%   background correction on the input probe intensities if TF is set to
%   true. Default is set to TRUE.
% 
%   EXPRESSION = AFFYGCRMA(...,'CORRCONST', RHO) sets the correlation
%   coefficient constant RHO of the bivariate-normal distribution of the
%   log non-specific binding background intensity for each PM and MM probe
%   pair. The value should be between 0 and 1. Default RHO = 0.7.
% 
%   EXPRESSION = AFFYGCRMA(...,'METHOD', METHOD) sets an estimate method to
%   estimate the signal. METHOD can be the faster, but ad hoc, Maximum
%   Likelihood Estimate 'MLE' (default), or 'EB', for the slower and more
%   formal Empirical-Bayes type estimator. 
% 
%   EXPRESSION = AFFYGCRMA(...,'TUNINGPARAM', TP) sets a tuning parameter
%   for the estimate method. This tuning parameter sets the lower bound of
%   signal values with positive probability. The value must be greater than
%   or equal to 0. Default is 5 for MLE or 0.5 for EB. See DOC GCRMA for
%   more details.
% 
%   EXPRESSION = AFFYGCRMA(...,'GSBCORR', TF) performs gene specific
%   binding (GSB) correction with probe affinity data if TF is set to true.
%   Default is TRUE. If there is no probe affinity information, this
%   option is ignored.
% 
%   EXPRESSION = AFFYGCRMA(...,'MEDIAN',true) takes the median of the
%   ranked values instead of the mean for normalization. For more
%   information, see the help for QUANTILENORM. 
% 
%   EXPRESSION = AFFYGCRMA(..., 'OUTPUT', TYPE), returns output expression
%   values in the form defined by TYPE. TYPE can be 'log', 'log2'
%   (default), 'log10', 'linear' scale, or a function handle. If a
%   function handle is passed, the output is transformed by the
%   transformation defined by the function.
% 
%   AFFYGCRMA(...,'SHOWPLOT',SP) displays optical adjusted log2 based MM
%   intensities vs. its affinities and LOWESS fit for computing NSB data of
%   the chip specified by index SP. When SP is TRUE, the plot of the first
%   chip in PMS is displayed. For more information, see the help for
%   GCRMABACKADJ.
%
%   AFFYGCRMA(...,'VERBOSE',false) turns off the display of preprocessing
%   progress. Default is true.
% 
%   Example:
%
%       % Process expression values from the CEL files in current directory
%       % and a CDF file in % a specified directory
%       expr = affygcrma('*','HG_U95Av2.CDF', 'HG-U95Av2_probe_tab'...
%                      'cdfpath', 'D:\Affymetrix\LibFiles\HGGenome');
% 
%   See also  AFFYPROBEAFFINITIES, AFFYRMA, CELINTENSITYREAD, GCRMA,
%   GCRMABACKADJ, MAFDR, MATTEST, QUANTILENORM, RMASUMMARY.

%   Copyright 2008-2009 The MathWorks, Inc.


%   AFFYGCRMA(...,'STRETCH', S) sets stretch correction constant. It can
%   be used for correcting results after MLE or EB background adjustments.
%   S can be a logical or a numerical number. No stretch correction If
%   S=FALSE (default). If S=TRUE, the correction factor is set to 1.15 for
%   MLE method and 1 for EB method. A numerical number will be the stretch
%   correction factor for the method used in background adjustment. 

%==Check input
bioinfochecknargin(nargin,2,mfilename);

%==Check if the first input is a structure
if isstruct(probeData) 
    if ~checkProbeStruct(probeData)
        error(message('bioinfo:affygcrma:InvalidProbeStructure'));
    end

    if nargin == 2 
        Seqs = CDFFile;
    elseif nargin > 2
        varargin = {Seqs, varargin{:}};
        Seqs = CDFFile;
    end
    
    if isstruct(Seqs) && isfield(Seqs, 'SequenceMatrix')
        Seqs = Seqs.SequenceMatrix;
    elseif isa(Seqs, 'uint8')
    else
        error(message('bioinfo:affygcrma:InvalidSeqMatrix'));
    end
else
    bioinfochecknargin(nargin,3,mfilename);
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
        [probeData, CDFFile] = celintensityread(probeData, CDFFile,...
                                     'celpath', inStruct.celpath,...
                                     'cdfpath', inStruct.cdfpath,...
                                     'verbose', inStruct.verbose,...
                                     'pmonly', false);
    catch ME
        bioinfoprivate.bioerrorrethrow(mfilename, ME);
    end 
end

if isstruct(Seqs) && isfield(Seqs, 'SequenceMatrix')
    Seqs = Seqs.SequenceMatrix;
end
    
%== Read Sequence file
if ~isa(Seqs, 'uint8')
    try
        if inStruct.verbose
            disp('Read sequence file.');
        end
        
        Seqs = affyprobeseqread(Seqs, CDFFile, 'seqpath', inStruct.seqpath,...
                                'cdfpath', inStruct.cdfpath, 'seqonly', true);
        Seqs = Seqs.SequenceMatrix;
    catch ME
        bioinfoprivate.bioerrorrethrow(mfilename, ME);
    end
end

%== AFFYPROBEAFFINITY
try
    if inStruct.verbose
        disp('Calculating probe affinities.');
    end
   
    [apm, amm] = affyprobeaffinities(Seqs,...
                        probeData.MMIntensities(:,inStruct.chipindex),...
                        'ProbeIndices', probeData.ProbeIndices);
catch ME
    bioinfoprivate.bioerrorrethrow(mfilename, ME);
end

% If not output assignment skipt the rest.
if nargout == 0
   return; 
end

%== GCRMA
try
    if isempty(inStruct.showplot)
        probeData.PMIntensities = gcrmabackadj(probeData.PMIntensities,...
            probeData.MMIntensities, apm, amm, ...
            'opticalcorr', inStruct.opticalcorr,...
            'corrconst',   inStruct.corrconst,...
            'method',      inStruct.method,...
            'tuningparam', inStruct.tuningparam,...
            'gsbcorr',     inStruct.gsbcorr,...
            'verbose',     inStruct.verbose,...
            'stretch',     inStruct.stretch);
    else
        probeData.PMIntensities = gcrmabackadj(probeData.PMIntensities,...
            probeData.MMIntensities,apm, amm, ...
            'opticalcorr', inStruct.opticalcorr,...
            'corrconst',   inStruct.corrconst,...
            'method',      inStruct.method,...
            'tuningparam', inStruct.tuningparam,...
            'gsbcorr',     inStruct.gsbcorr,...
            'showplot',    inStruct.showplot,...
            'verbose',     inStruct.verbose,...
            'stretch',     inStruct.stretch);
    end
    
catch ME
    bioinfoprivate.bioerrorrethrow(mfilename, ME);
end

%== quantilenorm
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

end % affygcrma

%------------------Hlper functions----------------
function inStruct = parse_inputs(varargin)
% Parse input PV pairs.

% Check for the right number of inputs
if rem(nargin,2)== 1
     error(message('bioinfo:affygcrma:IncorrectNumberOfArguments', mfilename))
end

% Allowed inputs
okargs = {'celpath', 'cdfpath','seqpath',...
          'chipindex','opticalcorr', 'method', 'corrconst',...
          'tuningparam', 'gsbcorr', 'median', 'output',...
          'showplot', 'verbose','stretch'};

% Defaults
inStruct.celpath = pwd;         % celintensityread
inStruct.cdfpath = pwd;          % celintensityread
inStruct.seqpath = pwd;         % affyprobeseqread
inStruct.chipindex = 1;         % gcrma
inStruct.opticalcorr = true;    % gcrma
inStruct.corrconst = 0.7;        % gcrma
inStruct.method = 'mle';        % gcrma
inStruct.tuningparam = [];      % gcrma
inStruct.gsbcorr = true;        % gcrma
inStruct.stretch = false;       % gcrma
inStruct.median = false;        % quantilenorm
inStruct.output = @log2;        % rmasummary
inStruct.showplot = [];         % affyprobeaffinities
inStruct.verbose = true;        % all

for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1 % cel path
            inStruct.celpath = pval;
        case 2 % cdfpath
            inStruct.cdfpath = pval;
        case 3 % seqpath
            inStruct.seqpath = pval;
        case 4 % chipidx
            inStruct.chipidx = pval;
        case 5 % opticalcorr
            inStruct.opticalcorr = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 6 % method
            if ischar(pval)
                okmethods = {'mle', 'eb'};
                nm = find(strncmpi(pval, okmethods, numel(pval))); %#ok
                
                if isempty(nm)
                    error(message('bioinfo:affygcrma:UnknownMethodName'));
                else
                    inStruct.method = pval;
                end
                
            else
                error(message('bioinfo:affygcrma:MethodNameNotValid'));
            end    
        case 7 % corrconst
            if ~isnumeric(pval)
                error(message('bioinfo:affygcrma:CorrConstMustBeNumeric'));
            elseif pval < 0 || pval > 1
                error(message('bioinfo:affygcrma:CorrConstMustBe0To1'));
            end
            
            inStruct.corrconst = pval;
        case 8 % tunningparam
            if ~isnumeric(pval) || ~isscalar(pval)
                error(message('bioinfo:affygcrma:TuningParamNotSingleNumericValue'));
            elseif (pval < 0 )
                error(message('bioinfo:affygcrma:BadTuningParameterValue'));
            else
                inStruct.tuningparam = pval;
            end
        case 9 % gsbcorr
            inStruct.gsbcorr = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 10  % median
            inStruct.median = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 11 % output
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
                        error(message('bioinfo:affygcrma:incompatibleOutputScale'));
                    case 5
                        inStruct.output = @(x)x; %No op
                end
            elseif isa(outputtype, 'function_handle')
                inStruct.output = outputtype;
            else
                error(message('bioinfo:affygcrma:OutputTypeNotFunctionHandle'));
            end
        case 12 % showplot
            inStruct.showplot = pval;
        case 13 % verbose
            inStruct.verbose = bioinfoprivate. opttf(pval,okargs{k},mfilename);
        case 14 % stretch
            inStruct.stretch = pval;   
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
         isfield(probeStruct,'PMIntensities') && ...
         isfield(probeStruct,'MMIntensities');
end % checkProbeStruct function
