function bowtie(varargin)
%BOWTIE maps short reads to a reference sequence using the Burrows-Wheeler transform.
%
%   BOWTIE(INDEXBASENAME, READS, OUTPUTFILENAME) aligns the reads specified
%   in READS to the indexed reference specified with INDEXBASENAME and
%   writes the results to the BAM formatted file OUTFILENAME. INDEXBASENAME
%   is a string containing the path and basename of the BOWTIE index files.
%   READS is either a string or a cell string indicating one or more FASTQ
%   formatted files with the input reads. The file name extension .bam is 
%   is automatically added to OUTFILENAME when missing.
%
%   BOWTIE(..., 'BAMFILEOUTPUT', FALSE) outputs alignment results in a SAM
%   formatted file. The file name extension .sam is automatically added to
%   OUTFILENAME when missing. BAMFILEOUTPUT defaults to true.
%
%   BOWTIE(..., 'PAIRED', FALSE) paired-read alignment is performed using
%   the odd elements in READS as the upstream mates and the even elements
%   in READS as the downstream mates.
%
%   BOWTIE(..., 'BOWTIEOPTIONS', OPTIONS) specifies BOWTIE additional
%   options. Type bowtie('--help') for available options. 
%
%   Example:
%
%     % Please see example in BOWTIEBUILD help on how to create the BOWTIE
%     % index files for the E. Coli genome. Some pre-built index files for
%     % model organisms can be downloaded directly from the BOWTIE
%     % repository http://bowtie-bio.sourceforge.net/.
%
%     % An example of E. Coli short reads are provided by the FASTQ file
%     % ecoli100.fq: 
%     fastqfile = which('ecoli100.fq')
%  
%     % Align short reads in ecoli100.fq to the index with basename 'ECOLI':
%     bowtie('ECOLI', fastqfile , 'ecoli100.bam')
%
%     % Access the mapped reads using BioMap:
%     bm = BioMap('ecoli100.bam')
%
%   More information about the BOWTIE algorithm (Version 0.12.7) can be
%   found at http://bowtie-bio.sourceforge.net/.
%
%   See also BAMINFO, BIOMAP, BOWTIEBUILD, FASTAINFO, FASTQINFO, SAMREAD,
%   SAMINFO.

%   Copyright 2012 The MathWorks, Inc.

if ispc 
    error(message('bioinfo:bowtie:UnsupportedPlatform',upper(mfilename)));
end
if nargin > 1
    if nargin < 3
        error(message('bioinfo:bowtie:IncorrectNumberOfArguments',upper(mfilename)));
    end

     % Parse optional PVPs and/or set defaults
    [bowtieOptions, BAMFileOutput, pairedInput] = parse_inputs(varargin{4:end});

    % Get INDEXBASENAME
    indexBaseName = varargin{1};
    
    % Get READS and uniformize to cellstr
    reads = varargin{2};
    if ischar(reads) && isrow(reads)
        reads = {reads};
    elseif ~iscellstr(reads)
        error(message('bioinfo:bowtie:InvalidReads'))
    end
   
    % Validate the extension of the output file.
    [outfilePath, outfileName, outfileExt] = fileparts(varargin{3});
    if ~isempty(outfileExt)
        if BAMFileOutput && ~strcmp(outfileExt, '.bam')
            error(message('bioinfo:bowtie:BAMFileExtensionRequired'));
        elseif ~BAMFileOutput && ~strcmp(outfileExt, '.sam')
            error(message('bioinfo:bowtie:SAMFileExtensionRequired'));
        end
    end
    samFileName = fullfile(outfilePath,[outfileName,'.sam']);
    
    % Prepare input file list for bowtie.
    if pairedInput
        if rem(numel(reads),2)==1 || isempty(reads)
            error(message('bioinfo:bowtie:InvalidPairedReads'))
        end
        % prepare two comma separated lists
        reads = { regexprep(sprintf(',%s',reads{1:2:end}),'^,','-1'),... 
                  regexprep(sprintf(',%s',reads{2:2:end}),'^,','-2') };
    else
        % prepare a single comma separated list
        reads = regexprep(sprintf('%s,',reads{:}),',$','');
    end
    
    % Prepare argument list for bowtie.
    % bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]    
    args = cat(2, {'-S', '--quiet'}, bowtieOptions, indexBaseName, reads, samFileName);
    
else
    % Undocummented: If there is one input argument, it is passed in bowtie
    % options as if they were typed at the command line in a shell.
    narginchk(1,1);
    argsCell = textscan(varargin{1}, '%s');
    args = argsCell{1}';
end

% Preface the argument list with the string 'bowtie' to select a call to
% bowtie.
args = cat(2, {'bowtie'}, args);

try   
    bioinfoprivate.bowtie_mex(args);
catch e
    infoStrStart = regexp(e.message,'Info:','once');
    if isempty(infoStrStart)
        throw(e);
    else
        disp(e.message(infoStrStart:end))
    end
end

% Postprocessing. Since bowtie does not support output in BAM format (only
% SAM) we use an internal utility to convert to BAM and order the file if
% so requested by the user.
if nargin > 1 && BAMFileOutput
    % If we were called with a MathWorks signature (rather than just with a
    % string of bowtie command line options) then we have the option of
    % providing users with sorted BAM files. The output name of the BAM
    % file will be the same as that specified for SAM.
    [pathstr, name, ~] = fileparts(samFileName);
    bioinfoprivate.bamaccessmex('sam2bam', samFileName, fullfile(pathstr,[name,'.bam']));
    
    % This next call to bamaccessmex automatically adds the .bam extension.
    bioinfoprivate.bamaccessmex('bamsort', fullfile(pathstr,[name,'.bam']),fullfile(pathstr,name));
end
end

function [bowtieOptions, BAMFileOutput, pairedInput] = parse_inputs(varargin)
% Parse input PV pairs.

% default values
pairedInput = false;
BAMFileOutput = true;
bowtieOptions = {};

if rem(nargin, 2) ~= 0
    error(message('bioinfo:bowtie:IncorrectNumberOfArguments', mfilename));
end

possible_arguments = {'BowtieOptions', 'BAMFileOutput', 'Paired'};

for j=1:2:nargin-1
    [arg_index, parameter_value] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, possible_arguments, mfilename);
    
    switch(arg_index)
        case 1  % BowtieOptions
            if ~ischar(parameter_value) || isempty(parameter_value)
                error(message('bioinfo:bowtie:OptionsMustBeString'));
            end
            parsedOptions = textscan(parameter_value, '%s');
            bowtieOptions = parsedOptions{1}';
         case 2 % BAMFileOutput
            BAMFileOutput = bioinfoprivate.opttf(parameter_value,possible_arguments{arg_index},mfilename);
         case 3 % Paired
            pairedInput = bioinfoprivate.opttf(parameter_value,possible_arguments{arg_index},mfilename);
    end
end
end
