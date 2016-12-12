function bowtiebuild(varargin)
%BOWTIEBUILD generates an index using Burrows-Wheeler transform
%
%   BOWTIEBUILD(INPUT, INDEXBASENAME) builds an index with name
%   INDEXBASENAME using the reference sequence(s) in INPUT. INDEXBASENAME
%   is a string containing the path and basename of the resulting BOWTIE
%   index files. INPUT is either a string or a cell string indicating one
%   or more FASTA formatted files with the reference sequences to be
%   indexed.
%
%   BOWTIEBUILD(..., 'BOWTIEBUILDOPTIONS', OPTIONS) specifies BOWTIEBUILD
%   additional options. Type bowtiebuild('--help') for available options.
%
%   Note: Some pre-built index files for model organisms can be downloaded
%   directly from the BOWTIE repository http://bowtie-bio.sourceforge.net/.
%
%   Example:
%
%     % Download the E. coli genome from NCBI:
%     getgenbank('NC_008253','tofile','NC_008253.fna','SequenceOnly',true)
%
%     % Build a BOWTIE index with the basename ECOLI:
%     bowtiebuild('NC_008253.fna','ECOLI')
%
%   More information about the BOWTIE algorithm (Version 0.12.7) can be
%   found at http://bowtie-bio.sourceforge.net/.
%
%   See also BAMINFO, BIOMAP, BOWTIE, FASTAINFO, FASTQINFO, SAMREAD,
%   SAMINFO. 

%   Copyright 2012 The MathWorks, Inc.

if ispc 
    error(message('bioinfo:bowtie:UnsupportedPlatform',upper(mfilename)));
end

if nargin > 1
    
    % Parse optional PVPs and/or set defaults
    bowtieBuildOptions = parse_inputs(varargin{3:end});
    
    % Get INPUT and uniformize to comma separated list
    inputReference = varargin{1};
    if iscellstr(inputReference)
        inputReference = regexprep(sprintf('%s,',inputReference{:}),',$','');
    elseif ~ischar(inputReference) || ~isrow(inputReference)
        error(message('bioinfo:bowtie:InvalidReference'))
    end
    
    % Get INDEXBASENAME
    indexBaseName = varargin{2};
    
    % Prepare argument list for bowtie.
    % bowtie-build [options]* <reference_in> <ebwt_outfile_base>  
    args = cat(2, {'--quiet'}, bowtieBuildOptions, inputReference, indexBaseName);
   
else
    % Undocummented: If there is one input argument, it is passed in bowtie
    % options as if they were typed at the command line in a shell.
    narginchk(1,1);
    argsCell = textscan(varargin{1}, '%s');
    args = argsCell{1}';
end

% Preface the argument list with the string 'bowtie_build' to select a call
% to bowtie.
args = cat(2, {'bowtie_build'}, args);

try
    bioinfoprivate.bowtie_mex(args);
catch e
    throw(e);    
end
end

function [bowtieBuildOptions] = parse_inputs(varargin)
% Parse input PV pairs.

% default value
bowtieBuildOptions = {};

if rem(nargin, 2) ~= 0
    error(message('bioinfo:bowtie:IncorrectNumberOfArguments', mfilename));
end

possible_arguments = {'BowtieBuildOptions'};

for j=1:2:nargin-1
    [arg_index, parameter_value] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, possible_arguments, mfilename);
    
    switch(arg_index)
        case 1  % BowtieBuildOptions
            if ~ischar(parameter_value) || isempty(parameter_value)
                error(message('bioinfo:bowtie:OptionsMustBeString'));
            end
            parsedOptions = textscan(parameter_value, '%s');
            bowtieBuildOptions = parsedOptions{1}';
    end
end
end
