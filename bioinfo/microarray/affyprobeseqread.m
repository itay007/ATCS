function S = affyprobeseqread(seqfile, cdf, varargin)
%AFFYPROBESEQREAD reads a data file describing the probe sequences on an
% Affymetrix GeneChip.
% 
%   S = AFFYPROBESEQREAD(SEQFILE, CDF) reads the probe sequence data file,
%   SEQFILE, of an Affymetrix GeneChip, and creates a structure, S. CDF can
%   be a structure or file name for CDF file. AFFYPROBESEQREAD can read a
%   tabulator-separated file with one row per probe, and a FASTA file with
%   one header for each probe. The returned structure with these fields
%       ProbeSetIDs 
%       ProbeIndices 
%       SequenceMatrix (int) - A=1, C=2, G=3, T=4, None = 0
%   The ProbeSetIDs are ordered as in CDF file, and ProbeIndices are
%   ordered according to the ProbesetIDs. For some probes that don't have
%   sequence information, the correspondent rows in SequenceMatrix should
%   contain all zeros.
% 
%   AFFYPROBESEQREAD(..., 'SEQPATH', SEQPATH) allows you to specify the
%   directory where the SEQFILE is stored.
% 
%   AFFYPROBESEQREAD(..., 'CDFPATH', CDFPATH) allows you to specify the
%   directory where the CDF file is stored.
% 
%   AFFYPROBESEQREAD(..., 'SEQONLY', TF) returns a structure with only a
%   SequenceMatrix field if TF set to true. Default is FALSE.
% 
%   Notes: 
%       1) The column names in the tabular file are: Probe set ID, Probe X,
%       Probe Y, probe interrogation position Probe Sequence, Target
%       strandedness. The file does not contain information of the chip
%       type. Please make sure that the file matches the CDF chip type
%       before using this function. 
%       2) The probe sequence FASTA file contains information of the chip
%       type, probeset ids, probe x, y, sequences in each probeset, and the
%       interrogation_position. 
%   
%   Example:
% 
%       % Read in a sequence fasta file, and the array's CDF library file
%       S = affyprobeseqread('HG-U95A_probe_fasta', 'HG_U95A.CDF');
%       
%       % Read in a sequence file in a specified directory, and a CDF 
%       % structure already in workspace. Your file may be in a different
%       % location.
%      
%       S = affyprobeseqread('HG-U95A_probe_tab',hgu95aCDFStruct,...
%                      'seqpath','C:\Affymetrix\SequenceFiles\HGGenome');
%       % To get the sequences of the first probeset in S.SequenceMatrix
%       seq = int2nt(S.SequenceMatrix(1:20, :))
% 
%   Sequence files and CDF library files are available from
%   http://www.affymetrix.com/support/technical/byproduct.affx
%   
%   See also AFFYGCRMA, AFFYPREPROCESSDEMO, AFFYPROBEAFFINITIES, AFFYREAD,
%   CELINTENSITYREAD, GCRMA, GCRMABACKADJ, PROBELIBRARYINFO,
%   PROBESETLOOKUP, PROBESETVALUES.
% 
%   GeneChip and Affymetrix are registered trademarks of Affymetrix, Inc.

% Copyright 2006-2012 The MathWorks, Inc.


bioinfochecknargin(nargin,2,mfilename);

S =[];

if ~ischar(seqfile) || (iscellstr(seqfile)&& length(seqfile) > 1) 
    error(message('bioinfo:affyprobeseqread:SequenceFileNameNotAString', seqfile));
end

cdf_struct = [];
cdffile =[];
if isstruct(cdf) && isfield(cdf,'ChipType')
    cdf_struct = cdf;
elseif ischar(cdf) || (iscellstr(cdf)&& length(cdf) == 1)
    cdffile = cdf;
    if strcmpi(computer,'SOL64')
        error(message('bioinfo:affyprobeseqread:ReadCDFFilesNoSolaris'));
    end
else
    name = inputname(2);
    error(message('bioinfo:affyprobeseqread:NotACDFStructureORFileName', name));
end

% Initialization
seqpath = pwd; 
cdfpath = '';
seqonly = false;

% Deal with the various inputs
if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:affyprobeseqread:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'seqpath', 'cdfpath', 'seqonly'};
    for j=1:2:nargin-3
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs, numel(pname)));
        if isempty(k)
            error(message('bioinfo:affyprobeseqread:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:affyprobeseqread:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % seqpath
                    if ~ischar(pval)
                        error(message('bioinfo:affyprobeseqread:SeqPathNotAString', upper( pname )));
                    end
                    seqpath = pval;
                case 2 % cdfpath
                    if ~ischar(pval)
                        error(message('bioinfo:affyprobeseqread:CDFPathNotAString', upper( pname )));
                    end
                    cdfpath = pval;
                case 3 % seqonly flag
                    seqonly = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            end
        end
    end
end

% Check seqfile is a file
fullseqfile = [seqpath, filesep, seqfile];
if  (exist(fullseqfile,'file') || exist(fullfile(pwd,fullseqfile),'file'))
    % continue
else
    error(message('bioinfo:affyprobeseqread:CannotReadInput', fullseqfile));
end

chip_type = [];
probe_ids = [];
probe_xs = [];
probe_seqs = [];

%-----------------------------------------------
% Guess the file format, let FASTAREAD check it first
try
    seq_struct = fastaread(fullseqfile);
    [chip_type, probe_ids, probe_xs, probe_seqs] = handlefastatext(seq_struct);
    fastaflag = true;
catch theException
    % When fastaread errors and it is determined that the file is readable
    % but it is not fasta formatted, error is not given since later we try
    % to read it as just a tabular sequence file. 
    eid = theException.identifier;
    if strcmpi(eid, 'bioinfo:fastaread:FastaNotValid')
        fastaflag = false;
    else
        msgId = 'bioinfo:affyprobeseqread:InvalidFastaFile';
        newException =  MException(msgId,'%s',getString(message(msgId,fullseqfile)));
        throw(addCause(newException,theException))
    end
end

if ~fastaflag
    % Tabular sequence file
    try
        fid = fopen(fullseqfile,'rt');
        % Pass  the empty lines
        pos = ftell(fid);
        while isempty(fgetl(fid))
            pos = ftell(fid);
        end
        fseek(fid, pos, -1);

        % Check for Probe at the beginning of the first line
        tline = fgetl(fid);
        if isempty(regexp(tline, '^Probe','once'))
            tabseqflag = false;
            fclose(fid);
        else
            tabseqflag = true;
            c = textscan(fid, '%s%d%*d%*d%s%*[^\n]');
            fclose(fid);
            probe_ids = c{1};
            probe_xs = c{2};
            probe_seqs = c{3};
        end
    catch theException %#ok<NASGU>
        error(message('bioinfo:affyprobeseqread:CannotReadInput', fullseqfile));
    end  
    if ~tabseqflag
        error(message('bioinfo:affyprobeseqread:UnknowSequenceFileFormat', fullseqfile));
    end
end

% Check CDF 
if isempty(cdf_struct) && ~isempty(cdffile)
    % Read cdf file
    fullcdf = [cdfpath, filesep, cdffile];
    try
        cdf_struct = affyread(fullcdf);
    catch theException
        msgId = 'bioinfo:affyprobeseqread:InvalidCDFFile';
        newException =  MException(msgId,'%s',getString(message(msgId,fullcdf)));
        throw(addCause(newException,theException))
    end    
end

% Fasta file contains chip_type information, compare with CDF chiptype
if ~isempty(chip_type)
    chip_type = chip_type{:};
    if ~strcmpi(chip_type, cdf_struct.ChipType)
        warning(message('bioinfo:affyprobeseqread:ChipTypeNameNotMatch', chip_type, cdf_struct.ChipType));
    end
end

% Build output structure
nProbeSets = cdf_struct.NumProbeSets;
nProbes = sum([cdf_struct.ProbeSets.NumPairs]);
probesetIDs = {cdf_struct.ProbeSets.Name}';

% Loop through the probesets to get sequences for each probe, and conver the
% sequences into a nProbes x 25 matrix 
seqmatrix = zeros(nProbes, 25, 'uint8');
probeIndices = zeros(nProbes,1, 'uint8');

% Unique probe_ids so less strcmpi time
[~, I_first] = unique(probe_ids, 'first');
[unique_probe_ids, I_last] = unique(probe_ids, 'last');

probeCount = 0;    
for i = 1:nProbeSets
    numPairs = cdf_struct.ProbeSets(i).NumPairs;
    thePairs = cdf_struct.ProbeSets(i).ProbePairs;
    theId = probesetIDs{i};
    
    id = strcmpi(theId, unique_probe_ids);
    
    probex = probe_xs(I_first(id):I_last(id));
    probeseq = probe_seqs(I_first(id):I_last(id));
    
    if isempty(probeseq) % if the probe set sequences are not provided
        seqmatrix(probeCount+1 : probeCount + numPairs, :) = 0;
    elseif size(probeseq, 1) ~= numPairs
        cdf_probeX = thePairs(:,3);
        for j = 1:numPairs
            xloc = find(probex == cdf_probeX(j));
            if isempty(xloc)
                seqmatrix(probeCount+j, :) = 0;
            else
                seqmatrix(probeCount+j, :) = nt2int(char(probeseq(xloc)));
            end
        end
    else
        % Order probeseq with the order in CDF file
        [dummy, loc] = ismember(probex, thePairs(:,3));
        probeseq(loc) = probeseq;

        seqmatrix(probeCount+1 : probeCount + numPairs, :) = nt2int(char(probeseq));
    end

    probeIndices(probeCount+1 : probeCount + numPairs) = (0:numPairs-1)';
    probeCount = probeCount + numPairs;
end

if ~seqonly
    S.ProbesetIDs = probesetIDs;
    S.ProbeIndices = probeIndices;
end

S.SequenceMatrix = seqmatrix;

%------------------------------------------------------------
function [chip_type, probe_ids, probe_xs, probe_seqs] = handlefastatext(seq_struct)
% Parsing the fasta sequence file

probe_headers = {seq_struct.Header};
probe_seqs = {seq_struct.Sequence}';

N = length(probe_headers);

probe_ids = cell(N,1);
probe_xs = zeros(N,1);

for i = 1:N
    [chip_type, probe_ids(i), probe_xs(i)]= strread(probe_headers{i},...
                     '%*s%s%s%d%*[^\n]', 'delimiter', ':');
end







