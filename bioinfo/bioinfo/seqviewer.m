function varargout = seqviewer(seq, varargin)
%SEQVIEWER visualize biological sequences.
%
%   SEQVIEWER is an interactive tool for viewing biological sequences.
%
%   SEQVIEWER(SEQ) loads a sequence SEQ into the GUI. SEQ can be a structure
%   with a field Sequence, a character array, or a filename with an
%   extension of .gbk, .gpt, .fasta, .fa, or .ebi.
%
%   SEQVIEWER(...,'ALPHABET',ALPHA) opens a sequence of alphabet ALPHA.
%   SEQVIEWER uses 'AA' as default except when all symbols in the sequence
%   are in {'A' 'C' 'G' 'T' '-'}, then it uses 'NT'. Use 'AA' to force for
%   an amino acid sequence.
%
%   Examples:
%       S = getgenbank('M10051')
%       seqviewer(S)
%
%       % open two sequences in a FASTA file
%       seqviewer('hexaNT.fasta')
%
%   See also AA2NT, AACOUNT, AMINOLOOKUP, BASECOUNT, BASELOOKUP,
%   DIMERCOUNT, EMBLREAD, FASTAREAD, FASTAWRITE, GENBANKREAD, GENETICCODE,
%   GENPEPTREAD, GETEMBL, GETGENBANK, GETGENPEPT, NT2AA, PROTEINPLOT, 
%   SEQALIGNVIEWER, SEQCOMPLEMENT, SEQDISP, SEQRCOMPLEMENT, SEQREVERSE,
%   SEQSHOWORFS, SEQSHOWWORDS, SEQWORDCOUNT.

% Copyright 2003-2012 The MathWorks, Inc.

% This function is the gateway to the sequence viewer.

%   SEQVIEWER(...,'ACC',accession) accession can be
%   genbank,genpept, or embl.

%   SEQVIEWER(...,'FROMVIEWER',true) a flag to indicate that the sequence viewer
%   is calling this function.

% SEQVIEWER passes a structure to svMatlab with these possible fields:
%       Alpha - the type of sequence {'AMINO', 'NT'}
%       FileType -  what the source is {'genbank','genpept','embl','otherstruct','char'}
%       PreambleText - the preamble info from get genbankread, etc...
%       LocusName - the name for the window header
%       Definition - the name for the sequence in the viewer
%       Features - any features that might be
%       OriginalData - the original variable sent to this function
%       Sequence - The sequence
 
if ~jvm_available
    error(message('bioinfo:seqviewer:JavaComponentMissing'));
end

% All code is in a try-catch block because if the function was called by
% the seqviewer GUI, then the error message is just handed over as an output
% argument instead of throwing an error and terminating.
try
    
    %defaults
    msg = '';
    varn_empty = ''; % default name for empty input
    varn = 'Sequence'; %default name for char arrays that are sequences
    alpha = '';
    download = '';
    from_viewer = false;
    
    if(nargin == 0) %start seqviewer with no sequences
        % If there is a viewer present no need to open a blank page.
        if ~(seqviewer.SVUtil.svMatlab('manage_bioinfodesktop', 'state'))
            seqviewer.SVUtil.svMatlab('create_ui',build_empty(varn_empty));
        end
        return;
    elseif nargin == 1 && ischar(seq) && strncmpi(seq, 'close', 5)
        % Undocumented option to close viewer programmatically
        seqviewer.SVUtil.svMatlab('manage_bioinfodesktop', 'close');
        return;
        
    elseif nargin > 2
        nvarargin = numel(varargin);
        
        if rem(nvarargin,2) == 1 || nargin > 5
            error(message('bioinfo:seqviewer:NumberOfArgumentsIncorrect', mfilename));
        end
        
        okargs = {'alphabet', 'acc', 'fromviewer'};
        for j=1:2:nvarargin-1
            [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
            switch(k)
                case 1 % alphabet
                    if  ~bioinfoprivate.optAlphabet(pval,okargs{k}, mfilename)
                        alpha = 'nt';
                    end
                case 2  % acc
                    [~,download] = bioinfoprivate.optPartialMatch(pval,{'embl','genpept','genbank'},okargs{k}, mfilename);
                case 3  % FROMVIEWER
                    if isempty(bioinfoprivate.opttf(pval))
                        error(message('bioinfo:seqviewer:FromViewerOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                    from_viewer = bioinfoprivate.opttf(pval,okargs{k}, mfilename);
            end
        end
    end
    
    if nargin == 2
        varname = varargin{1};
    else
        varname = inputname(1);
    end
    
    if(isempty(varname)) % get the default variable name
        varname = varn;
    end
    
    if isstruct(seq) % is a struct
        N = numel(seq);
        names = varname;
        
        for i = 1:N
            if(N > 1)
                names = sprintf('%s (%d)', varname, i);
            end
            handle_struct(seq(i), names, alpha);
        end
    elseif(ischar(seq) && size(seq,1)==1 && (size(fullfile(pwd,seq), 2) < 1000) && (exist(seq,'file') || exist(fullfile(pwd,seq),'file')))
        handle_file(seq, varname, alpha);
    elseif((ischar(seq) || (isnumeric(seq) && isreal(seq)))&& isvector(seq)) % char array seq or accession
        handle_chararray(seq, varname, alpha, download);
    else
        error(message('bioinfo:seqviewer:UnKnownSequencetype'))
    end
catch ME
    msg = ME.message;
    msgid = ME.identifier;
    if ~from_viewer
        if(~seqviewer.SVUtil.svMatlab('manage_bioinfodesktop', 'state')) % no viewer present
            rethrow(ME);
        else % show error in viewer
            seqviewer.SVUtil.svMatlab('create_ui', msg);
        end
    end
end

if nargout > 0
    varargout{1} = msg;
end

if nargout == 2
    varargout{2} = msgid;
end

end



function open_seqview(s)
% Open s in the viewer
% s - structure to be passed to viewer
    seqviewer.SVUtil.svMatlab('create_ui',s);
    h = seqviewer.SVUtil.svMatlab('get_seqviewer');
    st = h.sequenceIsOpen('SequenceTool');
    if(st)
        h.sequenceClose('SequenceTool');
    end
end


function handle_struct(seq, varname, alpha)
% Handle the input structure case, and open the sequence in the struct
% seq       - input sequence structure
% varname   - variable name to create sequence ids.
% alpha     - input sequence alphabet type

fnames = fieldnames(seq);
%mandatory genbank fields
gbankfields = {'LocusName','Definition','Accession','Version','GI',...
    'LocusSequenceLength','LocusNumberofStrands','LocusTopology', ...
    'LocusMoleculeType','LocusGenBankDivision','LocusModificationDate'};
%mandatory embl fields
emblfields ={'Identification','Accession','SequenceVersion','DateCreated', ...
    'Description','Keyword','OrganismSpecies','OrganismClassification', ...
    'Reference','Comments', 'Feature','BaseCount', 'Sequence'};
%struct with Header field
headerfields = {'Header'};
%check for Sequence
if(~strcmp(fnames,'Sequence'))
    error(message('bioinfo:seqviewer:NoSequenceFound'));
end
if(all(ismember(gbankfields,fnames)))%check for genbank/genpept
    if(strcmp(guess_seq_type(seq.Sequence),'NT'))
        tempacc = strread([seq.Accession(1:3) regexprep(seq.Accession(4:end), '_','')],'%s','delimiter',' ');
        if ~isfield(seq,'PreambleText')
            try
                pretemp = getgenbank(tempacc{1},'preambletext',true);
            catch allExceptions %#ok<NASGU>
                pretemp = seq;
                pretemp.PreambleText = '';
            end
        else
            pretemp = seq;
        end
        s = build_genbank(pretemp,'genbank',seq);
    else
        tempacc = strread([seq.Accession(1:3) regexprep(seq.Accession(4:end), '_','')],'%s','delimiter',' ');
        if ~isfield(seq,'PreambleText')
            try
                pretemp = getgenpept(tempacc{1},'preambletext',true);
            catch allExceptions %#ok<NASGU>
                pretemp = seq;
                pretemp.PreambleText = '';
            end
        else
            pretemp = seq;
        end
        s = build_genpept(pretemp,'genbank',seq);
    end
elseif(all(ismember(emblfields,fnames))) %check for embl
    temp = getembl(seq.Accession,'preambletext',true);
    s = build_embl(temp,'embl',seq);
else
    if ismember(headerfields,fnames)
        varname = seq.Header;
    end
    
    %find out what sequence type is
    if(isempty(alpha))
        alpha = guess_seq_type(seq.Sequence);
    end
    
    s = build_otherstruct(seq,varname, alpha, 'otherstruct',seq);
end

open_seqview(s)

end

function handle_file(seq, varname, alpha)
%% Handle the input is filename case, and open the sequence in the file
% seq       - input filename
% varname   - variable name to create sequence ids.
% alpha     - input sequence alphabet type

try
    file_path = seq;
    [data, fileformat] = bioinfoprivate.seqread(seq, 'preambletext', true);
    if(exist(fullfile(pwd,seq),'file'))
        file_path = fullfile(pwd,seq);
        seq = file_path;
    end

    N = numel(data);

    for i = 1:N
        if(N > 1)
            seq = sprintf('%s (%d)', file_path, i);
        end

        switch fileformat
            case {'GENBANK'}
                s = build_genbank(data(i),'genbankfile',seq);
            case {'GENPEPT'}
                s = build_genpept(data(i), 'genpeptfile',seq);
            case {'FASTA', 'TEXT'}
                if(isempty(alpha))
                    alpha = guess_seq_type(data(i).Sequence);
                end
                s = build_fasta(data(i),varname, alpha,'fastafile',seq);
            case {'EMBL'}
                s = build_embl(data(i),'emblfile',seq);
            otherwise
                error(message('bioinfo:seqviewer:UnKnownFiletype', seq))
        end

        open_seqview(s);
    end
catch ME1
    if strfind(ME1.identifier,'seqread')
        bioinfoprivate.bioerrorrethrow('seqviewer',ME1)
    else
        rethrow(ME1)
    end
end
end

function handle_chararray(seq, varname, alpha, download)
% Handle the input char array case, and open the sequence, or is accession
% number called from the viewer
%
% seq           - input sequence, or accession number
% varname       - variable name to create sequence ids.
% alpha         - input sequence alphabet type
% download      - should download the sequence from NCBI of the given
%                 accession number

if size(seq,2) == 1
    seq = seq';
end
if(isempty(download))
    if(~isempty(strfind(seq,'.')) && length(seq)<=32) % if there is a period and it is short, guess that this is a file
        error(message('bioinfo:seqviewer:UnknownFile', seq))
    end
    if(isempty(alpha))
        alpha = guess_seq_type(seq);
    end
    if(isnumeric(seq) && strcmp(alpha,'AMINO'))
        tempseq  = int2aa(seq);
        s = build_char(tempseq, varname, alpha, 'char',tempseq);
    elseif(isnumeric(seq) && strcmp(alpha,'NT'))
        tempseq  = int2nt(seq);
        s = build_char(tempseq, varname, alpha, 'char',tempseq);
    else
        tempseq = seq;
        s = build_char(tempseq, varname, alpha, 'char',seq);
    end
elseif(strcmpi(download,'genbank'))
    try
        acc = getgenbank(seq,'preambletext',true);
        s = build_genbank(acc,  'genbank',acc);
        s.OriginalData = rmfield(s.OriginalData,'PreambleText');
    catch theException
        % if there is an underscore then try without it.
        if ~isempty(strfind(seq,'_'))
            try
                tempacc = [seq(1:3) strrep(seq(4:end), '_','')];
                acc = getgenbank(tempacc,'preambletext',true);
                s = build_genbank(acc,  'genbank',acc);
                s.OriginalData = rmfield(s.OriginalData,'PreambleText');
            catch allExceptions
                rethrow(theException)
            end
        else
            rethrow(theException)
        end
    end
elseif(strcmpi(download,'genpept'))
    try
        acc = getgenpept(seq,'preambletext',true);
        s = build_genpept(acc,  'genpept',acc);
        s.OriginalData = rmfield(s.OriginalData,'PreambleText');
    catch theException
        % if there is an underscore then try without it.
        if ~isempty(strfind(seq,'_'))
            try
                tempacc = [seq(1:3) strrep(seq(4:end), '_','')];
                acc = getgenpept(tempacc,'preambletext',true);
                s = build_genpept(acc,  'genpept',acc);
                s.OriginalData = rmfield(s.OriginalData,'PreambleText');
            catch allExceptions
                rethrow(theException)
            end
        else
            rethrow(theException)
        end
    end
elseif(strcmpi(download,'embl'))
    acc = getembl(seq,'preambletext',true);
    s = build_embl(acc,'embl',acc);
    s.OriginalData = rmfield(s.OriginalData,'PreambleText');
end
open_seqview(s);
end

function s = build_otherstruct(acc,name,alpha,type,orig)  %#ok<INUSL>
s.FileType = type;
s.LocusName = name;
s.Accession = '';
s.Definition = name;
s.Features = '';
s.Sequence = orig.Sequence;
s.OriginalData  = orig;
s.Alpha = alpha;
s.PreambleText = '';
end

function s = build_char(acc,name,alpha,type,orig) %#ok
s.FileType = type;
s.LocusName = name;
s.Accession = '';
s.Definition = name;
s.Features = '';
s.Sequence = orig;
s.OriginalData  = orig;
s.Alpha = alpha;
s.PreambleText = '';
end

function s = build_genbank(acc,type,orig)
s.FileType = type;
s.LocusName = acc.LocusName;
s.Accession = acc.Accession;
s.Definition = acc.Definition;
s.Features = acc.Features;
if(strcmp(type,'genbankfile'))
    s.Sequence = acc.Sequence;
else
    s.Sequence = orig.Sequence;
end
s.OriginalData  = orig;
s.Alpha = 'NT';
s.PreambleText = acc.PreambleText;
end

function s = build_genpept(acc, type,orig)
s.FileType = type;
s.LocusName = acc.LocusName;
s.Accession = acc.Accession;
s.Definition = acc.Definition;
s.Features = acc.Features;
if(strcmp(type,'genpeptfile'))
    s.Sequence = acc.Sequence;
else
    s.Sequence = orig.Sequence;
end
s.OriginalData  = orig;
s.Alpha = 'AMINO';
s.PreambleText = acc.PreambleText;
end

function s = build_fasta(acc, name, alpha, type,orig) %#ok
s.FileType = type;
s.LocusName = acc(1).Header;
s.Accession = '';
s.Definition = acc(1).Header;
s.Features = '';
s.Sequence = acc(1).Sequence;
s.OriginalData  = orig;
s.Alpha = alpha;
s.PreambleText = acc(1).Header;
end

function s = build_embl(acc, type, orig)
s.FileType = type;
s.LocusName = acc.Accession;
s.Accession = acc.Accession;
s.Definition = acc.Description;
s.Features = acc.Feature;
if(strcmp(type,'emblfile'))
    s.Sequence = acc.Sequence;
else
    s.Sequence = orig.Sequence;
end
s.OriginalData  = orig;
s.Alpha = 'NT';
s.PreambleText = acc.PreambleText;
end

function s = build_empty(acc)
s.FileType = 'empty';
s.LocusName = acc;
s.Accession = '';
s.Definition = acc;
s.Features = '';
s.Sequence = '';
s.OriginalData  = '';
s.Alpha = 'AMINO';
s.PreambleText ='';
end

function t = guess_seq_type(seq)
%find out what sequence type is
if(~bioinfoprivate.isnt(seq))
    t = 'AMINO';
else
    t = 'NT';
end
end

function result = jvm_available
persistent sRlt;
if isempty(sRlt)
    sRlt=usejava('jvm') & usejava('awt') & usejava('swing');
end
result = sRlt;
end
