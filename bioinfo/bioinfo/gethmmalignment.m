function align=gethmmalignment(accessnum,varargin)
%GETHMMALIGNMENT retrieves multiple sequences aligned to a PFAM profile HMM.
%
%   SEQS = GETHMMALIGNMENT(KEY) returns a structure containing information
%   for a multiple alignment of sequences to the profile Hidden Markov Model
%   KEY in the PFAM database. KEY must be a valid string to query the PFAM
%   database, some databases use the PFAM version number appended to the
%   accession KEY.
%
%   SEQS = GETHMMALIGNMENT(NUM) uses the numeric value in NUM to figure out
%   the appropriate accession KEY to query the database and retrieve the
%   information into a structure. For example; if NUM = 2, GETHMMALIGNMENT
%   will retrieve the multiple alignment for the family 'PF00002'. 
%
%   SEQS = GETHMMALIGNMENT(..., 'TYPE', TYPE) returns only the alignments
%   used to generate the HMM model when TYPE is 'seed', and returns all
%   alignments that hit the model if TYPE is 'full'. The default TYPE is
%   'full'. 
%
%   SEQS = GETHMMALIGNMENT(..., 'MIRROR', MIRROR) selects a specific web
%   database. Options are 'Sanger' (default) or 'Janelia'. Other mirror
%   sites can be reached by passing the complete URL to FASTAREAD. Note:
%   these mirrors are maintained separately, therefore slight variations
%   may be found. For more information:
%   http://pfam.sanger.ac.uk/ and http://pfam.janelia.org/.
%
%   SEQS = GETHMMALIGNMENT(...,'IGNOREGAPS',true) removes any gap symbol
%   ('-' or '.') from the sequences. Default is false.
%
%   Examples:
%       % These both retrieve a multiple alignment of the sequences used to
%       % train the HMM profile model for global alignment to the
%       % 7-transmembrane receptor, secretin family (pfam key = PF00002).
%       pfamalign = gethmmalignment(2, 'type', 'seed')
%       
%       pfamalign = gethmmalignment('PF00002', 'type', 'seed')
%       
%       % View the alignment in the multiple alignment viewer
%       multialignviewer(pfamalign)
%
%   See also FASTAREAD, GETHMMPROF, GETHMMTREE, HMMPROFDEMO,
%   MULTIALIGNREAD, SEQALIGNVIEWER, PFAMHMMREAD.

% Copyright 2003-2012 The MathWorks, Inc.



if ~usejava('jvm')
    error(message('bioinfo:gethmmalignment:NeedJVM', mfilename));
end


mirrorOptions = {'sanger','janelia'};
selectedMirror = 1; %default option
ignoreGaps = false; %#ok<NASGU> varargs are passed through to private function

% find if one of the input args is a mirror selector
if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:gethmmalignment:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'database','type','mode','tofile','ignoregaps','mirror'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:gethmmalignment:UnknownParameter', pname));
        elseif length(k)>1
            error(message('bioinfo:gethmmalignment:AmbiguousParameter', pname)) ;
        elseif k==6 % mirror
            selectedMirror = strmatch(lower(pval),mirrorOptions); 
            if isempty(selectedMirror)
                error(message('bioinfo:gethmmalignment:InvalidMirror'))
            end
        end
    end
end

switch selectedMirror
    case 1
        % call getsangerdata with database set to 'align'
        align = getsangerdata(accessnum,'database','align',varargin{:});
    case 2
        % call getjaneliapfamdata with database set to 'align'
        align = getjaneliapfamdata(accessnum,'database','align',varargin{:});
end
