function model=gethmmprof(accessnum,varargin)
%GETHMMPROF retrieves profile Hidden Markov Models from the PFAM database.
%
%   MODEL = GETHMMPROF(KEY) searches for the PFAM family accession KEY in
%   the PFAM database and returns a structure containing information for
%   the HMM profile. KEY must be a valid string to query the PFAM database,
%   some databases use the PFAM version number appended to the accession
%   KEY.
%
%   MODEL = GETHMMPROF(NUM) uses the numeric value in NUM to figure out the
%   appropriate accession KEY to query the database and retrieve the
%   information into a structure. For example; if NUM = 2, GETHMMPROF will
%   retrieve the model for the family 'PF00002'. 
%
%   MODEL = GETHMMPROF(...,'TOFILE',FILENAME) saves the data returned from
%   the PFAM database in the file FILENAME.
%
%   MODEL = GETHMMPROF(...,'MODE',MODE) returns the global alignment model
%   when MODE='ls' and the local alignment model when MODE='fs'. Default is
%   'ls'.
%
%   MODEL = GETHMMPROF(...,'MIRROR',MIRROR) selects a specific web
%   database. Options are 'Sanger' (default) or 'Janelia'. Other mirror
%   sites can be reached by passing the complete URL to PFAMHMMREAD. Note:
%   these mirrors are maintained separately, therefore slight variations
%   may be found. For more information:
%   http://pfam.sanger.ac.uk/ and http://pfam.janelia.org/.
%   
%
%   Examples:
%
%       hmmmodel  = gethmmprof(2)
%       
%       hmmmodel  = gethmmprof('7tm_2')
%
%       % These both retrieve an HMM profile model for global alignment to
%       % the 7-transmembrane receptor, secretin family (pfam key = PF00002).
%
%   See also GETHMMALIGNMENT, GETHMMTREE, HMMPROFALIGN, HMMPROFDEMO,
%   HMMPROFSTRUCT, PFAMHMMREAD, SHOWHMMPROF.

% Copyright 2003-2008 The MathWorks, Inc.



if ~usejava('jvm')
    error(message('bioinfo:gethmmprof:NeedJVM', mfilename));
end

mirrorOptions = {'sanger','janelia'};
selectedMirror = 1; %default option

% find if one of the input args is a mirror selector
if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:gethmmprof:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'database','type','mode','tofile','mirror'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:gethmmprof:UnknownParameter', pname));
        elseif length(k)>1
            error(message('bioinfo:gethmmprof:AmbiguousParameter', pname)) ;
        elseif k==5 % mirror
            selectedMirror = strmatch(lower(pval),mirrorOptions);
            if isempty(selectedMirror)
                error(message('bioinfo:gethmmprof:InvalidMirror'))
            end
        end
    end
end

switch selectedMirror
    case 1
        % call getsangerdata with database set to 'hmm'
        model = getsangerdata(accessnum,'database','hmm',varargin{:});
    case 2
        % call getjaneliapfamdata with database set to 'hmm'
        model = getjaneliapfamdata(accessnum,'database','hmm',varargin{:});
end
