function tree=gethmmtree(accessnum,varargin)
%GETHMMTREE retrieves phylogenetic tress from the PFAM database.
%
%   TREE = GETHMMTREE(KEY) searches for the PFAM family accession KEY in
%   the PFAM database and returns an object containing a phylogenetic tree
%   representative of such protein family. KEY must be a valid string to
%   query the PFAM database, some databases use the PFAM version number
%   appended to the accession KEY.
%
%   TREE = GETHMMTREE(NUM) uses the numeric value in NUM to figure out the
%   appropriate accession KEY to query the database and retrieve the
%   information into a PHYTREE object. For example; if NUM = 2, GETHMMTREE
%   will retrieve the tree for the family 'PF00002'. 
%
%   TREE = GETHMMTREE(...,'TOFILE',FILENAME) saves the data returned from
%   the PFAM database in the file FILENAME.
%
%   TREE = GETHMMTREE(...,'TYPE',TYPE) returns only the tree of the
%   alignments used to generate the HMM model when TYPE is 'seed', and
%   returns the tree of all alignments that hit the model if TYPE is
%   'full'. The default TYPE is 'full'. 
% 
%   Examples:
%
%       tree  = gethmmtree(2)
%       
%       tree  = gethmmtree('PF00002')
%
%       These both retrieve the phylogenetic tree built from the multiple
%       aligned sequences used to train the HMM profile model for global
%       alignment to the 7-transmembrane receptor, secretin family (pfam
%       key = PF00002).
%
%   See also GETHMMALIGNMENT, GETHMMPROF, PHYTREEREAD.

% Copyright 2003-2005 The MathWorks, Inc.



if ~usejava('jvm')
    error(message('bioinfo:gethmmtree:NeedJVM', mfilename));
end

if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:gethmmtree:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'type','tofile'};
    for j=1:2:nargin-2
        pname = varargin{j};
      %  pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:gethmmtree:UnknownParameter', pname));
        elseif length(k)>1
            error(message('bioinfo:gethmmtree:AmbiguousParameter', pname)) ;
        end
    end
end

tree = getsangerdata(accessnum,'database','tree',varargin{:});

