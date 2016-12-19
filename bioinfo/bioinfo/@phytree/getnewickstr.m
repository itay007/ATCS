function str = getnewickstr(tr,varargin)
%GETNEWICKSTR creates a NEWICK formatted string.
%
%   STR = GETNEWICKSTR(TREE) returns the NEWICK formatted string of the
%   phylogenetic tree object TREE. 
%
%   GETNEWICKSTR(...,'DISTANCES',false) excludes the distances from the
%   output. Default is true.
%
%   GETNEWICKSTR(...,'BRANCHNAMES',true) includes the branch names into the
%   output. Default is false. 
%
%   The NEWICK tree format is found at: 
%        http://evolution.genetics.washington.edu/phylip/newicktree.html
%
%   Example:
%
%      seqs = int2nt(ceil(rand(10)*4));      % some random sequences
%      dist = seqpdist(seqs,'alpha','nt');   % pairwise distances
%      tree = seqlinkage(dist);              % construct phylogenetic tree
%      str  = getnewickstr(tree)             % get the NEWICK string
%     
%   See also PHYTREE, PHYTREEREAD, PHYTREEWRITE, PHYTREEVIEWER, SEQLINKAGE,
%   PHYTREE/GET, PHYTREE/GETBYNAME, PHYTREE/GETCANONICAL. 

%  Undocumented:
%   GETSTR(...,'MULTILINE',true) introduces 'new line' characters for a
%   multi-line output. This option is used by PHYTREEWRITE. Default is
%   false.

% Copyright 2003-2012 The MathWorks, Inc.


if numel(tr)~=1
     error(message('bioinfo:phytree:getstr:NoMultielementArrays'));
end

numBranches = size(tr.tree,1);
numLeaves = numBranches + 1;
numLabels = numBranches + numLeaves; 

% set defaults
writeDistances = true;
writeBranchNames = false;
multiLine = false;

nvarargin = numel(varargin);
if nvarargin
    if rem(nvarargin,2)
        error(message('bioinfo:phytree:getstr:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'multiline','distances','branchnames'};
    for j=1:2:nvarargin
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:phytree:getstr:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:phytree:getstr:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % multi-lines
                    multiLine = bioinfoprivate.opttf(pval);
                    if isempty(multiLine)
                        error(message('bioinfo:phytree:getstr:multiLineOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 2 % write distances
                    writeDistances = bioinfoprivate.opttf(pval);
                    if isempty(writeDistances)
                        error(message('bioinfo:phytree:getstr:writeDistancesOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 3 % write branch names
                    writeBranchNames = bioinfoprivate.opttf(pval);
                    if isempty(writeBranchNames)
                        error(message('bioinfo:phytree:getstr:writeBranchNamesOptionNotLogical', upper( char( okargs( k ) ) )));
                    end                    
            end
        end
    end
end

% add escape characters for special characters valid in the Newick format
tr.names = regexprep(tr.names,'[\(\),;:]','\\\\$0');

for i=1:numLabels-1;
    if (i<=numLeaves || writeBranchNames)
        if writeDistances 
            namedist{i} = [tr.names{i} ':' num2str(tr.dist(i))]; %#ok
        else
            namedist{i} = tr.names{i};
        end
    elseif writeDistances && ~writeBranchNames
        namedist{i} = [':' num2str(tr.dist(i))]; %#ok
    else 
        namedist{i} = '';
    end
end

if writeBranchNames
    namedist{numLabels} = [tr.names{numLabels} ';'];
else
    namedist{numLabels} = ';';
end

for i=1:numBranches
    if tr.tree(i,1) > numLeaves
        t1 = branchstr{tr.tree(i,1)};
    else
        t1 = namedist{tr.tree(i,1)};
    end
    if tr.tree(i,2) > numLeaves
        t2 = branchstr{tr.tree(i,2)};
    else
        t2 = namedist{tr.tree(i,2)};
    end
    branchstr{i+numLeaves} = ...
        [ '(\n' t1 ',\n' t2 ')\n' , namedist{i+numLeaves} ]; %#ok
end

str = sprintf(branchstr{numLabels});

if ~multiLine
    str = strrep(str,sprintf('\n'),'');
end

