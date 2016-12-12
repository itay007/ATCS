function phytreewrite(filename,varargin)
%PHYTREEWRITE saves a phylogenetic tree object as a NEWICK format file.
%
%   PHYTREEWRITE(FILENAME, TREE) writes the contents of a phylogenetic tree
%   TREE to file FILENAME using the NEWICK format. TREE must be an PHYTREE
%   object.
%
%   PHYTREEWRITE(TREE) writes the contents of a phylogenetic tree TREE to
%   file 'output.tree' in the current directory.
%
%   PHYTREEWRITE(...,'DISTANCES',false) excludes the distances from the
%   output. Default is true.
%
%   PHYTREEWRITE(...,'BRANCHNAMES',false) excludes the branch names from
%   the output. Default is true.
%
%   Example:
%
%      seqs = int2nt(ceil(rand(10)*4));      % some random sequences
%      dist = seqpdist(seqs,'alpha','nt');   % pairwise distances
%      tree = seqlinkage(dist);              % construct phylogenetic tree
%      phytreewrite('mytree.tree',tree)      % write it to NEWICK file
%
%   See also PHYTREE, PHYTREE/GETNEWICKSTR, PHYTREEREAD, PHYTREEVIEWER,
%   SEQLINKAGE.

%  Undocumented:
%   PHYTREEWRITE(...,'GUI',true) opens a standard get file dialog box to
%   select the target filename. FILENAME (if provided) is ignored. Default
%   is false. 

% Copyright 2003-2012 The MathWorks, Inc.


enableGui = false;
writeDistances = true;
writeBranchNames = true;

if isa(filename,'phytree')
    tr = filename; 
    filename = 'output.tree';
else
    if ~ischar(filename)
        error(message('bioinfo:phytreewrite:FilenameMustBeString'));
    end
    tr = varargin{1};
    varargin(1) = [];
end    

% check input object
if ~isa(tr,'phytree')
        error(message('bioinfo:phytreewrite:IncorrectInputType'));
end
if numel(tr)~=1
     error(message('bioinfo:phytreewrite:NoMultielementArrays'));
end

nvarargin = numel(varargin);
if  nvarargin 
    if rem(nvarargin,2)
        error(message('bioinfo:phytreewrite:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'gui','distances','branchnames'};
    for j=1:2:nvarargin
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:phytreewrite:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:phytreewrite:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % gui
                    enableGui = bioinfoprivate.opttf(pval);
                    if isempty(enableGui)
                        error(message('bioinfo:phytreewrite:enableGuiOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 2 % write distances
                    writeDistances = bioinfoprivate.opttf(pval);
                    if isempty(writeDistances)
                        error(message('bioinfo:phytreewrite:writeDistancesOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 3 % write branch names
                    writeBranchNames = bioinfoprivate.opttf(pval);
                    if isempty(writeBranchNames)
                        error(message('bioinfo:phytreewrite:writeBranchNamesOptionNotLogical', upper( char( okargs( k ) ) )));
                    end                    
            end
        end
    end
end

if enableGui
    [filename, pathname] = uiputfile('*.tree','Save Phylogenetic tree as');
    if ~filename
        return;
    end
    filename = [pathname, filename];
end

sout = getnewickstr(tr,'dis',writeDistances,'bran',writeBranchNames,'multi',1);

fid = fopen(filename,'wt') ;
if fid == (-1)
    if enableGui
        errordlg(sprintf('Could not open file %s.', filename))
    else
        error(message('bioinfo:phytreewrite:CouldNotOpenFile', filename));
    end
else
    fprintf(fid,'%s',sout) ;
    fclose(fid) ;
end
