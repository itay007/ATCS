function pdbstruct=getpdb(pdbID,varargin)
%GETPDB retrieves sequence information from the Protein Data Bank.
%
%   PDBSTRUCT = GETPDB(PDBID) searches for the ID in the Protein Data Bank
%   (PDB) database and returns a structure containing information for the
%   protein.
%
%   PDBSTRUCT = GETPDB(...,'TOFILE',FILENAME) saves the data returned from
%   the database in the file FILENAME.
%
%   PDBSTRUCT = GETPDB(...,'SEQUENCEONLY',true) returns just the protein
%   sequence. If the PDB file contains only one sequence then this will be
%   returned as a character array. If more than one sequence is found, then
%   these will be returned in a cell array.
%
%   Example:
%
%       pdbstruct = getpdb('2DHB')
%
%   This retrieves the structure information for horse deoxyhemoglobin
%   (PDB ID 2DHB).
%
%   See also GETEMBL, GETGENBANK, GETGENPEPT, MOLVIEWER, PDBDISTPLOT,
%   PDBREAD.

%   Reference:
%   H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H.
%   Weissig, I.N. Shindyalov, P.E. Bourne: The Protein Data Bank. Nucleic
%   Acids Research, 28 pp. 235-242 (2000)

%   Copyright 2002-2012 The MathWorks, Inc.
%

%   Mirror sites are no longer supported after Jan 1st 2006.
%   PDBSTRUCT = GETPDB(...,'MIRRORSITE',MIRROR) allows you to choose a
%   mirror site for the PDB database. The default site is the San Diego
%   Supercomputer Center, http://www.rcsb.org/pdb. Set MIRROR to
%   http://rutgers.rcsb.org/pdb to use the Rutgers University Site or
%   http://nist.rcsb.org/pdb for the National Institute of Standards and
%   Technology site. Follow this link for a full list of PDB mirror sites:
%   http://www.rcsb.org/pdb/static.do?p=general_information/mirror_sites/in
%   dex.html

if ~usejava('jvm')
    error(message('bioinfo:getpdb:NeedJVM', mfilename));
end

tofile = false;
seqonly = false;
PDBsite = 'http://www.rcsb.org/pdb';

if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:getpdb:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'tofile','mirrorsite','sequenceonly'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:getpdb:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:getpdb:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1    % tofile
                    if ischar(pval)
                        tofile = true;
                        filename = pval;
                    end
                case 2    % mirrorsite
                    warning(message('bioinfo:getpdb:MirrorSiteNotSupported'))
                    
                case 3  % sequenceonly
                    seqonly = bioinfoprivate.opttf(pval);
                    if isempty(seqonly)
                        error(message('bioinfo:getpdb:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
            end
        end
    end
end


% error if ID isn't a string
if ~ischar(pdbID)
    error(message('bioinfo:getpdb:NotString'))
end

% get sequence from pdb.fasta if SEQUENCEONLY is true, otherwise full pdb
if seqonly == true
    searchurl = [PDBsite '/downloadFile.do?fileFormat=FASTA&compression=NO&structureId=' pdbID];
    [~, pdb] = fastaread(searchurl);
else
    searchurl = [PDBsite '/downloadFile.do?fileFormat=pdb&compression=NO&structureId=' pdbID];
    
    % get the html file that is returned as a string
    s=urlread(searchurl);
    
    % replace the html version of &
    s=strrep(s,'&amp;','&');
    
    % Find first line of the actual data
    start = strfind(s,'HEADER');
    
    if isempty(start)
        % search for text indicating that there weren't any files found
        notfound=regexp(s,'The file you requested does not exist.','once');
        
        % string was found, meaning no results were found
        if ~isempty(notfound),
            error(message('bioinfo:getpdb:PDBIDNotFound', pdbID)) ;
        end
        error(message('bioinfo:getpdb:PDBIDAccessProblem', pdbID));
    end
    
    [~, endOfFile] = regexp(s,'\nEND\s.*\n');
    
    % shorten string, to search for uid info
    s=s(start:endOfFile);
    
    %make each line a separate row in string array
    pdbdata = char(strread(s,'%s','delimiter','\n','whitespace',''));
    
    %pass to PDBREAD to create structure
    pdb=pdbread(pdbdata);
    
    
end

if nargout
    pdbstruct = pdb;
    if ~seqonly
        % add URL
        pdbstruct.SearchURL = searchurl;
    end
else
    if seqonly || ~usejava('desktop')
        disp(pdb);
    else
        disp(pdb);
        disp([char(9) 'SearchURL: <a href="' searchurl '"> ' pdbID ' </a>']);
    end
    
end

%  write out file
if tofile == true
    fid = fopen(filename,'wt') ;
    if fid == (-1)
        error(message('bioinfo:getpdb:CouldNotOpenFile', filename));
    else
        rows = size(pdbdata,1);
        
        for rcount=1:rows-1,
            fprintf(fid,'%s\n',pdbdata(rcount,:));
        end
        fprintf(fid,'%s',pdbdata(rows,:));
        fclose(fid);
    end
end
end
