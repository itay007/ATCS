function blastformat(varargin)
%BLASTFORMAT create local BLAST database
%
%   BLASTFORMAT calls a local version of the NCBI FORMATDB executable file
%   and formats the input sequences as a local, blastable database
%   according to the commands specified in the input parameters. These can
%   be given in parameter/value pairs, in a string containing valid
%   FORMATDB commands or in a combination of the two. 
%
%   BLASTFORMAT(..., 'INPUTDB', DBFILENAME) specifies the name of the fasta
%   file containing the set of sequences to be formatted as blastable
%   database. Corresponds to the FORMATDB option -i.
%
%   BLASTFORMAT(..., 'FORMATPATH', FULLPATH) specifies the full path to the
%   FORMATDB executable file, including the name and extension of
%   the executable file. Default is the system path.
%
%   BLASTFORMAT(..., 'TITLE', TITLENAME) specifies the title for the local
%   database. Corresponds to the FORMATDB option -t. Default is the input
%   file name.
%
%   BLASTFORMAT(..., 'LOG', LOGNAME) specifies the name of the log file.
%   Corresponds to the FORMATDB option -l. Default is formatdb.log.
%
%   BLASTFORMAT(..., 'PROTEIN', TF) specifies whether the sequences to be
%   formatted are proteins or not. Corresponds to FORMATDB option -p.
%   Default is true.
%
%   BLASTFORMAT(..., 'FORMATARGS', FARGS) specifies a NCBI FORMATDB command
%   string FARGS, containing one or more instances of -x and the option
%   value associated with it. See FORMATDB help for a complete list of the
%   options available.
%
%   Examples:
%   The FASTA files with the E. coli genome used in the following examples
%   can be downloaded from:
%   ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_536_uid58531/
% 
%   Note: For easy identification, please rename the downloaded FASTA files
%   NC_008253.fna to ecoli.nt and NC_008253.faa to ecoli.aa.
%   
%   % Example 1
%   % Create a local blastable database and give it a title.
%   blastformat('inputdb', 'ecoli.nt', 'protein', 'false', 'title', 'myecoli_nt');
%   
%   % Example 2
%   % Create a local blastable database and specify a title and a log file.
%   blastformat('inputdb', 'ecoli.aa', 'formatargs', '-t myecoli_aa -l ecoli_aa.log');
%
%   % Example 3
%   % Create a local blastable database using the FORMATDB command syntax
%   % (FORMATDB must be in your system path)
%   blastformat('-i ecoli.nt -p F -l ecoli.log');
%  
% See also BLASTLOCAL, BLASTNCBI, BLASTREAD, BLASTREADLOCAL, GETBLAST.

% References
% http://www.ncbi.nlm.nih.gov/blast/docs/formatdb.html
% ftp://ftp.ncbi.nih.gov/genomes/Bacteria/

% Copyright 2007-2012 The MathWorks, Inc.


formatPath = []; % if not set, assume the executable is in path
cmd = [];        % blast command string

if nargin == 0
    error(message('bioinfo:blastformat:NotEnoughInput'));

elseif nargin > 1
    if rem(nargin,2) == 1
        error(message('bioinfo:blastformat:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'inputdb', 'title', 'log', 'protein', 'formatargs', 'formatpath'};
    
    for j=1:2:nargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs, numel(pname)));
        if isempty(k)
            error(message('bioinfo:blastformat:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:blastformat:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % inputdb
                    if ischar(pval) && exist(pval, 'file')
                        cmd = [cmd ' -i ' pval]; %#ok<AGROW>
                    elseif isstruct(pval) % not documented
                        tfilename = tempname;
                        fastawrite(tfilename, pval);
                        cmd = [cmd ' -i ' tfilename]; %#ok<AGROW>
                    else
                         error(message('bioinfo:blastformat:InvalidInput')); 
                    end
                case 2  % title
                    if ischar(pval)
                        cmd = [cmd ' -t ' pval]; %#ok<AGROW>
                    else
                        error(message('bioinfo:blastformat:InvalidTitle')); 
                    end
                case 3  % log
                    if ischar(pval)
                        cmd = [cmd ' -l ' pval]; %#ok<AGROW>
                    else
                        error(message('bioinfo:blastformat:InvalidLogName')); 
                    end
                case 4  % protein
                    prot = bioinfoprivate.opttf(pval);
                    if isempty(prot)
                        error(message('bioinfo:blastformat:InvalidProteinChoice')); 
                    elseif ~prot
                        cmd = [cmd ' -p F']; %#ok<AGROW>
                    end 
                case 5  % formatargs
                    if ischar(pval)
                        opt = regexprep(pval, '[\.|:|/|\\]', ''); %ignore these symbols in file paths when counting words
                        L = length(regexp(opt, '\w+')); 
                        pairs = length(regexp(opt, '\s*-\w\s+\w+\s*')); % number of options set properly 
                        if rem(L,2) == 0 && pairs == L/2
                            cmd = [cmd ' ' pval]; %#ok<AGROW>
                        else
                            error(message('bioinfo:blastformat:InvalidCommandString', pval));
                        end
                    else
                         error(message('bioinfo:blastformat:InvalidFormatcmdOption'));
                    end
                case 6 % formatpath
                    if exist(pval, 'dir') % does not include the executable name
                        error(message('bioinfo:blastformat:IncompletePath'));
                    elseif exist(pval, 'file') 
                        [xdir,xfile] = fileparts(pval);
                        formatPath = fullfile(xdir,xfile);
                    else
                        error(message('bioinfo:blastformat:InvalidPath'));
                    end
            end
        end
    end

else % nargin == 1
    cmd = [' ' varargin{:}];
end

%=== run FORMATDB 
if isempty(formatPath)
        [formatStatus, formatResult]= system(['formatdb ' cmd]);
else
    % Quote formatPath as it might contain spaces
    [formatStatus, formatResult] = system(['"' formatPath '"' cmd]);
end

if formatStatus ~= 0 
     error(message('bioinfo:blastformat:formatdbError', formatResult));
end

