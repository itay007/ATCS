function results = blastlocal(varargin)
%BLASTLOCAL perform search on local BLAST database to create BLAST report
%
%   BLASTLOCAL calls a local version of the NCBI BLASTALL executable file
%   and searches the input sequence(s) against a local, blastable database
%   according to the commands specified in the input parameters. These can
%   be given in parameter/value pairs, in a string containing valid
%   BLASTALL commands or in a combination of the two. 
% 
%   DATA = BLASTLOCAL(...) returns the BLAST search results in DATA, a
%   MATLAB structure or array of structures (if multiple query sequences),
%   containing fields corresponding to BLAST keywords.
%
%   BLASTLOCAL(..., 'INPUTQUERY', QUERY) specifies the fasta file
%   containing the query sequence(s). Corresponds to the BLASTALL option
%   -i. 
%
%   BLASTLOCAL(..., 'PROGRAM', PNAME) specifies the program to be used in
%   the search. Valid programs are: blastn, blastp, blastx, tblastn,
%   tblastx. Default is blastp. Corresponds to BLASTALL option -p.
%
%   BLASTLOCAL(..., 'DATABASE', DB) specifies the local database to search.
%   The database must have been formatted using the FORMATDB program.
%   Default is a local version of the nr database in the current directory.
%   Corresponds to the BLASTALL option -d. 
%
%   BLASTLOCAL(..., 'BLASTPATH', FULLPATH) specifies the full path to the
%   BLASTALL executable file, including the name and extension of the
%   executable file. Default is the system path.
%
%   BLASTLOCAL(..., 'EXPECT', EVALUE) specifies the statistical
%   significance threshold for matches against the database sequences.
%   Default is 10. Corresponds to the BLASTALL option -e. 
%
%   BLASTLOCAL(..., 'FORMAT', MFORMAT) specifies the alignment format of
%   the BLAST search results. Corresponds to the BLASTALL option -m. Valid
%   formats are: 
%         0 - pairwise (default)
%         1 - query-anchored, showing identities
%         2 - query-anchored, no identities
%         3 - flat query-anchored, showing identities
%         4 - flat query-anchored, no identities
%         5 - query-anchored, no identities, blunt ends
%         6 - flat query-anchored, no identities, blunt ends
%         8 - tabular
%         9 - tabular with comment lines
%
%   BLASTLOCAL(..., 'TOFILE', OUTFILE) specifies the name of the file to
%   which to save the contents of the BLAST report. Corresponds to the
%   BLASTALL option -o. 
%
%   BLASTLOCAL(..., 'FILTER', TF) specifies whether the filter DUST
%   (for blastn) or SEG (for other programs) is to be used. Default is
%   true.  Corresponds to the BLASTALL option -F. 
%
%   BLASTLOCAL(..., 'GAPOPEN', GO) specifies the cost of opening a gap.
%   Default is -1.  Corresponds to the BLASTALL option -G. 
%
%   BLASTLOCAL(..., 'GAPEXTEND', GE) specifies the cost of extending a gap.
%   Default is -1.  Corresponds to the BLASTALL option -E. 
% 
%   BLASTLOCAL(..., 'BLASTARGS', BARGS) specifies a NCBI BLASTALL command
%   string BARGS, containing one or more instances of -x and the option
%   value associated with it. See BLASTALL help for a complete list of the
%   options available.
%
%   Examples:
%      
%   % Retrieve query sequence (E. coli threonine operon) from GenBank  
%   S = getgenbank('M28570');
%
%   % Write sequence into a fasta file, using the accession number as header
%   S.Header = S.Accession;
%   fastawrite('query_nt.fa', S);
% 
%   % Please see the examples in BLASTFORMAT help on how to create the
%   % local blastable database ecoli.nt.
% 
%   % Example 1
%   % Search the sequence against the local blastable database ecoli.nt
%   rel_nt = blastlocal('inputquery', 'query_nt.fa', ...
%                       'database', 'ecoli.nt',...
%                       'program', 'blastn')
%
%   % The same search can be performed using BLASTALL syntax
%   rel_nt = blastlocal('-i query_nt.fa -d ecoli.nt -p blastn')
%
%   % Example 2 
%   % Search a nucleotide sequence against a local blastable amino acid database
%   rel_aa = blastlocal('inputquery', 'query_nt.fa', ...
%                        'database', 'ecoli.aa', ...
%                        'program', 'blastx')
%
%   % Example 3
%   % Search a nucleotide sequence and write the results into a file in tabular format
%   blastlocal('inputquery', 'query_nt.fa', ...
%              'database', 'ecoli.nt', 'tofile', 'myecoli_nt.txt', ...
%              'blastargs', '-p blastn -m 8')
%
% See also BLASTFORMAT, BLASTNCBI, BLASTREAD, BLASTREADLOCAL, GETBLAST.

% References:
% http://www.ncbi.nlm.nih.gov/blast/docs/blastall.html

% Copyright 2007-2012 The MathWorks, Inc. 



blastPath = '';     % if not set, assume the executable is in path
cmd = '';           % blast command string
mopt = 0;           % default output format
tofile = '';        % default is not to write to a file
queryFlag = false;  % =1 if query is set by user (without query BLASTALL hang)
tfilename = '';     % temporary fasta file for writing the queries
progFlag = false;   % =1 if program is set by user 

if nargin == 0
    error(message('bioinfo:blastlocal:NotEnoughInput'));

elseif nargin > 1
    if rem(nargin,2) == 1
        error(message('bioinfo:blastlocal:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'inputquery', 'program', 'database', 'expect', 'format', ...
        'tofile', 'filter', 'gapopen', 'gapextend',  'blastargs', 'blastpath'};
    okprogs = {'blastp', 'blastn', 'blastx', 'tblastx', 'tblastn'};

    for j=1:2:nargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs, numel(pname)));
        if isempty(k)
            error(message('bioinfo:blastlocal:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:blastlocal:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % input query
                    queryFlag = true;
                    if ischar(pval) &&  exist(pval, 'file')
                        cmd = [cmd ' -i ' pval]; %#ok<AGROW>
                    elseif isstruct(pval) % not documented
                        tfilename = [tempname '.fasta'];
                        fastawrite(tfilename, pval);
                        cmd = [cmd ' -i ' tfilename]; %#ok<AGROW>
                    else
                        error(message('bioinfo:blastlocal:InvalidQuery'));
                    end
                        
                case 2  % program
                    if any(strcmpi(pval, okprogs))
                        cmd = [cmd ' -p ' pval]; %#ok<AGROW>
                        progFlag = true;
                    else
                         error(message('bioinfo:blastlocal:InvalidProgram'));
                    end
                case 3  % database
                    if ischar(pval) % do not check for existence of db, since necessary files depend on program
                        cmd = [cmd ' -d ' pval]; %#ok<AGROW>
                    else
                         error(message('bioinfo:blastlocal:InvalidDatabaseName'));
                    end
                case 4  % expect
                    if isnumeric(pval) && isreal(pval)
                        cmd = [cmd ' -e ' num2str(pval)]; %#ok<AGROW>
                    else
                        error(message('bioinfo:blastlocal:InvalidEvalue'));
                    end
                case 5  % format
                    if isnumeric(pval) && isvector(pval) && (pval >= 0 && pval < 12) ...
                            && (pval - floor(pval)) == 0
                            cmd = [cmd ' -m ' num2str(pval)]; %#ok<AGROW>
                            mopt = pval; 
                    else
                        error(message('bioinfo:blastlocal:InvalidFormatOut'));
                    end
                case 6  % tofile
                    if ischar(pval)
                        cmd = [cmd ' -o ' pval]; %#ok<AGROW>
                        tofile = pval;
                    else
                        error(message('bioinfo:blastlocal:InvalidOutputName')); 
                    end
                case 7  % filter
                    filt = bioinfoprivate.opttf(pval);
                    if isempty(filt)
                        error(message('bioinfo:blastlocal:InvalidFilterChoice')); 
                    elseif ~filt
                        cmd = [cmd ' -F F']; %#ok<AGROW>
                    end
                case 8  % gapopen
                    if (pval - floor(pval)) == 0
                        cmd = [cmd ' -G ' num2str(pval)]; %#ok<AGROW>
                    else
                         error(message('bioinfo:blastlocal:InvalidGapOpen')); 
                    end
                case 9  % gapextend
                    if (pval - floor(pval)) == 0
                        cmd = [cmd ' -E ' num2str(pval)]; %#ok<AGROW>
                    else
                         error(message('bioinfo:blastlocal:InvalidGapClose')); 
                    end
                case 10 % blastargs
                    if ischar(pval)
                        
                        L = length(regexp(pval, '\S+'));
                        pairs = length(regexp(pval, '\s*-\w\s+\S+')); % number of option pairs

                        if rem(L,2) == 0 && pairs == L/2
                            cmd = [cmd ' ' pval]; %#ok<AGROW>
                            options = strread(pval, '%s');
                            optm = options(find(strcmp(options, '-m'),1)+1); % option m
                            opto = options(find(strcmp(options, '-o'),1)+1); % option o
                            optp = options(find(strcmp(options, '-p'),1)+1); % option p
                            opti = options(find(strcmp(options, '-i'),1)+1); % option i
                            
                            if ~isempty(optm)
                                mopt = str2double(optm);
                            end
                            if ~isempty(opto)
                                tofile = char(opto);
                            end
                            if ~isempty(optp)
                                progFlag = true;
                            end
                            if ~isempty(opti)
                                queryFlag = true;
                            end
                        else
                            error(message('bioinfo:blastlocal:InvalidBlastString', pval));
                        end
                    else
                         error(message('bioinfo:blastlocal:InvalidBlastcmdOption'));
                    end
                case 11 % blastPath
                    if ischar(pval)
                        if exist(pval, 'dir') % does not include the executable name
                            error(message('bioinfo:blastlocal:IncompletePath'));
                        elseif exist(pval, 'file')
                            [xdir,xfile] = fileparts(pval);
                            blastPath = fullfile(xdir,xfile);
                        else
                            error(message('bioinfo:blastlocal:InvalidPath'));
                        end
                    end
            end
        end
    end
    if ~queryFlag
        error(message('bioinfo:blastlocal:MissingQuery'));%blastall hangs w/o query
    end
    if ~progFlag
        cmd = [cmd ' -p blastp']; % assume blastp as default
    end
else % nargin == 1
    cmd = [' ' varargin{:}]; % add a space before BLASTALL command string
    options = strread(cmd, '%s');
    optm = options(find(strcmp(options, '-m'), 1)+1); % option m
    opto = options(find(strcmp(options, '-o'), 1)+1); % option o
    opti = options(find(strcmp(options, '-i'), 1)+1); % option i
    
    if isempty(opti)
       error(message('bioinfo:blastlocal:MissingQuery')); % blastall hangs w/o query
    end
        
    if ~isempty(optm)
        mopt = str2double(optm);
    end
    if ~isempty(opto)
        tofile = char(opto);
    end
    
end

if isempty(tofile) && nargout==0
   return; 
end

%=== run blast 
if isempty(blastPath) % BLASTALL must be in path
    [blastStatus, blastResult]= system(['blastall ' cmd]);
else
    % Quote blastPath as it might contain spaces
    [blastStatus, blastResult] = system(['"' blastPath '"' cmd]);
end

%=== Delete temp file
if tfilename
    delete(tfilename);
end

%=== parse results 
if blastStatus ~= 0
    error(message('bioinfo:blastlocal:blastError', blastResult));
else
    if nargout > 0
        if ~isempty(tofile) % need to first open the file with the report
            fid = fopen(tofile);
            blastResult = fread(fid,'*char')';
            fclose(fid);
        end
        try
            results = blastreadlocal(blastResult, mopt);
        catch le
            if strncmpi(le.identifier,'bioinfo:',8) % did we catch one of our errors?
                rethrow(le)
            elseif isempty(tofile)
                error(message('bioinfo:blastlocal:UnableToReadAskSave'));
            else
                error(message('bioinfo:blastlocal:UnableToRead', tofile))
            end
        end
    end
end


