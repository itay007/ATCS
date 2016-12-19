function emblout=getembl(accessnum,varargin)
%GETEMBL retrieves sequence information from the EMBL-EBI database.
%
%   GETEMBL(ACCESSNUM) searches for the accession number in the EBI
%   database (http://www.ebi.ac.uk/embl), and returns a structure
%   containing information about the sequence.
%
%   GETEMBL(...,'TOFILE',FILENAME) returns the data as a structure and
%   saves data in the file FILENAME in the EMBL-EBI data format.
%
%   GETEMBL(...,'SEQUENCEONLY',true) returns only the sequence
%   information and none of the metadata.
%
%   For more details about the EMBL database, see
%       http://www.ebi.ac.uk/embl/Documentation/index.html
%
%   Examples:
%
%       % Retrieve the data for rat liver apolipoprotein A-I.
%       emblout = getembl('X00558')
%
%       % Now retrieve just the sequence information for the same protein.
%       seq = getembl('X00558','SequenceOnly',true)
%
%   See also EMBLREAD, GETGENBANK, GETGENPEPT, GETPDB, SEQVIEWER.

% Copyright 2002-2012 The MathWorks, Inc.



sequenceOnly=false;
savefile=false;
getPreambleText = false;
if nargin > 1
    if rem(nargin,2)== 0
        error(message('bioinfo:getembl:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'tofile','sequenceonly','preambletext'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname))); 
        if isempty(k)
            error(message('bioinfo:getembl:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:getembl:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % save file data
                    filename=pval;
                    savefile=true;
                case 2  % sequenceOnly
                    sequenceOnly = bioinfoprivate.opttf(pval);
                    if isempty(sequenceOnly)
                        error(message('bioinfo:getembl:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 3  % 'preambletext'
                    getPreambleText = bioinfoprivate.opttf(pval);
            end
        end
    end
end

% convert accessnum to a string if it is a number
if isnumeric(accessnum)
    accessnum = num2str(accessnum);
end

% error if accessnum isn't a string
if ~ischar(accessnum)
    error(message('bioinfo:getembl:NotString'))
end

accessnum = strtrim(accessnum);
if any(isspace(accessnum))
    error(message('bioinfo:getembl:InvalidAccession'))
end

%URL that is used to call database
%see http://www.ebi.ac.uk/cgi-bin/emblfetch for more information
retrieveURL = ['http://www.ebi.ac.uk/cgi-bin/dbfetch?db=EMBL&id=' accessnum '&style=raw']; 
temp = urlread(retrieveURL);

%search for returned text indicating that the accession number was not
%found in any files
if any(strfind(temp,'ERROR')) || any(strfind(temp, 'No entries found'))
    error(message('bioinfo:getembl:AccessionNotFound', accessnum));
end

%make each line a separate row in string array
embldata = char(strread(temp,'%s','delimiter','\n','whitespace',''));

%pass to EMBLREAD to create structure
if getPreambleText
    emblout=emblread(embldata,'sequenceOnly',sequenceOnly,'preambletext',true);
else
    emblout=emblread(embldata,'sequenceOnly',sequenceOnly);
end

%  write out file
if savefile == true
    fid = fopen(filename,'wt') ;
    if fid == (-1)
        error(message('bioinfo:getembl:CouldNotOpenFile', filename));
    else
        rows = size(embldata,1);

        for rcount=1:rows-1,
            fprintf(fid,'%s\n',embldata(rcount,:));
        end

        fprintf(fid,'%s',embldata(rows,:));
        fclose(fid);
    end
end

if sequenceOnly == false % emblout is an structure
    for i = 1:numel(emblout)    
        emblout(i).RetrieveURL = retrieveURL;
    end
    % in case there is no output assigned we show the hyperlinks in the
    % command window
    if ~nargout && usejava('desktop')
        for i = 1:numel(emblout)    
            disp(rmfield(emblout(i),'RetrieveURL'));
            retrieveURL = ['http://www.ebi.ac.uk/cgi-bin/dbfetch?db=EMBL&id=' accessnum ];
            disp([char(9) 'RetrieveURL: <a href="' retrieveURL '">' emblout(i).Accession '</a>']);
        end
        clear emblout
    end
end

    
   
