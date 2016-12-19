function data=emblread(embltext,varargin)
%EMBLREAD read in EMBL-EBI format data files.
%
%   DATA = EMBLREAD(FILENAME) imports the EMBL-EBI format data from file
%   FILENAME and creates a structure DATA with fields corresponding to the
%   EMBL-EBI two-character line type code. Each line type code is stored as
%   a separate element of the structure.
%
%   FILENAME can also be a URL or a MATLAB character array that contains
%   the text of an EMBL format file. 
%
%   DATA = EMBLREAD('FILENAME.EXT','SEQUENCEONLY',true) reads only the
%   sequence information.
% 
%   Based on release 107 of EMBL-Bank flat file format.
% 
%   Remarks:
%   [1] Topology information was not included in EMBL flat files before
%   release 87 of the database. When reading a file created before release
%   87, EMBLREAD returns an empty Identification.Topology field.
%   [2] The entry name is no longer displayed in the ID line of the EMBL
%   flat files in release 87. When reading a file created in release 87,
%   EMBLREAD returns the accession number in the Identification.EntryName
%   field.
%
%   Examples:
%
%       % Get EMBL data and save it to a file.
%       emblout = getembl('X00558','TOFILE','X00558.ebi')
%
%       % In subsequent MATLAB sessions you can use emblread to access the
%       % local copy from disk instead of accessing it from the EMBL site.
%       data = emblread('X00558.ebi')
%
%   See also FASTAREAD, GENBANKREAD, GETEMBL, SEQVIEWER.

% Copyright 2002-2012 The MathWorks, Inc.



[SeqOnly, getPreambleText] = parse_inputs(varargin{:});

if ~ischar(embltext) && ~iscellstr(embltext)
    error(message('bioinfo:emblread:InvalidInput'));
end

% guess what type of input we have -- string, file or URL
if iscellstr(embltext)
    % do not mess with it, just put it to char and try to decipher in the try-catch trap
    % embltext=embltext;
elseif size(embltext,1)==1 %it is char, lets check if it has an url or a file before try to decipher it
    if  ~isempty(strfind(embltext(1:min(10,end)), '://'))
        % must be a URL
        if (~usejava('jvm'))
            error(message('bioinfo:emblread:NoJava'))
        end
        try
            embltext = urlread(embltext);
        catch allExceptions
            error(message('bioinfo:emblread:CannotReadURL', embltext));
        end
        embltext = strrep(embltext,'&amp;','&');
        % make each line a separate row in a string array
        embltext = textscan(embltext,'%s','delimiter','\n','whitespace','');
        embltext = embltext{1};
    elseif (exist(embltext,'file') || exist(fullfile(pwd,embltext),'file'))
        fid = fopen(embltext,'r');
        embltext = textscan(fid,'%s','delimiter','\n','whitespace','');
        embltext = embltext{1};
        fclose(fid);
    else
        embltext = textscan(embltext,'%s','delimiter','\n','whitespace','');
        embltext = embltext{1};
    end
end

if ischar(embltext)
    embltext = cellstr(embltext);
end
% If the input is a string of EMBL data then words ID must be present
if size(embltext,1)==1 || isempty(strfind(embltext{1},'ID'))
    error(message('bioinfo:emblread:NonMinimumRequiredFields'))
end

if ~iscellstr(embltext)
    error(message('bioinfo:emblread:InvalidFormat'));
end

refcount=0;
% remove XX lines
embltext(~cellfun('isempty',regexp(embltext,'^XX','once'))) = [];

% Create the output structure
numrecords = sum(~cellfun('isempty',regexp(embltext,'^ID')));
data(numrecords).Identification.EntryName = '';

%line number
ln=1;
recordcount=1;

while 1

    lnPreambleTextStarts = ln;

    %ID-Identification
    idLine = textscan(strtrim(embltext{ln}(6:end)),'%s','delimiter',';');
    % check if ID line is formatted as EMBL Release 87
    if numel(idLine{1})==7
        formatIsR87 = true;
        data(recordcount).Identification.EntryName = idLine{1}{1};
        data(recordcount).Identification.Version = strrep(idLine{1}{2},'SV ','');
        data(recordcount).Identification.Topology = idLine{1}{3};
        data(recordcount).Identification.Molecule = idLine{1}{4};
        data(recordcount).Identification.DataClass = idLine{1}{5};
        data(recordcount).Identification.Division = idLine{1}{6};
        data(recordcount).Identification.SequenceLength = strrep(idLine{1}{7},'.','');
    elseif numel(idLine{1})==4
        % then we try read the flat format before Release 87
        formatIsR87 = false;
        [tmp1,tmp2] = strread(idLine{1}{1},'%s%s');
        data(recordcount).Identification.EntryName = tmp1{1};
        data(recordcount).Identification.Version = ''; % fill it later
        data(recordcount).Identification.Topology = ''; % info not available before Release 87
        data(recordcount).Identification.Molecule = idLine{1}{2};
        data(recordcount).Identification.DataClass = tmp2{1};
        data(recordcount).Identification.Division = idLine{1}{3};
        data(recordcount).Identification.SequenceLength = strrep(idLine{1}{4},'.','');
    else
        error(message('bioinfo:emblread:CannotInterpretIDline'));
    end
    ln=ln+1;

    %AC-Accession number
    data(recordcount).Accession=deblank(embltext{ln}(6:end));
    data(recordcount).Accession = strrep(data(recordcount).Accession,';','');
    ln=ln+1;

    %SV-Sequence Version
    if formatIsR87 == true;
        data(recordcount).SequenceVersion = ...
            [ data(recordcount).Identification.EntryName '.' ...
            data(recordcount).Identification.Version ];
    else
        tmp = deblank(embltext{ln}(6:end));
        data(recordcount).SequenceVersion = tmp;
        ln=ln+1;
        data(recordcount).Identification.Version = ...
            cell2mat(regexp(tmp,'\w*.(\w*)','tokens','once'));
    end
    
    %PR-Project (INSDC project identifier)
    if strcmp(embltext{ln}(1:2),'PR')
        data(recordcount).Project=deblank(embltext{ln}(6:end));
        data(recordcount).Project = strrep(data(recordcount).Project,';','');
        ln=ln+1;
    end
    
    %DT-Date
    data(recordcount).DateCreated=deblank(embltext{ln}(6:end));
    ln=ln+1;
    data(recordcount).DateUpdated=deblank(embltext{ln}(6:end));
    ln=ln+1;

    %DE-Description
    [data(recordcount).Description, ln] = extractfield(embltext,ln,'DE');

    %KW-Keyword
    [data(recordcount).Keyword, ln] = extractfield(embltext,ln,'KW');

    %OS-Organism Species
    [data(recordcount).OrganismSpecies, ln] = extractfield(embltext,ln,'OS');

    %OC-Organism Classification
    [data(recordcount).OrganismClassification, ln] = extractfield(embltext,ln,'OC');

    %OG-Organelle
    [data(recordcount).Organelle, ln] = extractfield(embltext,ln,'OG');

    %Reference
    while matchstart(embltext{ln},'RN') %ref name
        refcount=refcount+1;

        %RN-Reference name
        data(recordcount).Reference{refcount}.Number=deblank(embltext{ln}(6:end));
        ln=ln+1;

        %RC-Reference Comment
        [data(recordcount).Reference{refcount}.Comment, ln] = extractfield(embltext,ln,'RC');

        %RP-Reference Position
        [data(recordcount).Reference{refcount}.Position, ln] = extractfield(embltext,ln,'RP');

        %RX-Reference Database identifier
        data(recordcount).Reference{refcount}.MedLine = '';
        data(recordcount).Reference{refcount}.PubMed = '';
        [databaseText, ln] = extractfield(embltext,ln,'RX');
        databaseText = lower(databaseText);
        medlineRow = strmatch('medline',databaseText);
        if medlineRow
            [junk, medlineText] = strtok(databaseText(medlineRow,:)); 
            medlineText(medlineText == '.') = '';
            medlineText(medlineText == ' ') = '';
            data(recordcount).Reference{refcount}.MedLine = medlineText;
        end

        pubmedRow = strmatch('pubmed',databaseText);
        if pubmedRow
            [junk, pubmedText] = strtok(databaseText(pubmedRow,:)); 
            pubmedText(pubmedText == '.') = '';
            pubmedText(pubmedText == ' ') = '';
            data(recordcount).Reference{refcount}.PubMed = pubmedText;
        end

        %RG-Reference Group
        [data(recordcount).Reference{refcount}.Group, ln] = extractfield(embltext,ln,'RG');

        %RA-Reference Authors
        [data(recordcount).Reference{refcount}.Authors, ln] = extractfield(embltext,ln,'RA');

        %RT-Reference Title
        [data(recordcount).Reference{refcount}.Title, ln] = extractfield(embltext,ln,'RT');

        %RL-Reference Location
        [data(recordcount).Reference{refcount}.Location, ln] = extractfield(embltext,ln,'RL');

    end

    %DR-Database Cross Reference
    [data(recordcount).DatabaseCrossReference, ln] = extractfield(embltext,ln,'DR');

    lnPreambleTextEnds = ln;

    %CC-free text Comment
    [data(recordcount).Comments, ln] = extractfield(embltext,ln,'CC');

    %AH Assembly header
    [data(recordcount).Assembly, ln] = extractfield(embltext,ln,'A');

    %FH Features header
    [data(recordcount).Feature, ln] = extractfield(embltext,ln,'F');

    % Need to strip empty lines from Feature field for compatibility with
    % other tools
    emptyFeatureRow = cellfun('isempty',cellstr(data(recordcount).Feature));
    data(recordcount).Feature(emptyFeatureRow,:) = '';
    
    oth = false;
    %SQ sequence header
    while matchstart(embltext{ln},'SQ')
        temp=embltext{ln}(14:end);
        [numbp,temp]=strtok(temp,'BP'); temp= temp(4:end);%#ok
        [numa,temp]=strtok(temp,'A'); temp = temp(3:end);%#ok
        [numc,temp]=strtok(temp,'C'); temp = temp(3:end);%#ok
        [numg,temp]=strtok(temp,'G'); temp = temp(3:end);%#ok
        [numt,temp]=strtok(temp,'T'); %#ok
        if strfind(temp,'other'),
            temp = temp(3:end);
            numo = strtok(temp,'o');
            oth = true;
        end

        data(recordcount).BaseCount.BP= str2double(numbp);
        data(recordcount).BaseCount.A = str2double(numa);
        data(recordcount).BaseCount.C = str2double(numc);
        data(recordcount).BaseCount.G = str2double(numg);
        data(recordcount).BaseCount.T = str2double(numt);
        data(recordcount).BaseCount.Other = [];
        if oth
            data(recordcount).BaseCount.Other = str2double(numo);
        end

        ln=ln+1;
    end

    %Read in sequence until end of entry, //
    data(recordcount).Sequence=[];
    startLn = ln;
    while ~matchstart(embltext{ln},'//') && ~matchstart(embltext{ln},'ID')
        ln=ln+1;
    end
    if ln > startLn
        data(recordcount).Sequence = char(embltext{startLn:ln});
        data(recordcount).Sequence = data(recordcount).Sequence(:,6:72);
        data(recordcount).Sequence = data(recordcount).Sequence';
        data(recordcount).Sequence = data(recordcount).Sequence(:);
        data(recordcount).Sequence = data(recordcount).Sequence';
        data(recordcount).Sequence = strrep(data(recordcount).Sequence,' ','');
    end
    if getPreambleText
        data(recordcount).PreambleText = embltext{lnPreambleTextStarts:lnPreambleTextEnds-1};
    end

    while  (ln < size(embltext,1)) && ~matchstart(embltext{ln},'ID')
        ln=ln+1;
    end

    if matchstart(embltext{ln},'ID') %see if line is part of next record
        recordcount = recordcount+1;
        continue
    else
        if SeqOnly == true
            if recordcount==1
                data=data(recordcount).Sequence;
            else
                data={data.Sequence};
            end
        end
        return
    end
end

%-------------------------------------------------------------------------%
function tf = matchstart(string,pattern)
%MATCHES start of string with pattern

tf = ~isempty(regexp(string,['^',pattern],'once'));
%-------------------------------------------------------------------------%
function [data, outLine] = extractfield(embltext,ln,lineID)
%Extracts a field from the embltext cellstr

startLn = ln;
while matchstart(embltext{ln},lineID)
    ln=ln+1;
end
data = char(embltext{startLn:ln-1});
data = data(:,6:end);
outLine = ln;

%-------------------------------------------------------------------------%
function [SeqOnly, getPreambleText] = parse_inputs(varargin)
% Parse the varargin parameter/value inputs

% Check that we have the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:emblread:IncorrectNumberOfArguments', mfilename));
end

% set defaults
SeqOnly = false;
getPreambleText = false;
okargs = {'sequenceonly','preambletext'};

for j=1:2:nargin
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1 % sequence
            SeqOnly = bioinfoprivate.opttf(pval,okargs{k},mfilename);
        case 2  % 'preambletext'
            getPreambleText = bioinfoprivate.opttf(pval,okargs{k},mfilename);
    end
end

