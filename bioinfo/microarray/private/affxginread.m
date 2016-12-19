function ginStruct = affxginread(filename,pathdir,libdir) %#ok
% AFFXGINREAD reads in GIN files that contain gene identifiers
%
%   GinStruct = AFFXGINREAD(FILE) reads an Affymetrix gene identifier file
%   FILE and creates a structure GinStruct.
%
%   See also AFFYREAD, AGFEREAD, CELINTENSITYREAD, GPRREAD,
%   PROBELIBRARYINFO, PROBESETLINK, PROBESETLOOKUP, PROBESETPLOT,
%   PROBESETVALUES, SPTREAD.
%
%   GeneChip and Affymetrix are registered trademarks of Affymetrix, Inc.

% Copyright 2004-2007 The MathWorks, Inc.



fullFileName = '';
if exist(filename,'file')
    % we have the file
    fullFileName = filename;
else
    if nargin > 1
        fullFileName = fullfile(libdir,filename);
    end
end

fid = fopen(fullFileName,'r');
if fid < 0
    error(message('bioinfo:affyread:BadFile', fullFileName));
end
% Read in the header lines
line = fgetl(fid);
try
    if strncmpi('version',line,7) == 1  % versioned files (version 2 and higher?)

        version = strread(line,'%*s%d','delimiter','=');
        line = fgetl(fid);
        name = strread(line,'%*s%s','delimiter','=');
        line = fgetl(fid);
        numURLs = strread(line,'%*s%d','delimiter','=');

        % allocate space for the references and URLS.
        urls = cell(numURLs,1);
        refs = cell(numURLs,1);
        for i = 1: numURLs
            line = fgetl(fid);
            [refs(i) , urls(i)] = strread(line,'%s%*s%s',1,'delimiter',';');
        end

        % Read in the rest of the file skipping the line with column headings.
        out = textscan(fid,'%d%s%s%s%s%*[^\n]','headerlines',1,'delimiter','\t');
        % close the file
        fclose(fid);

        % Set output structure
        ginStruct.Name = name{:};
        ginStruct.Version = version;
        ginStruct.ProbeSetName = out{:,4};
        ginStruct.ID = out{:,2};
        ginStruct.Description = out{:,5};
        if numURLs == 1
            ginStruct.SourceNames = refs{:};
            ginStruct.SourceURL = urls{:};
            ginStruct.SourceID = 1;
        elseif numURLs == 0
            ginStruct.SourceNames = '';
            ginStruct.SourceURL = '';
            ginStruct.SourceID = 0;
        else
            ginStruct.SourceNames = refs;
            ginStruct.SourceURL = urls;
            [dummy, ginStruct.SourceID ] = ismember(out{:,3}, refs);%#ok
        end
    else  % old files e.g. Hu6800.GIN
        name = strread(line,'%*s%s','delimiter','=');
        line = fgetl(fid);
        refname = strread(line,'%*s%s','delimiter','=');
        line = fgetl(fid);
        [dummy,urls] = strtok(line,'='); %#ok
        urls = urls(2:end);
        ginStruct.Name = name{:};
        ginStruct.Version = 1;
        out = textscan(fid,'%d%s%s%s%s%*[^\n]','headerlines',1,'delimiter','\t');
        fclose(fid);
        ginStruct.ProbeSetName = out{:,3};
        ginStruct.ID = out{:,2};
        ginStruct.Description = out{:,4};
        ginStruct.SourceNames =refname{:};
        ginStruct.SourceURL = urls;
        ginStruct.SourceID = 1;
    end
catch
    error(message('bioinfo:affxginread:badGINfile', filename));
end
