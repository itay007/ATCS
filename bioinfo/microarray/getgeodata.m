function [geoStruct,str] =getgeodata(accessnum,varargin)
% GETGEODATA retrieves Gene Expression Omnibus (GEO) data. 
%
%   GEOSTRUCT = GETGEODATA(ACCESSNUM) searches for the accession number in
%   the Gene Expression Omnibus (GEO) database, and returns a structure
%   containing information for the object. 
%
%   GEOSTRUCT = GETGEODATA(...,'TOFILE',FILENAME) saves the data returned
%   from the database in the file FILENAME.
%
%   Note that currently Sample (GSM), DataSet (GDS), Series (GSE) and
%   Platform (GPL) records are supported. 
%
%   Example:
%   
%          % Get a sample file
%          geoSample = getgeodata('GSM1768')
%
%          % Get a data set and save it to a file
%          geoDataSet = getgeodata('GDS2602','tofile','gds2602.txt')
%
%          % Get a series data matrix
%          geoSeries = getgeodata('GSE11287')
%
%          % Get a Platform record
%          geoSeries = getgeodata('GPL74')
%
%   See http://www.ncbi.nlm.nih.gov/About/disclaimer.html for information
%   about using the GEO database.
%
%   See also GEOSERIESREAD, GEOSOFTREAD, GETGENBANK, GETGENPEPT.


% Copyright 2003-2012 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);

if ~usejava('jvm')
    error(message('bioinfo:getgeodata:NeedJVM', mfilename));
end
tofile = false;


if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:getgeodata:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'tofile',''};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:getgeodata:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:getgeodata:AmbiguousParameterName', pname));
        else
            switch(k)
                case  1  % tofile
                    if ischar(pval)
                        tofile = true;
                        filename = pval;
                    end
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
    error(message('bioinfo:getgeodata:NotString'))
end

% create the url that is used
% see
%    http://www.ncbi.nlm.nih.gov/entrez/query/static/linking.html
% for more information
searchurl = sprintf('http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?form=text&acc=%s&view=full',accessnum);

% get the html file that is returned as a string
try
    str=urlread(searchurl);
catch theException
    error(message('bioinfo:getgeodata:URLProblem'));
end

% search for text indicating that there weren't any files found
notfound=strfind(str,'No items found');

% string was found, meaning no results were found
if ~isempty(notfound),
    error(message('bioinfo:getgeodata:BadAccessionNumber', accessnum)) ;
end

% GDS datasets
gdsBrowser = strfind(str,'gds_browse.cgi');
if ~isempty(gdsBrowser)
    % We expect the file on the FTP site to be called ACCNUM.soft.gz
    accessnum = upper(accessnum);
    softFileGZName = sprintf('%s.soft.gz',accessnum);
    softFileName = sprintf('%s.soft',accessnum);
    % Make a temporary directory for copying the file
    tempDirName = tempname;
    mkdir(tempDirName);
    % open FTP connection to NCBI and CD to GDS directory and get the file
    ftpConnection = ftp('ftp.ncbi.nih.gov');
    try
        cd(ftpConnection,'pub/geo/DATA/SOFT/GDS'); %#ok<MCCD>
        mget(ftpConnection,softFileGZName,tempDirName);
        close(ftpConnection);
    catch theErr
        % try to close the connection if there is an error
        close(ftpConnection);
        rethrow(theErr);
    end
    % gunzip the file and remove the .gz file
    localCopy = fullfile(tempDirName,softFileGZName);
    gunzip(localCopy);
    delete(localCopy);
    % Read the file
    geoStruct = geosoftread(fullfile(tempDirName,softFileName));
    if tofile
        copyfile(fullfile(tempDirName,softFileName),filename);
    end
    % Clean up
    delete(fullfile(tempDirName,softFileName));
    rmdir(tempDirName)
    return
end

%GSE series
gseFlag = strncmpi(str,'^SERIES',7);
if gseFlag ~= 0
     % We expect the file on the FTP site ftp://ftp.ncbi.nih.gov
     % in pub/geo/DATA/SeriesMatrix/%ACCNUM%
     % named ACCNUM_series_matrix.txt.gz
    accessnum = upper(accessnum);
    seriesFileGZName = sprintf('%s*_series_matrix.txt.gz',accessnum);
    seriesFileName = sprintf('%s_series_matrix.txt',accessnum);
    % Make a temporary directory for copying the file
    tempDirName = tempname;
    mkdir(tempDirName);
    % open FTP connection to NCBI and CD to GDS directory and get the file
    ftpConnection = ftp('ftp.ncbi.nih.gov');
    try
        cd(ftpConnection,sprintf('pub/geo/DATA/SeriesMatrix/%s',accessnum)); %#ok<MCCD>
        mget(ftpConnection,seriesFileGZName,tempDirName);
        close(ftpConnection);
    catch theErr
        % try to close the connection if there is an error
        close(ftpConnection);
        eid = theErr.identifier;
        if strcmpi(eid,'MATLAB:ftp:NoSuchDirectory')
            error(message('bioinfo:getgeodata:NoSuchDirectory'));
        else
            rethrow(theErr)
        end
    end
    % gunzip the file and remove the .gz file
    localCopy = fullfile(tempDirName,seriesFileGZName);
    gunzip(localCopy);
    delete(localCopy);
    if ~exist(fullfile(tempDirName,seriesFileName),'file')
        dirInfo = dir([tempDirName filesep '*.txt']);
        if numel(dirInfo) > 1
            error(message('bioinfo:getgeodata:MultipleGSE', tempDirName));
        else
        ftpPath = sprintf('ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/%s',accessnum);
        error(message('bioinfo:getgeodata:CannotExtractGSE', ftpPath));
        end
    end
    % Read the file
    geoStruct = geoseriesread(fullfile(tempDirName,seriesFileName));
    if tofile
        copyfile(fullfile(tempDirName,seriesFileName),filename);
    end
    % Clean up
    delete(fullfile(tempDirName,seriesFileName));
    rmdir(tempDirName)
    return
end

str = strrep(str,char(0),' ');
geoStruct = geosoftread(str);

%  write out file?
if tofile == true
    writefile = 'Yes';
    % check to see if file already exists
    if exist(filename,'file')
        % use dialog box to display options
        writefile=questdlg(sprintf('The file %s already exists. Do you want to overwrite it?',filename), ...
            '', ...
            'Yes','No','Yes');
    end

    switch writefile,
        case 'Yes',
            if exist(filename,'file')
                fprintf('File %s overwritten.',filename);
            end
            savedata(filename,str);
        case 'No',
            fprintf('File %s not written.',filename);
    end

end

function savedata(filename,str)

fid=fopen(filename,'w');

rows = size(str,1);

for rcount=1:rows-1
    fprintf(fid,'%s\n',str(rcount,:));
end
fprintf(fid,'%s',str(rows,:));

fclose(fid);

