function info = mzxmlinfo(filename,varargin)
%MZXMLINFO Information about mzXML file.
%   INFO = MZXMLINFO(FILENAME) returns a structure whose fields contain
%   information about a mzXML file.  FILENAME is a string containing a file
%   name, or a path and a file name, of an mzXML file that conforms to the
%   mzXML 2.1 specification or earlier specifications.  INFO is a structure
%   with the following fields:
%
%       Filename - Name of the file
%
%       FileModDate - Modification date of the file
%
%       FileSize - Size of the file in bytes
%
%   Optional mzXML Schema Attributes:
%
%       NumberOfScans - Number of scans in the file
%
%       StartTime - Run start time
%
%       EndTime - Run end time
%
%       DataProcessingIntensityCutoff - Minimum intensity value for an m/z
%
%       DataProcessingCentroided - Indicates data is centroided (T or F)
%
%       DataProcessingDeisotoped - Indicates data is deisotoped (T or F)
%
%       DataProcessingChargeDeconvoluted - Indicates data is deconvoluted
%       (T or F)
%
%       DataProcessingSpotIntegration - For LC-MALDI experiments, peaks
%       eluting over multiple spots that have been integrated to a single
%       spot. (T or F)
%
%       NOTE: If optional attributes are not in the mzXML file the field
%       value will be set to 'N/A'.
%
%   MZXMLINFO(...,'NUMOFLEVELS',T/F) scans the entire mzXML file to
%   determine the number of MS levels.  This will return an additional
%   field, NumberOfMSLevels, in INFO that indicates the number of scan
%   levels in the mzXML file.  Default is false.
%
%   Example:
%
%       % Get mzXML file information. 
%       info = mzxmlinfo('results.mzxml')
%
%       % Along with file information, also get number of MS levels.
%       info = mzxmlinfo('results.mzxml','numoflevels',true)
%
%   Note that the file results.mzxml is not provided. Sample files can be
%   found at http://sashimi.sourceforge.net/repository.html.
%
%   See also JCAMPREAD, MZCDF2PEAKS, MZCDFINFO, MZCDFREAD, MZXML2PEAKS,
%   MZXMLREAD, TGSPCINFO, TGSPCREAD.

% Copyright 2008-2012 The MathWorks, Inc.


msLevelFlag = 0;
bioinfochecknargin(nargin,1,mfilename);

if  nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:mzxmlinfo:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'numoflevels'};
    for j=1:2:nargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:mzxmlinfo:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:mzxmlinfo:AmbiguousParameterName', pname));
        else
            if bioinfoprivate.opttf(pval, okargs{k}, mfilename);
                msLevelFlag = pval;
            end
        end
    end
end

% Check if file exists
if ~exist(filename,'file')
    error(message('bioinfo:mzxmlinfo:invalidFilename', filename));
end

%check if filename contains full path, if not get it
if isempty(regexp(filename,filesep,'once'))
    filename = which(filename);
end

% Grab first few lines of file to check if XML file and mzXML format
fid = fopen(filename,'rt');
str = fread(fid,120,'*char')';
fclose(fid);

% Check if XML file
if isempty(regexp(str,'<\?xml','once'))
    error(message('bioinfo:mzxmlinfo:missingXMLdeclaration', filename, filename));
end

% Check if mzXML format
if isempty(regexp(str,'<mzXML|<msRun','once'))
    error(message('bioinfo:mzxmlinfo:notValidMZXMLFile', filename));
end


% Instantiate filtered parser
filter = com.mathworks.toolbox.bioinfo.util.xml.mzXMLInfoFilter();
biostax = com.mathworks.toolbox.bioinfo.util.xml.BioStAX();
errorCallback = handle(biostax.getNotifyFailedParserCallback());
handle.listener(errorCallback,'delayed',{@notifyFailedParserCallback});

parser = biostax.createStreamFilterParser(filename,filter);

if isempty(parser)
    error(message('bioinfo:mzxmlinfo:FailedParser', filename));    
end


%Create output structure
info = struct('Filename',[],...
    'FileModDate',[],...
    'FileSize',[],...
    'NumberOfScans','N/A',...
    'StartTime','N/A',...
    'EndTime','N/A',...
    'DataProcessingIntensityCutoff','N/A',...
    'DataProcessingCentroided','N/A',...
    'DataProcessingDeisotoped','N/A',...
    'DataProcessingChargeDeconvoluted','N/A',...
    'DataProcessingSpotIntegration','N/A');


% Get filename, filesize and date
fInfo = dir(filename);

info.Filename = fInfo.name;
info.FileModDate = fInfo.date;
info.FileSize = fInfo.bytes;

if msLevelFlag
    info.NumberOfMSLevels = 0;
end

%Parse mzXML file
while parser.hasNext()
    if parser.isStartElement()
        switch char(parser.getLocalName())

            % Get scanCount, startTime and endTime from msRun element
            case 'msRun'
                if parser.isStartElement()
                    attrNum = parser.getAttributeCount();
                    for ind =0:attrNum-1
                        if parser.getAttributeName(ind).toString().equals('scanCount')
                            info.NumberOfScans = double(java.lang.Integer(parser.getAttributeValue(ind)));
                        elseif parser.getAttributeName(ind).toString().equals('startTime')
                            info.StartTime = char(parser.getAttributeValue(ind));
                        elseif parser.getAttributeName(ind).toString().equals('endTime')
                            info.EndTime = char(parser.getAttributeValue(ind));
                        end
                    end
                end

                % Get dataProcessing steps (centroided, deisotoped or
                % chargeDeconvoluted)
            case 'dataProcessing'

                if parser.isStartElement()
                    attrNum = parser.getAttributeCount();
                    for ind =0:attrNum-1
                        if parser.getAttributeName(ind).toString().equals('intensityCutoff')
                            info.DataProcessingIntensityCutoff = double(java.lang.Integer(parser.getAttributeValue(ind)));
                        elseif parser.getAttributeName(ind).toString().equals('centroided')
                            if parser.getAttributeValue(ind).equals('1')
                                info.DataProcessingCentroided = 'true';
                            else
                                info.DataProcessingCentroided = 'false';
                            end

                        elseif parser.getAttributeName(ind).toString().equals('deisotoped')
                            if parser.getAttributeValue(ind).equals('1')
                                info.DataProcessingDeisotoped = 'true';
                            else
                                info.DataProcessingDeisotoped = 'false';
                            end
                        elseif parser.getAttributeName(ind).toString().equals('chargeDeconvoluted')
                            if parser.getAttributeValue(ind).equals('1')
                                info.DataProcessingChargeDeconvoluted = 'true';
                            else
                                info.DataProcessingChargeDeconvoluted = 'false';
                            end
                        elseif parser.getAttributeName(ind).toString().equals('spotIntegration')
                            if parser.getAttributeValue(ind).equals('1')
                                info.DataProcessingChargeSpotIntegration= 'true';
                            else
                                info.DataProcessingChargeSpotIntegration = 'false';
                            end
                        end
                    end
                end

                %Close parser and return if NUMOFLEVELS = false
                if ~msLevelFlag
                    parser.close();
                    return;
                end
                
            %Get scanLevel from scan element
            case 'scan'
                if parser.isStartElement()
                    if double(java.lang.Integer(parser.getAttributeValue(1))) > info.NumberOfMSLevels
                        info.NumberOfMSLevels = info.NumberOfMSLevels+1;
                    end
                end

        end

    else
    end
    parser.next();
end

%close parser
parser.close();

end % MZXMLINFO

%-----------------------------------------------------------
% Error message callbacks
%-----------------------------------------------------------
function notifyFailedParserCallback(hsrc,hevt)%#ok

     fileFailed.loaded = hevt.JavaEvent.loaded;
     fileFailed.errormsg = hevt.JavaEvent.errormsg;
     fileFailed.filename = hevt.JavaEvent.filename;
   
end
