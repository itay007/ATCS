function theURL = probesetlink(affyStruct,theID,varargin)
% PROBESETLINK links to NetAffx Web site
%
%   PROBESETLINK(AFFYSTRUCT,PS) displays information from the NetAffx Web
%   site about probe set ID from the CHP or CDF structure AFFYSTRUCT.
%   PS can be the index of the probe set or the probe set name.
%
%   URL = PROBESETLINK(AFFYSTRUCT,PS) returns the URL for the information.
%
%   PROBESETLINK(...,'SOURCE',true) links to the data source (e.g., GenBank,
%   Flybase) for the probe set, PS.
%
%   PROBESETLINK(...,'BROWSER',true) displays the information in the system
%   Web browser.
%
%   URL = PROBESETLINK(...,'NODISPLAY',true) returns the URL but does not
%   open a browser.
%
%   Note that the NetAffx web site requires you to register and provide a
%   username and password.
%
%   Example:
%       chpStruct = affyread('Ecoli-antisense-121502.chp',...
%                            'C:\Affymetrix\LibFiles\Ecoli_ASv2');
%       probesetlink(chpStruct,'argG_b3172_at');
%
%   See also AFFYDEMO, AFFYREAD, CELINTENSITYREAD, PROBELIBRARYINFO,
%   PROBESETLOOKUP, PROBESETPLOT, PROBESETVALUES.

%   Affymetrix and NetAffx are registered trademarks of Affymetrix, Inc.

% Copyright 2003-2010 The MathWorks, Inc.


bioinfochecknargin(nargin,2,mfilename);

browserFlag = true;
sourceFlag = false;
displayFlag = true;

% Now check that the AFFYSTRUCT is a struct and CDF or CEL structs
if  ~isstruct(affyStruct)
    error(message('bioinfo:probesetlink:CdfStructNotStruct'));
end

if ~isfield(affyStruct,'Name') || ~isfield(affyStruct,'ChipType') || ~isfield(affyStruct,'ProbeSets')
    error(message('bioinfo:probesetlink:UnsupportedAffyStruct', mfilename));
end

%Check that PS input is not a cell array of strings or vector
if iscellstr(theID)
    error(message('bioinfo:probesetlink:BadIDInputString'));
end

if isnumeric(theID)&&numel(theID)>1
    error(message('bioinfo:probesetlink:BadIDInputNumeric'));
end

% get the ID
if ischar(theID)
    ID = find(strncmp(theID,{affyStruct.ProbeSets.Name},numel(theID)));
    if isempty(ID)
        error(message('bioinfo:probesetlink:UnknownProbeName', theID));
    elseif length(ID)>1
        warning(message('bioinfo:probesetlink:AmbiguousProbeName', theID));
        ID = ID(1);
    end
else
    ID = theID;
end

% deal with the various inputs
if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:probesetlink:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'browser','source','nodisplay'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:probesetlink:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:probesetlink:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % browser
                    browserFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 2  % source
                    sourceFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                  case 3  % no display
                    displayFlag = ~bioinfoprivate.opttf(pval,okargs{k},mfilename);
            end
        end
    end
end

% Extract the data from the big struct
theChip = affyStruct.ChipType;
theName = affyStruct.ProbeSets(ID).Name;

%Grab Source URL using PROBESETLOOKUP
if sourceFlag
    try
        probeSetStruct = probesetlookup(affyStruct,theName); 
        link = probeSetStruct.SourceURL;
    catch theException
        throw(theException);
    end

else
    theName = strrep(theName,'#','%23');
    link = sprintf('https://www.affymetrix.com/LinkServlet?array=%s&probeset=%s',...
        theChip,theName);
end
% open the link in a browser
if displayFlag
    if browserFlag
        web(link,'-browser');
    else
        web(link);
    end
end
% set the output if required
if nargout > 0
    theURL = link;
end
