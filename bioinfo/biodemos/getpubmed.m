function pmstruct = getpubmed(searchterm,varargin)
% GETPUBMED Search PubMed database and write results to MATLAB structure
% 
%
% PMSTRUCT = GETPUBMED(SEARCHTERM) searches NCBI PubMed db with SEARCHTERM,
% retrieves citation and abstract info, and writes the search results to 
% PMSTRUCT, a MATLAB structure.
%
% SEARCHTERM is a string of text that is valid for any PubMed search.
%
% For help creating PubMed searches, including valid symbols and tags, see 
% http://www.ncbi.nlm.nih.gov/books/bv.fcgi?rid=helppubmed.chapter.pubmedhelp 
% For example, a PubMed search for all publications related to invasive
% brain cancer could be 'Brain Cancer AND Invasion'.
%
% NOTE: SEARCHTERM cannot contain spaces. Use '+' to replace spaces, 
% e.g. Brain+Cancer+AND+Invasion. 
%
% PMSTRUCT is a MATLAB structure or structure array containing information
% about the article(s) found with SEARCHTERM.
% PMSTRUCT contains the following fields:
%   PubMedID
%   PublicationDate
%   Title
%   Abstract
%   Authors
%   Citation 
%
% GETPUBMED(...,'NUMBEROFRECORDS',MAXNUM) limits the search by
% returning up to the specified MAXNUM of records. Default is 50 records.
%
% GETPUBMED(...,'DATEOFPUBLICATION',PUBDATE) limits the search by date(s) 
% of publication. PUBDATE can be a year, year/month, year/month/day or a 
% range in the following formats: 
%
%   Format:                            Example:
%   YYYY[DP]                           2008[DP]
%   YYYY/MM[DP]                        2008/01[DP]
%   YYYY/MM/DD[DP]                     2008/01/01[DP]
%   YYYY/MM/DD[DP]:YYYY/MM/DD[DP]      2008/1/1[DP]:2008/5/30[DP]
%
% Default is an empty string, which sets no limitations on publication
% date.
%
% Example:   
%
%   % Find all articles related to invasive brain cancer
%   out = getpubmed('brain+cancer+AND+invasion');
%
%   % Find articles related to clinical trials for breast cancer with BRCA1
%   % mutation.  Limit the number articles returned to 10.
%   out = getpubmed('breast+cancer+AND+BRCA1+AND+Clinical+Trial[PT]',...
%   'NumberOfRecords',10)
%
%   % NOTE: [PT] is the tag for Publication Type.  
%   % For more information on PubMed search tags, see
%   % http://www.ncbi.nlm.nih.gov/books/bv.fcgi?rid=helppubmed.section.pubmedhelp.Search_Field_Descrip#pubmedhelp.Publication_Type_PT
%
%   % Find articles that are related to bird flu and published between 
%   % January 1, 2008 and March 1, 2008. 
%   out = getpubmed('bird+flu',...
%   'DateOfPublication','2008/1/1[DP]:2008/3/1[DP]');
%

%   Copyright 2008-2012 The MathWorks, Inc.


% Error checking for required input SEARCHTERM
if(nargin<1)
    error(message('bioinfo:getpubmed:NotEnoughInputArguments'));
end

% Set default settings for property name/value pairs, 
% 'NUMBEROFRECORDS' and 'DATEOFPUBLICATION'
maxnum = 50; % NUMBEROFRECORDS default is 50
pubdate = ''; % DATEOFPUBLICATION default is an empty string

% Parsing the property name/value pairs 
num_argin = numel(varargin);
for n = 1:2:num_argin
    arg = varargin{n};
    switch lower(arg)
        
        % If NUMBEROFRECORDS is passed, set MAXNUM
        case 'numberofrecords'
            maxnum = varargin{n+1};
        
        % If DATEOFPUBLICATION is passed, set PUBDATE
        case 'dateofpublication'
            pubdate = varargin{n+1};          
            
    end     
end

% You access the PubMed database through a search URL, which submits a
% search term and options, and then returns the search results in a
% specified format. The search URL is comprised of a base URL and defined
% parameters.

% Create base URL for PubMed db site 
baseSearchURL = 'http://www.ncbi.nlm.nih.gov/sites/entrez?cmd=search';

% Create five defined parameters that GETPUBMED will use: db (database),
% term (search term), report (report type, such as MEDLINE), format (format
% type, such as text), and dispmax (maximum number of records to display).

% Set db parameter to pubmed
dbOpt = '&db=pubmed';

% Set term parameter to SEARCHTERM and PUBDATE
% (Default PUBDATE is '')
termOpt = ['&term=',searchterm,'+AND+',pubdate];

% Set report parameter to medline
reportOpt = '&report=medline';

% Set format parameter to text
formatOpt = '&format=text';

% Set dispmax to MAXNUM
% (Default MAXNUM is 50)
maxOpt = ['&dispmax=',num2str(maxnum)];

% Create search URL
searchURL = [baseSearchURL,dbOpt,termOpt,reportOpt,formatOpt,maxOpt];

% An example search URL is: 
% http://www.ncbi.nlm.nih.gov/sites/entrez?cmd=search&db=pubmed&term=bird+f
% lu+AND+2008/1/1[DP]:2008/3/1[DP]&report=medline&format=text

% For more information on building URLs to access Entrez databases, see 
% http://www.ncbi.nlm.nih.gov/books/bv.fcgi?rid=helplinks.chapter.linkshelp


% Use the MATLAB function, URLREAD, to submit the search URL, retrieve the
% search results, and return the results, as text in the MEDLINE report
% type, in a character array.  
medlineText = urlread(searchURL);

% Use the MATLAB function, REGEXP, and regular expressions
% to parse and extract the info in the character array, medlineText,
% into hits, a cell array where each cell contains the MEDLINE-formatted
% text for one article.
% For instructions on using regular expressions in MATLAB, see 
% http://www.mathworks.com/access/helpdesk/help/techdoc/index.html?/access/
% helpdesk/help/techdoc/matlab_prog/f0-42649.html 
hits = regexp(medlineText,'PMID-.*?(?=PMID|</pre>$)','match');

% Instantiate the PMSTRUCT structure that will be returned by GETPUBMED
% to contain six fields.
pmstruct = struct('PubMedID','','PublicationDate','','Title','',...
                 'Abstract','','Authors','','Citation','');

% Loop through each article in hits and extract the PubMed ID,
% publication date, title, abstract, authors, and citation information.
% Place the information in PMSTRUCT, a MATLAB structure array.
for n = 1:numel(hits)
    pmstruct(n).PubMedID         = regexp(hits{n},'(?<=PMID- ).*?(?=\n)','match', 'once');
    pmstruct(n).PublicationDate  = regexp(hits{n},'(?<=DP  - ).*?(?=\n)','match', 'once');
    pmstruct(n).Title            = regexp(hits{n},'(?<=TI  - ).*?(?=PG  -|AB  -)','match', 'once');
    pmstruct(n).Abstract         = regexp(hits{n},'(?<=AB  - ).*?(?=AD  -)','match', 'once');
    pmstruct(n).Authors          = regexp(hits{n},'(?<=AU  - ).*?(?=\n)','match');
    pmstruct(n).Citation         = regexp(hits{n},'(?<=SO  - ).*?(?=\n)','match', 'once');
end


