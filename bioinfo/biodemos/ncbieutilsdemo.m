%% Accessing NCBI Entrez Databases with E-Utilities
% This example shows how to programmatically search and retrieve data from
% NCBI's Entrez databases using NCBI's Entrez Utilities (E-Utilities) with
% |urlread| and |urlwrite|.

%   Copyright 2007-2012 The MathWorks, Inc.


%%  Using NCBI E-Utilities to Retrieve Biological Data 
% E-Utilities (eUtils) are eight server-side programs (e.g. ESearch,
% ESummary, EFetch, etc.,) developed and maintained by NCBI for searching
% and retrieving data from most Entrez Databases. You access tools via URLs
% with a strict syntax of a specific base URL, a call to the eUtil's script
% and its associated parameters. For more details on eUtils, see
% <http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html
% E-Utilities Help>.

%% Searching Nucleotide Database with ESearch
% If you are interested in studying the H5N1 virus you might start your
% investigation by analyzing the genes sequenced from this virus.  The H5N1
% virus isolated in 1997 from a chicken in Hong Kong (A/Chicken/Hong
% Kong/915/97(H5N1)) could be a good starting point for your analysis. This
% particular virus jumped from chickens to humans, killing six people before
% the spread of the disease was brought under control by destroying all
% poultry in Hong Kong [1]. You can use ESearch to find the sequence data
% needed for your research. ESearch requires input of a database (|db|) and
% search term (|term|). Optionally, you can request for ESearch to store
% your search results on the NCBI history server through the |usehistory|
% parameter. For additional information on ESearch, see
% <http://eutils.ncbi.nlm.nih.gov/entrez/query/static/esearch_help.html
% ESearch Help>.   

baseURL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
eutil = 'esearch.fcgi?';
dbParam = 'db=nuccore';
termParam = '&term=A/chicken/Hong+Kong/915/97+OR+A/chicken/Hong+Kong/915/1997';
usehistoryParam = '&usehistory=y';
esearchURL = [baseURL, eutil, dbParam, termParam, usehistoryParam]

%%
% The |term| parameter can be any valid Entrez query (see
% <http://www.ncbi.nlm.nih.gov/books/bv.fcgi?rid=helpentrez.chapter.EntrezHelp
% Entrez Help>).  There cannot be any spaces in the URL, so parameters are
% separated by '&' and any spaces in a query term need to be replaced with
% '+' (e.g. 'Hong+Kong'). 
%
% You can use |urlread| to send the URL and return the results from ESearch
% as a character array. 

searchReport = urlread(esearchURL)

%%
% ESearch returns the search results in XML.  The report contains
% information about the query performed, which database was searched and
% UIDs (unique IDs) to the records that match the query.  If you use the
% history server, the report contains two additional IDs, |WebEnv| and
% |query_key|, for accessing the results. |WebEnv| is the location of the
% results on the server, and |query_key| is a number indexing the queries
% performed.  Since WebEnv and |query_key| are query dependent they will
% change every time the search is executed. Either the UIDs or |WebEnv| and
% |query_key| can be parsed out of the XML report then passed to other
% eUtils.  You can use |regexp| to do the parsing and store the tokens in
% the structure with fieldnames |WebEnv| and |QueryKey|.

ncbi = regexp(searchReport,...
    '<QueryKey>(?<QueryKey>\w+)</QueryKey>\s*<WebEnv>(?<WebEnv>\S+)</WebEnv>',...
    'names')

%% Getting GenBank(R) File Summaries with ESummary
% To get a quick overview of sequences that matched the query you can use
% ESummary.  ESummary retrieves a brief summary, or Document Summary
% (DocSum), for each record.  ESummary requires an input of which database
% to access and which records to retrieve, identified either by a list of
% UIDs passed through |id| parameter or by the |WebEnv| and |query_key|
% parameters.  ESummary returns a report in XML that contains the summary
% information for each record. See
% <http://eutils.ncbi.nlm.nih.gov/entrez/query/static/esummary_help.html
% ESummary Help> for additional information.  Use |urlwrite| with ESummary
% to perform the record summary retrieval and write out the XML report to a
% file.

tmpDirectory = tempdir;
summaryFname = fullfile(tmpDirectory,'summaryReport.xml');

urlwrite([baseURL...
    'esummary.fcgi?db=nuccore&WebEnv=',ncbi.WebEnv,...
    '&query_key=',ncbi.QueryKey],summaryFname);

%%
% You can create an XSL stylesheet to view information from the ESummary
% XML report in a web browser.  For more information on writing XSL
% stylesheets, see <http://www.w3.org/Style/XSL W3C(R) XSL>.  An XSL
% stylesheet was created for this example to view the sequence summary
% information and provide links to their full GenBank(R) files.  |Xslt| can
% be used to view the XML report in a Web browser from MATLAB(R).

xslt(summaryFname,'genbankSummary.xsl','-web');
%%
% <<ncbieutilsdemo_a.png>>
%% Retrieving Full GenBank Files with EFetch
% To perform the sequence analysis, you need to get the full
% GenBank record for each sequence.  EFetch retrieves full records from
% Entrez databases.  EFetch requires an input of a database and a list of
% UIDs or |WebEnv| and |query_key|.  Additionally, EFetch can return the
% output in different formats.  You can specify which output format (i.e.
% GenBank (|gb|), FASTA) and file format (i.e. text, ASN.1, XML) you want
% through the |rettype|  and |retmode|  parameters, respectively.  |Rettype|
% equals |gb| for GenBank file format and |retmode| equals |text| for this
% query.  |Genbankread| can be used directly with the EFetch URL to
% retrieve all the GenBank records and read them into a structure array.
% This structure can then be used as input to |seqviewer| to visualize the
% sequences. 
 
ch97struct = genbankread([baseURL...
    'efetch.fcgi?db=nuccore&rettype=gb&retmode=text&WebEnv=',ncbi.WebEnv,...
    '&query_key=',ncbi.QueryKey]);
seqviewer(ch97struct)
%%
% <<ncbieutilsdemo_b.png>>
%% Finding Links Between Databases with ELink
% It might be useful to have PubMed articles related to these genes
% records.  ELink provides this functionality.  It finds associations 
% between records within or between databases.  You can give ELink the
% query_key and WebEnv IDs from above and tell it to find records in the
% PubMed Database (|db| parameter) associated with your records from the
% Nucleotide (nuccore) Database (|dbfrom| parameter).  ELink returns an XML
% report with the UIDs for the records in PubMed.  These UIDs can be parsed
% out of the report and passed to other eUtils (e.g. ESummary). Use the
% stylesheet created for viewing ESummary reports to view the results of
% ELink.

elinkReport = urlread([baseURL...
    'elink.fcgi?dbfrom=nuccore&db=pubmed&WebEnv=', ncbi.WebEnv,...
    '&query_key=',ncbi.QueryKey]);

%%
% Extract the PubMed UIDs from the ELink report.
pubmedIDs = regexp(elinkReport,'<Link>\s+<Id>(\w*)</Id>\s+</Link>','tokens');
NumberOfArticles = numel(pubmedIDs)

%Put PubMed UIDs into a string that can be read by EPost URL.
pubmed_str = [];
for ii = 1:NumberOfArticles
    pubmed_str = sprintf([pubmed_str '%s,'],char(pubmedIDs{ii}));
end

%% Posting UIDs to NCBI History Server with EPost
% You can use EPost to posts UIDs to the history server.  It returns an XML
% report with a |query_key| and |WebEnv| IDs pointing to the location of
% the history server.  Again, these can be parsed out of the report and
% used with other eUtils calls.    

epostReport = urlread([baseURL 'epost.fcgi?db=pubmed&id=',pubmed_str(1:end-1)]);
epostKeys = regexp(epostReport,...
    '<QueryKey>(?<QueryKey>\w+)</QueryKey>\s*<WebEnv>(?<WebEnv>\S+)</WebEnv>','names')

%% Using ELink to Find Associated Files within the Same Database
% ELink can do "within" database searches.  For example, you can query for
% a nucleotide sequence within Nucleotide (nuccore) database to find
% similar sequences, essentially performing a BLAST search.  For "within"
% database searches, ELink returns an XML report containing the related
% records, along with a score ranking its relationship to the query record.
% From the above PubMed search, you might be interested in finding all
% articles related to those articles in PubMed.  This is easy to do with
% ELink. To do a "within" database search, set |db| and |dbfrom| to PubMed.
% You can use the |query_key| and |WebEnv| from the EPost call.

pm2pmReport = urlread([baseURL...
    'elink.fcgi?dbfrom=pubmed&db=pubmed&query_key=',epostKeys.QueryKey,...
    '&WebEnv=',epostKeys.WebEnv]);
pubmedIDs = regexp(pm2pmReport,'(?<=<Id>)\w*(?=</Id>)','match');
NumberOfArticles = numel(unique(pubmedIDs))

pubmed_str = [];
for ii = 1:NumberOfArticles
    pubmed_str = sprintf([pubmed_str '%s,'],char(pubmedIDs{ii}));
end

%%
% Use |urlwrite| with EFetch to retrieve full abstracts for all the
% articles and write out the returned XML report to a file. An XSL
% stylesheet was created for and is provided with this example for viewing
% the results of the EFetch query. The XML report can be transformed using
% the stylesheet and opened in a Web browser from MATLAB using |xslt|.

fullFname = fullfile(tmpDirectory,'H5N1_relatedArticles.xml');
urlwrite([baseURL 'efetch.fcgi?db=pubmed&retmode=xml&id=',...
    pubmed_str(1:end-1)],fullFname);
xslt(fullFname,'pubmedFullReport.xsl','-web');
%%
% <<ncbieutilsdemo_c.png>>
%% Using EGQuery to get a Global View of H5N1 Related-Records in Entrez
% To see what other Entrez databases contain information about the H5N1
% virus, use EGQuery.  EGQuery performs a text search across all available
% Entrez databases and returns the number of hits in each.  EGQuery accepts
% any valid Entrez text query as input through the |term| parameter.  

entrezSearch = urlread([baseURL,'egquery.fcgi?term=H5N1+AND+virus']);
entrezResults = regexp(entrezSearch,...
    '<MenuName>(?<DB>\w+\s*\w*)</MenuName>\s*<Count>(?<Count>\d+)</Count>',...
    'names');
entrezDBs = {entrezResults(:).DB};
dbCounts = str2double({entrezResults(:).Count});
%remove databases with no records
entrezDBs = entrezDBs(logical(dbCounts));
[dbCounts,sortInd] = sort(dbCounts(logical(dbCounts)));
entrezDBs = entrezDBs(sortInd);
numDBs = numel(entrezDBs);
barh(log10(dbCounts));
set(gca,'YLim',[.5 numDBs+.5]);
set(gca,'YTick',1:numDBs);
set(gca,'YTickLabel',entrezDBs);
xlabel('Log(Number of Records)');
title('Number of H5N1 Related-Records Per Entrez Database');

%% References
% [1] Cristianini,N., Hahn, M.W., Introduction to Computational Genomics :
% A Case Studies Approach, Cambridge University Press, 2007

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20NCBIEUTILSDEMO%20in%20Bioinformatics%20Toolbox%204.2 Provide feedback for this example.>*

displayEndOfDemoMessage(mfilename)



