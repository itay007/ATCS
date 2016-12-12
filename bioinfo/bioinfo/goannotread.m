function A = goannotread(file,varargin)
% GOANNOTREAD parses Gene Ontology Annotated files and returns annotations
%
%   ANNOTATION = GOANNOTREAD(FILE) converts the contents of FILE into an
%   array of structs, ANNOTATION. Files should have the structure specified
%   in http://www.geneontology.org/GO.format.annotation.shtml
%
%   GOANNOTATREAD complies with the GAF 1.0 and GAF 2.0 file formats.
%
%   GOANNOTREAD(...,'FIELDS',FIELDNAMES) allows you to specify the fields
%   to be read from the annotation file. FIELDNAMES can be a single string
%   or a cell array of strings. The default is to read all fields. Valid
%   fieldnames are:
%
%   Database
%	DB_Object_ID
%	DB_Object_Symbol
%	Qualifier
%	GOid
%	DBReference
%	Evidence
%	WithFrom
%	Aspect
%	DB_Object_Name
%	Synonym
%	DB_Object_Type
%	Taxon
%	Date
%	Assigned_by
%   Annotation_Extension
%   Gene_Product_ID
%
%   More detailed information about these fields can be found in the
%   reference page for this function and at
%   http://www.geneontology.org/GO.format.annotation.shtml
%
%   GOANNOTREAD(...,'ASPECT',ASPECT) allows you to specify which aspects
%   are read from the annotation file. ASPECT should be a character array
%   containing one of more characters corresponding to valid aspects. Valid
%   aspects are: P (biological process), F (molecular function) or C
%   (cellular component). The default is to read all aspects ('CFP').
%
%   Many annotation files can be found here:
%   http://www.geneontology.org/GO.current.annotations.shtml
%
%   Example:
%        % Note that the SGD file must first be downloaded and unzipped from
%        % http://www.geneontology.org/GO.current.annotations.shtml.
%        SGDGenes = goannotread('gene_association.sgd');
%        S = struct2cell(SGDGenes);
%        genes = S(3,:)'
%
%   See also CNSGENEEXPDEMO, GENEONT, GENEONT.GENEONT.GETANCESTORS,
%   GENEONT.GENEONT.GETDESCENDANTS, GENEONT.GENEONT.GETMATRIX,
%   GENEONT.GENEONT.GETRELATIVES, GENEONT.GENEONT.SUBSREF,
%   GENEONTOLOGYDEMO, NUM2GOID.

%   Copyright 2005-2012 The MathWorks, Inc.

% GAF 1.0
% COL CONTENT         MANDATORY CARDINALITY        DESCRIPTION
% 1) DB               mandatory  1  - the database contributing the gene_association file
% 2) DB_Object_ID     mandatory  1  - a unique identifier in DB for the item being annotated
% 3) DB_Object_Symbol mandatory  1  - symbol to which DB_Object_ID is matched
% 4) Qualifier 	      optional  >=0 - flags that modify the interpretation of an annotation
% 5) GOid             mandatory  1  - the GO identifier for the term attributed to the DB_Object_ID
% 6) DB:Reference	  mandatory >=1 - the unique identifier appropriate to DB for the authority for the attribution of the GOid to the DB_Object_ID
% 7) Evidence		  mandatory  1  - one of IMP, IGI, IPI, ISS, IDA, IEP, IEA, TAS,NAS, ND, IC, RCA
% 8) With (or) From   optional  >=0 - one of: DB:gene_symbol,DB:gene_symbol[allele_symbol],DB:gene_id,DB:protein_name,DB:sequence_id,GO:GO_id
% 9) Aspect			  mandatory  1  - one of P (biological process), F (molecular function) or C (cellular component)
%10) DB_Object_Name   optional  0,1 - name of gene or gene product  [white space allowed]
%11) Synonym          optional  >=0 - Gene_symbol
%12) DB_Object_Type	  mandatory  1  - one of gene, transcript, protein,protein_structure, complex
%13) Taxon			  mandatory 1,2 - taxonomic identifier(s)
%14) Date			  mandatory  1  - Date on which the annotation was made; format is YYYYMMDD
%15) Assigned_by	  mandatory  1  - The database which made the annotation

% GAF 2.0 adds two extra columns
%16) Annotation_Extension optional >=0 - References to other ontologies that enhance the annotation
%17) Gene_Product_ID      optional 0,1 - Unique identifier to variants of the gene or gene product

bioinfochecknargin(nargin,1,mfilename);
[fields,aspect] = parse_inputs(varargin{:});

fid = fopen(file);
if(fid == -1)
    error(message('bioinfo:goannotread:CannotOpenFile', file));
end
c = onCleanup(@()fclose(fid));

%Skip the header text and only keep the gaf-version (if found)
pos = 0;
gafver = 0;
line = fgetl(fid);
while (~isempty(regexp(line,'^!', 'once')) || ~isempty(regexp(line,'^#', 'once')))
    if numel(line)>1
        vs = regexpi(line(2:end),'^\s*gaf-version:\s*([\d\.]+)\s*','tokens','once');
        if ~isempty(vs)
            gafverFile = str2double(vs{1});
            gafver = floor(gafverFile); % this reader only differentiates between major versions
        end
    end
    pos = ftell(fid);
    line = fgetl(fid);
end

if gafver>=3
    error(message('bioinfo:goannotread:UnsupportedGAFversion',sprintf('%3.1f',gafverFile)))
end

%--% Check the validity of the format by looking at the number of columns
%    in the first line:

firstline = strread(line,'%s',1,'delimiter','\n');
columncount = length(regexp(firstline{1},'(.*?\t)')) + 1;

% Old formats may contain 14 or 7 columns, if this is the case gafver
% should have not been detected in the header and we still try to read them
% for backwards compatibility.

% For GAF 1.0 files we must count 15 columns (last column is mandatory)
if (gafver==1) && (columncount~=15)
    error(message('bioinfo:goannotread:CorruptedFileGAF10'))
end

% For GAF 2.0 files we must count 15,16 or 17 columns (last two columns are
% optional)
if (gafver==2) && ~any(columncount == [15 16 17])
    error(message('bioinfo:goannotread:CorruptedFileGAF20'))
end

% Anyfile that is unversioned (Old formats) and does not have 7,14,15,16
% columns cannot be read:
% Note: 16 is accepted because there are old unversioned files with all
% columns as in GAF 1.0 but with an extra \t at the end of the line.
if gafver==0 && ~any(columncount == [7 14 15 16])
    error(message('bioinfo:goannotread:NonAnnotationFile', file));
end

%--% Check that filtering options are in accordance with the GAF version
allAspects = isequal(aspect.validAspects,aspect.selectedAspects);
if (columncount==7)
    f = {'TaxonomyID','GeneID','GOid','Evidence','Qualifier','GOterm','PubMed'};
    goidFieldNum = 3;
    aspectFieldNum = 0; % aspect does not exist in this old format
    formatString = '%s%s%s%s%s%s%s';
    % For old formats no filtering is allowed:
    if ~allAspects || ~strcmpi(fields,'all')
        error(message('bioinfo:goannotread:FilteringNotAllowed'))
    end
else
    if (gafver==2) % GAF 2.0
        f = {'Database','DB_Object_ID','DB_Object_Symbol','Qualifier','GOid','DBReference','Evidence','WithFrom','Aspect','DB_Object_Name','Synonym','DB_Object_Type','Taxon','Date','Assigned_by','Annotation_Extension','Gene_Product_Form_ID'};
    else % GAF 1.0 and old 14|15 column format
        f = {'Database','DB_Object_ID','DB_Object_Symbol','Qualifier','GOid','DBReference','Evidence','WithFrom','Aspect','DB_Object_Name','Synonym','DB_Object_Type','Taxon','Date','Assigned_by'};
    end
    goidFieldNum = 5;
    aspectFieldNum = 9;
    returnAspect = true;
    if strcmpi(fields,'all')
        matches = true(size(f));
    else
        [~,loc] = ismember(lower(fields),lower(f));
        if any(loc == 0)
            badField = (loc == 0);
            error(message('bioinfo:goannotread:UnknownField', fields{ badField }));
        end
        matches = ismember(lower(f),lower(fields));
        if matches(goidFieldNum)
            goidFieldNum = sum(matches(1:goidFieldNum));
        else
            goidFieldNum = 0;
        end
        if ~allAspects
            if ~matches(aspectFieldNum)
                returnAspect = false;
            end
            matches(aspectFieldNum) = true;
            aspectFieldNum = sum(matches(1:aspectFieldNum));
        end
    end
    formatString = '';
    for counter = 1:numel(f)
        if matches(counter)
            formatString = [formatString '%s']; %#ok<AGROW>
        else
            formatString = [formatString '%*s']; %#ok<AGROW>
        end
    end
    
    f = f(matches);
    if ~returnAspect
        f(aspectFieldNum) = [];
    end
end

%-% Read -ALL- data from file
fseek(fid,pos,'bof');
theData = textscan(fid,formatString, 'delimiter', '\t');
d = [theData{:}];

%-% Filter annotations to unwanted aspects
if ~allAspects
    mask = ~ismember(upper(d(:,aspectFieldNum)),cellstr(aspect.selectedAspects'));
    if ~returnAspect % If aspect is not wanted in the output remove column
        d(:,aspectFieldNum) = [];
    end
    d(mask,:) =[];
end

%-% Change GO ids to numbers
if goidFieldNum >0
    goids = char(d(:,goidFieldNum));
    d(:,goidFieldNum) = num2cell(str2num(goids(:,4:10))); %#ok
end

%-% Output is a structure
A = cell2struct(d,f,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fields,aspect] = parse_inputs(varargin)
% Parse the varargin parameter/value inputs

% Check that we have the right number of inputs
if rem(nargin,2)== 1
    error(message('bioinfo:goannotread:IncorrectNumberOfArguments', mfilename));
end

% The allowed inputs
okargs = {'fields','aspect'};

% Set default values
fields = 'all';
aspect.validAspects  = 'CFP';
aspect.selectedAspects = aspect.validAspects;
% Loop over the values
for j=1:2:nargin
    % Lookup the pair
    [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % Fields -- check these later when we know the type of file
            if ~ischar(pval) && ~iscellstr(pval)
                error(message('bioinfo:goannotread:FieldsNotString'))
            end
            if isempty(pval)
                error(message('bioinfo:goannotread:EmptyField'))
            end
            if ischar(pval)
                pval = {pval};
            end
            fields = sort(upper(pval));
        case 2  % Aspect
            if ~ischar(pval)
                error(message('bioinfo:goannotread:AspectNotChar'))
            end
            aspect.selectedAspects = sort(upper(pval));
            if isempty(intersect(aspect.selectedAspects,aspect.validAspects)) ||...
                    ~isempty(setdiff(aspect.selectedAspects,aspect.validAspects))
                error(message('bioinfo:goannotread:InvalidAspect'))
            end
    end
end
