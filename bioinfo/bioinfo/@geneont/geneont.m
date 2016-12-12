function hh = geneont(varargin)
% GENEONT creates a Gene Ontology object from an 'OBO' formatted file.
%
%   GO = GENEONT('FILE', FILENAME) uses an 'OBO' formatted file to create
%   the Gene Ontology object.
%
%   GO = GENEONT without input arguments will search for the file
%   'gene_ontology.obo' in the current directory.
%
%   GO = GENEONT('LIVE', true) creates the Gene Ontology object from
%   http://www.geneontology.org/ontology/gene_ontology.obo which is the
%   most up to date version of the Gene Ontology Database.
%
%   GO = GENEONT('LIVE', true,'TOFILE', FILENAME) saves the contents of
%   http://www.geneontology.org/ontology/gene_ontology.obo to a local file.
%
%   Note that running this function will take a few minutes for the full
%   Gene Ontology Database from http://www.geneontology.org.
%
%   Example:
%       % Read the latest version of gene_ontology.obo from the Web
%       GO = geneont('LIVE',true);
%
%   See also GENEONT.GENEONT.GETANCESTORS, GENEONT.GENEONT.GETDESCENDANTS,
%   GENEONT.GENEONT.GETMATRIX, GENEONT.GENEONT.GETRELATIVES,
%   GENEONT.GENEONT.SUBSREF, GENEONTOLOGYDEMO, GOANNOTREAD, NUM2GOID.

%   References:
% [1]  http://www.geneontology.org/

%   Copyright 2002-2012 The MathWorks, Inc.


toFile = false;
live = false;
file = '';
%%% check arguments
if nargin > 0
    if rem(nargin,2) ~= 0
        error(message('bioinfo:geneont:geneont:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'live', 'file','tofile'};
    for j=1:2:nargin
        pname = varargin{j};
        pval = varargin{j+1};
        if(isstruct(pval))
            error(message('bioinfo:geneont:geneont:StructParamError'));
        end
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:geneont:geneont:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:geneont:geneont:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % live
                    live = bioinfoprivate.opttf(pval);
                    if isempty(live)
                        error(message('bioinfo:geneont:geneont:LiveNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 2  % file
                    if (ischar(pval) && size(pval,1)==1 && (size(fullfile(pwd,pval), 2) < 1000) && (exist(pval,'file') || exist(fullfile(pwd,pval),'file')))
                        file = pval;
                    else
                        error(message('bioinfo:geneont:geneont:FileNotFound', pval));
                    end
                case 3  % tofile
                    toFile = true;
                    saveFile = pval;
            end
        end
    end
end

if(live == true)
    try
        theURL = 'http://www.geneontology.org/ontology/gene_ontology.obo';
        source = theURL;
        obo = urlread(theURL);
    catch theException
        error(message('bioinfo:geneont:geneont:UrlreadProblem', theException.message));
    end
    if isempty(obo)
        error(message('bioinfo:geneont:geneont:EmptyOBOFromWeb', theURL));
    end
else
    if(~isempty(file))
        source = file;
        fid = fopen(file);
        if(fid == -1)
            error(message('bioinfo:geneont:geneont:CannotOpenFile', file));
        end
    else
        source = 'gene_ontology.obo';
        fid = fopen('gene_ontology.obo');
        if(fid == -1)
            error(message('bioinfo:geneont:geneont:CannotOpenDefaultFile'));
        end
    end
    obo = fread(fid,'*char')';
    fclose(fid);
end
if(toFile == true)
    fid = fopen(saveFile,'wt');
    if(fid == -1)
        warning(message('bioinfo:geneont:geneont:CannotOpenSaveFile', saveFile));
    else
        fwrite(fid,obo);
        fclose(fid);
    end
end

% make each line a cell string
obo = strread(obo,'%s','delimiter','\n','whitespace','');
obo(cellfun('isempty',obo))=[];
rowsWithTokens = find(cellfun('length',regexp(obo,'^\[.+\]')));
terms = strcmp(obo(rowsWithTokens),'[Term]');
typedef = strcmp(obo(rowsWithTokens),'[Typedef]');
rowsWithTerms = rowsWithTokens(terms);
rowsWithTypedef = rowsWithTokens(typedef);
lastRowsforTokens = [rowsWithTokens(2:end)-1;numel(obo)];
lastRowsforTerms = lastRowsforTokens(terms);
lastRowsforTypedef = lastRowsforTokens(typedef);

preamble = sprintf('%s\n',obo{1:rowsWithTokens(1)-1});
hk = hashkeywords;
keywords = full(hk(max(0,char(regexp(obo,'^\w{1,6}','match','once'))-'`')*[27^5;27^4;27^3;27^2;27;1]+1));
obo = regexprep(obo,'^\w*\: ',''); % once I have the keywords, remove them from obo

regionWithTerms = zeros(numel(obo)+1,1);
regionWithTerms(rowsWithTerms,1) = 1;
regionWithTerms(lastRowsforTerms+1,1) = regionWithTerms(lastRowsforTerms+1,1)-1;
regionWithTerms = cumsum(regionWithTerms(1:end-1));

regionWithTypedef = zeros(numel(obo)+1,1);
regionWithTypedef(rowsWithTypedef,1) = 1;
regionWithTypedef(lastRowsforTypedef+1,1) = regionWithTypedef(lastRowsforTypedef+1,1)-1;
regionWithTypedef = cumsum(regionWithTypedef(1:end-1));

GOacc = sscanf(char(obo(regionWithTerms & keywords==1))','GO:%d');
if numel(unique(GOacc))~=numel(GOacc)
    error(message('bioinfo:geneont:geneont:IDnotUnique', source))
end

IDregions =  zeros(numel(obo),1);
IDregions(rowsWithTerms+1,1) = GOacc;
IDregions(lastRowsforTerms+1,1) = -GOacc;
IDregions = cumsum(IDregions);

% build a vector with the relation types
relaTypes = (regionWithTerms & keywords==13)+0;
relList = ['is_a' obo(regionWithTypedef & keywords==1)'];
idx = regionWithTerms & keywords==14;
relaTypes(idx) = seqmatch(regexp(obo(idx),'^\w+','match','once'),relList);

% find relationships (now we only do is_a and part_of)
fromID1 = sscanf(char(regexp(obo(relaTypes==1),'GO:\d+','match','once'))','GO:%d')+1;
toID1 = IDregions(relaTypes==1)+1;
pof = find(strcmp(relList,'part_of'));
fromID2 = sscanf(char(regexp(obo(relaTypes==pof),'GO:\d+','match','once'))','GO:%d')+1;
toID2 = IDregions(relaTypes==pof)+1;

hh = geneont.geneont;  % create gene ontology object

% insert header info to the object
A = regexp(preamble, '\n?(.*?)\:\s(.*?)\n', 'tokens');
for m = 1:length(A)
    val = A{m}{2};
    switch A{m}{1}
        case {'format-version'} % format-version
            hh.format_version = val;
        case {'data-version'} % data-version
            hh.data_version = val;  
        case {'version'} % version
            hh.version = val;          
        case {'date'} % date
            hh.date = val;  
        case {'saved-by'} % saved-by
            hh.saved_by = val;
        case {'auto-generated-by'} % auto-generated-by
            hh.auto_generated_by = val;
        case {'subsetdef'} % subsetdef
            hh.subsetdef = [hh.subsetdef; cellstr(val)] ;
        case {'import'} % import
            hh.import = val;
        case {'synonymtypedef'} %synonymtypedef
            hh.synonymtypedef = val;
        case {'idspace'} % idspace
            hh.idspace = val;    
        case {'default-relationship-id-prefix'} % default-relationship-id-prefix
            hh.default_relationship_id_prefix = val;
        case {'id-mapping'} % id-mapping
            hh.id_mapping = val; 
        case {'remark'} % remark
            hh.remark = val;            
        case {'default-namespace'} % default-namespace
            hh.default_namespace = val;
        case {'typeref'} % typeref
            hh.typeref = val;
        otherwise %Handling Unrecognized Tags -> record the unrecognized tag for later serialization (recommended)
            hh.unrecognized_tag = [hh.Unrecognized_Tag;A{m}];
%             warning('bioinfo:geneont:geneont:UnknownField',...
%                 'Unknown field %s with value %s for header information', A{m}{1},val);
    end
end
hh.version = char(regexp(preamble, 'Revision:\s([\.\d]+)','tokens','once'));

m = GOacc(end)+1; % 10000000 is the maximum number of ID a GO can have

isarelationmatrix = sparse(fromID1,toID1,true,m,m);
isarelationmatrixi = sparse(toID1,fromID1,true,m,m);
hh.isarelationmatrix = isarelationmatrix;
hh.isarelationmatrixi = isarelationmatrixi;

partofrelationmatrix = sparse(fromID2,toID2,true,m,m);
partofrelationmatrixi = sparse(toID2,fromID2,true,m,m);
hh.partofrelationmatrix = partofrelationmatrix;
hh.partofrelationmatrixi = partofrelationmatrixi;

hh.hashID = sparse(GOacc+1,1,1:numel(GOacc),10000000,1);
hh.relationtype = relList;
hh.gonumlist = GOacc;

ontTypes =  get(findtype('OntologyTypes'),'Strings');
nT = length(GOacc);
T = handle(nan(nT,1));
for i = nT:-1:1
    u = rowsWithTerms(i)+2:lastRowsforTerms(i);
    data = obo(u);
    keys = keywords(u);
    T(i) = geneont.term(hh);
    h = T(i);
    id = GOacc(i);
    h.id = id;
    h.is_a = find(isarelationmatrix(:,id+1))-1;
    h.part_of = find(partofrelationmatrix(:,id+1))-1;
    h.is_a_i = find(isarelationmatrixi(:,id+1))-1;
    h.part_of_i = find(partofrelationmatrixi(:,id+1))-1;
    for m = 1:length(data)
        switch keys(m)
            case 2 % name
                h.name = data{m};
            case 3 % namespace
                q = find(strncmp(data{m},ontTypes ,3));
                if isempty(q) % unknown
                    q = 1;
                end
                switch q
                    case 1
                        h.ontology = 'unknown';
                    case 2
                        h.ontology = 'cellular component';
                    case 3
                        h.ontology = 'molecular function';
                    case 4
                        h.ontology = 'biological process';
                end
            case 4 % def
                if(isempty(h.definition))
                    h.definition = data{m};
                else
                    warning(message('bioinfo:geneont:geneont:noMoreThanOneDefPerTerm'));
                end
            case 5 % comment
                if(isempty(h.comment))
                    h.comment = data{m};
                else
                    warning(message('bioinfo:geneont:geneont:noMoreThanOneCommentPerTerm'));
                end
            case 6 % alt_id
                h.synonym = [h.synonym;{'alt_id',data{m}}];
            case 7 % synonym
                h.synonym = [h.synonym;{'synonym',data{m}}];
            case 8 % related_synonym
                h.synonym = [h.synonym;{'related_synonym',data{m}}];
            case 9 % exact_synonym
                h.synonym = [h.synonym;{'exact_synonym',data{m}}];
            case 10 % narrow_synonym
                h.synonym = [h.synonym;{'narrow_synonym',data{m}}];
            case 11 % broad_synonym
                h.synonym = [h.synonym;{'broad_synonym',data{m}}];
            case 12 % is_obsolete
                h.obsolete =  true;
        end
    end
end
hh.Terms = T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hk = hashkeywords
% HASHKEYWORDS prepares a hash table for keywords (it uses the ASCII code
% of the first six letters of every keyword)
powers = [27^5;27^4;27^3;27^2;27;1];
hk = sparse([],[],[],27^6,1);
hk(max(0,'id    '-'`')*powers+1) = 1;
hk(max(0,'name  '-'`')*powers+1) = 2;
hk(max(0,'namesp'-'`')*powers+1) = 3;
hk(max(0,'def   '-'`')*powers+1) = 4;
hk(max(0,'commen'-'`')*powers+1) = 5;
hk(max(0,'alt_id'-'`')*powers+1) = 6;
hk(max(0,'synony'-'`')*powers+1) = 7;
hk(max(0,'relate'-'`')*powers+1) = 8;
hk(max(0,'exact_'-'`')*powers+1) = 9;
hk(max(0,'narrow'-'`')*powers+1) = 10;
hk(max(0,'broad_'-'`')*powers+1) = 11;
hk(max(0,'is_obs'-'`')*powers+1) = 12;
hk(max(0,'is_a  '-'`')*powers+1) = 13;
hk(max(0,'relati'-'`')*powers+1) = 14;
hk(max(0,'subset'-'`')*powers+1) = 15;
hk(max(0,'xref_a'-'`')*powers+1) = 15;
hk(max(0,'xref_u'-'`')*powers+1) = 15;
hk(max(0,'use_te'-'`')*powers+1) = 15;
hk(max(0,'domain'-'`')*powers+1) = 15;
hk(max(0,'range '-'`')*powers+1) = 15;
hk(max(0,'is_cyc'-'`')*powers+1) = 15;
hk(max(0,'is_tra'-'`')*powers+1) = 15;
hk(max(0,'is_sym'-'`')*powers+1) = 15;


  


