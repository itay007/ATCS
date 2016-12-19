function [GOanc idx] = getancestors(GO,goid,varargin)
%GETANCESTORS finds ancestors of a GO term.
%
%  GETANCESTORS(GO,ID) returns the numeric IDs of the ancestors of the
%  term ID including ID in the Gene Ontology object GO. ID is a
%  non-negative integer, it can also contain a numeric vector with a set
%  of IDs.
%
%  GETANCESTORS(...,'HEIGHT', L) only searches upwards through L levels of
%  the Gene Ontology. L is a non-negative integer. Default is Inf.
%
%  GETANCESTORS(...,'RELATIONTYPE', T)  traverses the gene ontology by
%  following relationships of type T, options are 'is_a', 'part_of', or
%  'both' (default).
%
%  GETANCESTORS(...,'EXCLUDE', true) excludes the original queried term(s)
%  from the output, unless a term was reached from traversing the graph.
%  Default is false. 
%
%  [ANC, COUNTS] = GETANCESTORS(...) returns the number of times each
%  ancestor is reached. COUNTS is a column vector with the same number of
%  elements as terms in the GO object.
%
%  Example:
%
%   % Find all the ancestors of the 'response to DDT' GO term, display a
%   % report of the ancestors and show their relationships in a figure:
%   GO = geneont('LIVE', true)
%   ancestors = getancestors(GO,46680)
%   subontology = GO(ancestors)
%   rpt = get(subontology.terms,{'id','name'})
%   cm = getmatrix(subontology);
%   BG = biograph(cm, get(subontology.terms, 'name'));
%   view(BG)
%
%  See also GENEONT GENEONT.GENEONT.GETDESCENDANTS
%  GENEONT.GENEONT.GETMATRIX GENEONT.GENEONT.GETRELATIVES
%  GENEONT.GENEONT.SUBSREF GOANNOTREAD NUM2GOID

%   References:
% [1]  http://www.geneontology.org/

%   Copyright 2002-2005 The MathWorks, Inc.


relTypes = {'both','is_a','part_of'};
% set defaults
exclude = false;
height = inf;
relation = 0;

% check input arguments
if nargin > 1
    if rem(nargin,2) ~= 0
        error(message('bioinfo:geneont:getancestors:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'height','exclude','relationtype'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        if(isstruct(pval))
            error(message('bioinfo:geneont:getancestors:StructParamError'));
        end
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:geneont:getancestors:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:geneont:getancestors:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % height
                    if(pval >= 0 && isnumeric(pval))
                        height =  pval;
                    else
                        error(message('bioinfo:geneont:getancestors:IncorrectHeight', upper( char( okargs( k ) ) )));
                    end
                case 2  % exclude
                    exclude = bioinfoprivate.opttf(pval);
                    if isempty(exclude)
                        error(message('bioinfo:geneont:getancestors:ExcludeNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 3 % relation
                    if isnumeric(pval) && isscalar(pval) && pval>=0 && pval<=2
                        relation = round(pval);
                    elseif ischar(pval)
                        relation = strmatch(lower(pval),relTypes)-1;
                        if isempty(relation) 
                            error(message('bioinfo:geneont:getancestors:IncorrectRelationTypeString', pval))
                        end
                    else
                        error(message('bioinfo:geneont:getancestors:IncorrectRelationType'))
                    end
            end
        end
    end
end
if isa(goid,'geneont.term')
    terms = get(goid,'id');
    if iscell(terms)
        terms = cell2mat(terms);
    end
else
    terms = goid(:);
end

% remove any IDs that are invalid
null_terms = GO.hashID(terms+1)==0;
if any(null_terms)
    warning(message('bioinfo:geneont:getancestors:InvalidID', num2str( terms( null_terms ) )));
    terms(null_terms) = [];
end

if isempty(terms)
   GOanc = [];
   if nargout>0
       idx = zeros(size(GO.gonumlist));
   end
   return
end

if nargout>1
    switch relation
        case 0
            cm = double(GO.isarelationmatrix|GO.partofrelationmatrix);
        case 1
            cm = double(GO.isarelationmatrix);
        case 2
            cm = double(GO.partofrelationmatrix);
    end
    n = length(cm);
    buf = accumarray(terms+1,1,[n,1]);
    if exclude
        idx = zeros(n,1);
    else
        idx = buf;
    end
    while nnz(buf) && height>0
        buf = cm * buf;
        idx = idx + buf;
        height = height - 1;
    end
    GOanc = find(idx>0)-1;
    idx = idx(GO.gonumlist+1);
else
    switch relation
        case 0
            cm = GO.isarelationmatrix|GO.partofrelationmatrix;
        case 1
            cm = GO.isarelationmatrix;
        case 2
            cm = GO.partofrelationmatrix;
    end
    n = length(cm);
    buf = accumarray(terms+1,1,[n,1],@all);
    if exclude
        idx = false(n,1);
    else
        idx = buf;
    end
    while any(buf) && height>0
        buf = any(cm(:,buf),2);
        idx = idx | buf;
        height = height - 1;
    end
    GOanc = find(idx)-1;
end

