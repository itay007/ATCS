function [GOdec idx] = getdescendants(GO,goid,varargin)
%GETDESCENDANTS finds descendants of a GO term.
%
%  GETDESCENDANTS(GO,ID) returns the numeric IDs of the descendants of the
%  term ID including ID in the Gene Ontology object GO. ID is a
%  non-negative integer, it can also contain a numeric vector with a set
%  of IDs.
%
%  GETDESCENDANTS(...,'DEPTH', L) only searches downwards through L levels
%  of the Gene Ontology. L is a non-negative integer. Default is Inf.
%
%  GETDESCENDANTS(...,'RELATIONTYPE', T) traverses the gene ontology by
%  following relationships of type T, options are 'is_a', 'part_of', or
%  'both' (default).
%
%  GETDESCENDANTS(...,'EXCLUDE', true) excludes the original queried
%  term(s) from the output, unless a term was reached from traversing the
%  graph. Default is false. 
%
%  [DEC, COUNTS] = GETDESCENDANTS(...) returns the number of times each
%  descendant is reached. COUNTS is a column vector with the same number
%  elements as terms in the GO object.
%
%  Example:
%
%   % Find all the descendants of the 'aldo-keto reductase activity' GO term,
%   % display a report of the descendants and show their relationships in a
%   % figure:
%   GO = geneont('LIVE', true)
%   descendants = getdescendants(GO,4033)
%   subontology = GO(descendants)
%   rpt = [num2goid(cell2mat(get(subontology.terms,'id')))...
%          get(subontology.terms,'name')]';
%   disp(sprintf('%s --> %s \n',rpt{:}))
%   cm = getmatrix(subontology);
%   BG = biograph(cm, rpt(1,:));
%   view(BG)
%
%  See also GENEONT GENEONT.GENEONT.GETANCESTORS
%  GENEONT.GENEONT.GETMATRIX GENEONT.GENEONT.GETRELATIVES
%  GENEONT.GENEONT.SUBSREF GOANNOTREAD NUM2GOID

%   References:
% [1]  http://www.geneontology.org/

%   Copyright 2002-2005 The MathWorks, Inc.


relTypes = {'both','is_a','part_of'};
% set defaults
exclude = false;
depth = inf;
relation = 0;

% check  input arguments
if nargin > 1
    if rem(nargin,2) ~= 0
        error(message('bioinfo:geneont:getdescendants:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'depth','exclude','relationtype'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        if(isstruct(pval))
            error(message('bioinfo:geneont:getdescendants:StructParamError'));
        end
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:geneont:getdescendants:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:geneont:getdescendants:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % depth
                    if(pval >= 0 && isnumeric(pval))
                        depth =  pval;
                    else
                        error(message('bioinfo:geneont:getdescendants:IncorrectHeight', upper( char( okargs( k ) ) )));
                    end
                case 2  % exclude
                    exclude = bioinfoprivate.opttf(pval);
                    if isempty(exclude)
                        error(message('bioinfo:geneont:getdescendants:ExcludeNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 3 % relation
                    if isnumeric(pval) && isscalar(pval) && pval>=0 && pval<=2
                        relation = round(pval);
                    elseif ischar(pval)
                        relation = strmatch(lower(pval),relTypes)-1; %#ok
                        if isempty(relation) 
                            error(message('bioinfo:geneont:getdescendants:IncorrectRelationTypeString', pval))
                        end
                    else
                        error(message('bioinfo:geneont:getdescendants:IncorrectRelationType'))
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
    warning(message('bioinfo:geneont:getdescendants:InvalidID', num2str( terms( null_terms ) )));
    terms(null_terms) = [];
end

if isempty(terms)
   GOdec = [];
   if nargout>0
       idx = zeros(size(GO.gonumlist));
   end
   return
end

if nargout>1
    switch relation
        case 0
            cm = double(GO.isarelationmatrixi|GO.partofrelationmatrixi);
        case 1
            cm = double(GO.isarelationmatrixi);
        case 2
            cm = double(GO.partofrelationmatrixi);
    end
    n = length(cm);
    buf = accumarray(terms+1,1,[n,1]);
    if exclude
        idx = zeros(n,1);
    else
        idx = buf;
    end
    while nnz(buf) && depth>0
        buf = cm * buf;
        idx = idx + buf;
        depth = depth - 1;
    end
    GOdec = find(idx>0)-1;
    idx = idx(GO.gonumlist+1);
else
    switch relation
        case 0
            cm = GO.isarelationmatrixi|GO.partofrelationmatrixi;
        case 1
            cm = GO.isarelationmatrixi;
        case 2
            cm = GO.partofrelationmatrixi;
    end
    n = length(cm);
    buf = accumarray(terms+1,1,[n,1],@all);
    if exclude
        idx = false(n,1);
    else
        idx = buf;
    end
    while any(buf) && depth>0
        buf = any(cm(:,buf),2);
        idx = idx | buf;
        depth = depth - 1;
    end
    GOdec = find(idx)-1;
end

