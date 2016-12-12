function [GOind idx] = getrelatives(GO,terms,varargin)
%GETRELATIVES finds related terms to a GO term.
%
%  GETRELATIVES(GO,ID) returns the numeric IDs of the relatives of the term
%  ID including ID in the Gene Ontology object GO. ID is a non-negative
%  integer, it can also contain a numeric vector with a set of IDs. 
%
%  GETRELATIVES(...,'HEIGHT', H) includes terms that are related up to H
%  levels. H is a non-negative integer. Default is 1.
%
%  GETRELATIVES(...,'DEPTH', D) includes terms that are related down to D
%  levels. D is a non-negative integer. Default is 1.
%
%  GETRELATIVES(...,'LEVELS', L) includes terms that are related both, up
%  and down to L levels. L is a non-negative integer. When specified it
%  overrides 'HEIGHT' and 'DEPTH' values.
%
%  GETRELATIVES(...,'RELATIONTYPE', T)  traverses the gene ontology by
%  following relationships of type T, options are 'is_a', 'part_of', or
%  'both' (default).
%
%  GETRELATIVES(...,'EXCLUDE', true) excludes the original queried term(s)
%  from the output, unless a term was reached from traversing the graph.
%  Default is false.  
%
%  [REL, COUNTS] = GETRELATIVES(...) returns the number of times each
%  relative is reached. COUNTS is a column vector with the same number of
%  elements as terms in the GO object.
%
%  Example: 
%
%   % Find all the relatives of the 'mitochondrial membrane' GO term,
%   % display a report of the relatives and show their relationships in a
%   % figure:
%   GO = geneont('LIVE', true)
%   relatives = getrelatives(GO,31966,'levels',1)
%   subontology = GO(relatives)
%   rpt = get(subontology.terms,{'id','name'})
%   [cm acc rels] = getmatrix(subontology);
%   BG = biograph(cm, get(subontology.terms, 'name'));
%   BG.nodes(acc==31966).Color = [1 0 0];
%   view(BG)
%
%  See also GENEONT GENEONT.GENEONT.GETANCESTORS
%  GENEONT.GENEONT.GETDESCENDANTS GENEONT.GENEONT.GETMATRIX 
%  GENEONT.GENEONT.SUBSREF GOANNOTREAD NUM2GOID

% References:
% [1]  http://www.geneontology.org/

%   Copyright 2002-2005 The MathWorks, Inc.



relTypes = {'both','is_a','part_of'};
% set defaults 
exclude = false;
depth = 1;
height = 1;
relation = 0;
levels = [];

% check  input arguments
if nargin > 2
    if rem(nargin,2) ~= 0
        error(message('bioinfo:geneont:getrelatives:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'height','depth','exclude','relationtype','levels'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        if(isstruct(pval))
            error(message('bioinfo:geneont:getrelatives:StructParamError'));
        end
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:geneont:getrelatives:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:geneont:getrelatives:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % height
                    if(pval >= 0 && isnumeric(pval))
                        height =  pval;
                    else
                        error(message('bioinfo:geneont:getrelatives:IncorrectHeight', upper( char( okargs( k ) ) )));
                    end
                case 2  % depth
                    if(pval >= 0 && isnumeric(pval))
                        depth =  pval;
                    else
                        error(message('bioinfo:geneont:getrelatives:IncorrectDepth', upper( char( okargs( k ) ) )));
                    end
                case 3  % exclude
                    exclude = bioinfoprivate.opttf(pval);
                    if isempty(exclude)
                        error(message('bioinfo:geneont:getrelatives:ExcludeNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 4 % relation
                    if isnumeric(pval) && iscalar(pval) && pval>=0 && pval<=2
                        relation = round(pval);
                    elseif ischar(pval)
                        relation = strmatch(lower(pval),relTypes)-1; 
                        if isempty(relation) 
                            error(message('bioinfo:geneont:getrelatives:IncorrectRelationTypeString', pval))
                        end
                    else
                        error(message('bioinfo:geneont:getrelatives:IncorrectRelationType'))
                    end     
                case 5  % levels
                    if(pval >= 0 && isnumeric(pval))
                        levels =  pval;
                    else
                        error(message('bioinfo:geneont:getrelatives:IncorrectLevels', upper( char( okargs( k ) ) )));
                    end
            end
        end
    end
end

% remove any IDs that are invalid
null_terms = GO.hashID(terms+1)==0;
if any(null_terms)
    warning(message('bioinfo:geneont:getrelatives:InvalidID', num2str( terms( null_terms ) )));
    terms(null_terms) = [];
end

if isempty(terms)
   GOind = [];
   if nargout>0
       idx = zeros(size(GO.gonumlist));
   end
   return
end

if ~isempty(levels)
    if levels==0 % quick shortcut for just getting the counts
        GOind = unique(terms);
        idx = accumarray(full(GO.hashID(terms+1)),1,[numel(GO.gonumlist) 1]);
        return
    else
        height = levels;
        depth = levels;
    end
end


if nargin == 1
    GOind = unique([GO.getdescendants(terms,'depth',depth,'relation',relation,'exclude',exclude);...
                    GO.getancestors(terms,'height',height,'relation',relation,'exclude',true)]);
else
    [GOdec idxDec] = GO.getdescendants(terms,'depth',depth,'relation',relation,'exclude',exclude);
    [GOanc idxAnc] = GO.getancestors(terms,'height',height,'relation',relation,'exclude',true);
    GOind = unique([GOdec;GOanc]);
    idx = idxDec+idxAnc;
end
