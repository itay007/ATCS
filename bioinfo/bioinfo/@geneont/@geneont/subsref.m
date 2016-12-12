function T = subsref(GO,s)
%SUBSREF Subscripted reference.
%
%   GSO = GO(IDS) returns a Gene Sub-Ontology object GSO containing only
%   the terms referenced by ID. IDS is a vector with non-negative integers.
%   IDS can also be a vector of GENEONT.TERM objects.
%
%   Example:
%
%      % Retrieve the gene subontology which contains specific terms:
%      t = [5575 5622 5623 5737 5840 30529 43226 43228 43229 43232 43234];
%      GSO = GO(t)
%
%   See also GENEONT GENEONT.GENEONT.GETANCESTORS
%   GENEONT.GENEONT.GETDESCENDANTS GENEONT.GENEONT.GETMATRIX 
%   GENEONT.GENEONT.GETRELATIVES GOANNOTREAD NUM2GOID

%   Copyright 2002-2005 The MathWorks, Inc.


% '()' overloads to integer subscripts, everything else left alone
switch s(1).type
    case '()'
        if numel(s(1).subs) ~= 1
            error(message('bioinfo:geneont:subsref:highDimensionalSubscripting'))
        end
        terms = gettermbyid(GO, s(1).subs{1});
        T = copy(GO);
        if numel(terms)==0
            set(T,'Terms',handle(zeros(0,1)))
            set(T,'gonumlist',[])
            set(T,'hashID',sparse([],[],[],10000000,1))
            tm = get(GO,'isarelationmatrix');
            m = length(tm);
            set(T,'isarelationmatrix',sparse([],[],false,m,m))
            set(T,'isarelationmatrixi',sparse([],[],false,m,m))
            set(T,'partofrelationmatrix',sparse([],[],false,m,m))
            set(T,'partofrelationmatrixi',sparse([],[],false,m,m))
        else
            set(T,'Terms',copy(terms))
            if numel(terms)==1
                gonumlist = get(terms,'id');
            else
                gonumlist = cell2mat(get(terms,'id'));
            end
            set(T,'gonumlist',gonumlist);
            set(T,'hashID',sparse(gonumlist+1,1,1:numel(gonumlist),10000000,1));
            tm = get(GO,'isarelationmatrix');
            m = length(tm);
            [i,j]=find(tm(gonumlist+1,gonumlist+1));
            set(T,'isarelationmatrix',sparse(gonumlist(i)+1,gonumlist(j)+1,true,m,m))
            tm = get(GO,'isarelationmatrixi');
            [i,j]=find(tm(gonumlist+1,gonumlist+1));
            set(T,'isarelationmatrixi',sparse(gonumlist(i)+1,gonumlist(j)+1,true,m,m))
            tm = get(GO,'partofrelationmatrix');
            [i,j]=find(tm(gonumlist+1,gonumlist+1));
            set(T,'partofrelationmatrix',sparse(gonumlist(i)+1,gonumlist(j)+1,true,m,m))
            tm = get(GO,'partofrelationmatrixi');
            [i,j]=find(tm(gonumlist+1,gonumlist+1));
            set(T,'partofrelationmatrixi',sparse(gonumlist(i)+1,gonumlist(j)+1,true,m,m))
        end
    otherwise
        T = builtin('subsref', GO, s(1));
end

% if this is not the last level of subsref, do the rest
if length(s)>1
  T = subsref(T, s(2:end));
end

% -------------------------------------------------------------------------
function T = gettermbyid(GO,id)
if isa(id,'geneont.term') % allowed to reference by the terms themselves
    if any(~ismember(id,GO.Terms))
        error(message('bioinfo:geneont:subsref:notValidTerm'));
    else
        T = id;
        return
    end
end
if  ~isnumeric(id)
    error(message('bioinfo:geneont:subsref:IDNotNumeric'));
end
id=id(:);
if any(~real(id+1) | id<0 | rem(id,1))
    error(message('bioinfo:geneont:subsref:NotAPositiveInteger'));
end
if any(id>9999999)
    error(message('bioinfo:geneont:subsref:IDtooLarge'));
end
hID = get(GO,'hashID');
idx = full(hID(id+1));
if ~all(idx)
    warning(message('bioinfo:geneont:subsref:IDnotFound'));
end
GOT = get(GO,'Terms');
T = GOT(nonzeros(idx));
