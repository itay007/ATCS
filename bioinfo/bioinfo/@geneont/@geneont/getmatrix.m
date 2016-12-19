function [cm t rels] = getmatrix(GO,varargin)
%GETMATRIX converts a Gene Ontology Object into a relationship matrix
%
%   [MATRIX, ID, RELATIONSHIP] = GETMATRIX(GO) converts the Gene Ontology
%   object GO into a matrix of relationship values. ID is a list of Gene
%   Ontology ids that correspond to the rows and columns of MATRIX. The
%   values in the matrix are indices of the relationship types in
%   RELATIONSHIP, usually 1 for 'is_a' and 2 for 'part_of'.
%
%   Example:
%
%      [MATRIX ID REL] = getmatrix(GO); 
%
%   See also GENEONT GENEONT.GENEONT.GETANCESTORS
%   GENEONT.GENEONT.GETDESCENDANTS GENEONT.GENEONT.GETRELATIVES
%   GENEONT.GENEONT.SUBSREF GOANNOTREAD NUM2GOID

%   References:
% [1]  http://www.geneontology.org/

%   Copyright 2002-2005 The MathWorks, Inc.



 %%% check arguments
if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:geneont:getmatrix:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'no_input_arguments'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        if(isstruct(pval))
             error(message('bioinfo:geneont:getmatrix:StructParamError'));
        end
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:geneont:getmatrix:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:geneont:getmatrix:AmbiguousParameterName', pname));
        else
            switch(k)
            end
        end
    end
end

cm = GO.relationmatrix(GO.gonumlist+1,GO.gonumlist+1);
t = GO.gonumlist;
rels = GO.relationtype;
