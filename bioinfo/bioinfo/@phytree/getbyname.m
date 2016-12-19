function sel = getbyname(tr,query,varargin)
%GETBYNAME Selects branches and leaves by name.
%
%   S = GETBYNAME(T,EXPRESSION) returns a logical vector S of size 
%   [NUMNODES x 1] indicating the node names of the phylogenetic tree T
%   that match the regular expression EXPRESSION regardless of case. 
%
%   Symbols than can be used in a matching regular expression are explained
%   in help REGEXP.
%
%   When EXPRESSION is a cell array of strings, GETBYNAME returns a matrix
%   where every column corresponds to every query in EXPRESSION.
%
%   S = GETBYNAME(T,STRING,'EXACT',true) finds only exact matches, ignoring
%   case. When STRING is a cell array of strings, GETBYNAME returns a
%   vector with indices for each node to the corresponding match in STRING,
%   using zero for nodes without a match.
%
%   Example:
%
%      % Load a phylogenetic tree created from a protein family:
%      tr = phytreeread('pf00002.tree');
%       
%      % Select all the 'mouse' and 'human' proteins:
%      sel = getbyname(tr,{'mouse','human'});
%      view(tr,any(sel,2));
%       
%   See also PHYTREE, PHYTREE/PRUNE, PHYTREE/SELECT, PHYTREE/GET.

% Copyright 2003-2009 The MathWorks, Inc.


if numel(tr)~=1
     error(message('bioinfo:phytree:getbyname:NoMultielementArrays'));
end

doExactMatch = false;

if  nargin > 2
    okargs = {'exact',''};
    for j=1:2:nargin-2
        pname = varargin{j};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:phytree:getbyname:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:phytree:getbyname:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1
                   if nargin == 3
                       doExactMatch = true;
                   else
                       doExactMatch = bioinfoprivate.opttf(varargin{j+1});
                       if isempty(doExactMatch)
                           error(message('bioinfo:phytree:getbyname:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                       end
                   end
            end %switch
        end %if
    end %for
end %if


numLabels = numel(tr.names);
if iscell(query)
    if doExactMatch
        sel = zeros(numLabels,1);
    else
        sel = false(numLabels,numel(query));
    end
    for ind = 1:numel(query)
        if doExactMatch
            sel(strcmpi(query{ind},tr.names)) = ind;
        else
            try
                regexpiOutput = regexpi(tr(:).names,query{ind});
            catch theException
                error(message('bioinfo:phytree:getbyname:IncorrectRegularExpression', theException.message));
            end
            sel(:,ind) = ~cellfun('isempty',regexpiOutput);
        end
    end
else % must be a single string of chars
    if doExactMatch
        sel = strcmpi(query,tr.names);
    else
        try
            regexpiOutput = regexpi(tr(:).names,query);
        catch theException
            error(message('bioinfo:phytree:getbyname:IncorrectRegularExpression', theException.message));
        end
        sel = ~cellfun('isempty',regexpiOutput);
    end
end
