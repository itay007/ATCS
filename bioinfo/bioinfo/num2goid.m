function IDS = num2goid(X)
%NUM2GOID Converts numbers to Gene Ontology IDs.
%
%   GOID = NUM2GOID(X) converts the numbers in X to strings with Gene
%   Ontology IDs. IDs are a 7 digit number preceded by the prefix 'GO:'. 
%
%   Example:
%
%      % Get the Gene Ontology IDs of the following numbers:
%      t = [5575 5622 5623 5737 5840 30529 43226 43228 43229 43232 43234];
%      ids = num2goid(t)
%
%   See also GENEONT GENEONT.GENEONT.GETANCESTORS
%   GENEONT.GENEONT.GETDESCENDANTS GENEONT.GENEONT.GETMATRIX 
%   GENEONT.GENEONT.GETRELATIVES GENEONT.GENEONT.SUBSREF GOANNOTREAD

%   Copyright 2002-2007 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);

if  ~isnumeric(X)
    error(message('bioinfo:num2goid:IDNotNumeric'));
end
if any(~isreal(X(:))) || any(X(:)<0) || any(rem(X(:),1))
    error(message('bioinfo:num2goid:NotAPositiveInteger'));
end

if any(X(:)>9999999)
    error(message('bioinfo:num2goid:IDtooLarge'));
end

IDS = cell(size(X));
for i = 1:numel(X)
    IDS{i} = sprintf('GO:%7d',X(i));
end
IDS = regexprep(IDS,' ','0');

