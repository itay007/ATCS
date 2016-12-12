function taxidRank = findTaxoRank (taxid, rankObj, parentObj, rTarget)
%FINDTAXORANK Helper function for METAGENOMICDEMO

%   Copyright 2007-2012 The MathWorks, Inc.

N = numel(taxid);
unclassified = true(N,1);
taxidRank = zeros(N,1);

queryTaxid = taxid;
index = 1:N;

% find the taxonomic assignment according to specified rank rTarget, by
% walking up the hierarchy given by parentObj

while(any(unclassified))
    unclassifiedRank = rankObj.Data(queryTaxid); % get rank of unclassified entries

    out = unclassifiedRank < rTarget;
    unclassified(out) = 0;
    taxidRank(index(out)) = 1;  % assign it to the root
    
    done = unclassifiedRank == rTarget;
    unclassified(done) = 0;
    taxidRank(index(done)) = queryTaxid(done); 

    queryTaxid = parentObj.Data(queryTaxid(unclassified));
    index = index(unclassified);
 
    if all(queryTaxid == 1)
        taxidRank(index) = 1; 
        unclassified = 0;
    else
          unclassified = true(numel(queryTaxid),1);
    end
end
    
    

