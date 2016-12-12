function mapTaxoFile(taxonomyFilenameIn, taxonomyFilenameOut, blockSize)
%MAPTAXOFILE Helper function for METAGENOMICDEMO

%   Copyright 2007-2012 The MathWorks, Inc.


if exist(taxonomyFilenameOut,'file')
    taxonomyFilenameInInfo = dir(which(taxonomyFilenameIn));
    taxonomyFilenameOutInfo = dir(which(taxonomyFilenameOut));
    if taxonomyFilenameOutInfo.datenum > taxonomyFilenameInInfo.datenum
       warning('bioinfo:mapTaxoFile:FileExists','Memory mapped file exists.')
       return
    end
end

fid1 = fopen(which(taxonomyFilenameIn),'rt'); % from NCBI TAXONOMY FTP site
if fid1<0
    error('bioinfo:mapTaxoFile:invalidFile','Cannot open input file.')
end
fid2 = fopen(taxonomyFilenameOut, 'w');% binary file used for mapping

%===create a map between gi numbers and taxids
curr = 1; % current gi to consider
while(~feof(fid1))
    
    data = textscan(fid1, '%d %d', blockSize);
    gi = data{1};
    taxa = data{2};
    gap = gi(1) - curr; 
    
    %=== missing gi numbers between blocks are assigned a taxid = -1
    if gap
        D = -1 * ones(gap, 1); 
        fwrite(fid2, D, 'int32');
    end
        
    %=== populate array D such that D(gi) = taxid of gi
    curr = gi(end) + 1;   % current gi position in the final list
    offset = min(gi) - 1; % starting gi in the current block
    N = max(gi) - offset; % number of gi's  to consider
    D = -1 * ones(N,1);   
    D(gi - offset) = taxa;
       
    %=== write array D into binary file
    fwrite(fid2, D, 'int32');
end
fclose all;
    
