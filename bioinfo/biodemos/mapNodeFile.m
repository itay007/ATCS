function mapNodeFile(nodeFilename,  parentFilenameOut, rankFilenameOut, blockSize)
%MAPNODEFILE Helper function for METAGENOMICDEMO

%   Copyright 2007-2012 The MathWorks, Inc.


if exist(parentFilenameOut,'file') && exist(rankFilenameOut,'file')
    nodeFilenameInfo = dir(which(nodeFilename));
    parentFilenameOutInfo = dir(which(parentFilenameOut));
    rankFilenameOutInfo = dir(which(rankFilenameOut));
    if parentFilenameOutInfo.datenum > nodeFilenameInfo.datenum && ...
       rankFilenameOutInfo.datenum > nodeFilenameInfo.datenum         
       warning('bioinfo:mapNodeFile:FileExists','Memory mapped file exists.')
       return
    end
end

rank = {'superkingdom', 'kingdom', 'subkingdom', ...
    'superphylum', 'phylum', 'subphylum', ...
    'superclass', 'class', 'subclass', ...
    'superorder', 'order', 'suborder', ...
    'superfamily', 'family', 'subfamily',... 
    'supergenus', 'genus', 'subgenus',...
    'superspecies', 'species', 'subspecies', 'norank'};
code = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};


fid1 = fopen(which(nodeFilename),'rt');    % from NCBI TAXONOMY FTP site
if fid1<0
    error('bioinfo:mapNodeFile:invalidFile','Cannot open input file.')
end
fid2 = fopen(parentFilenameOut, 'w');% binary file used for mapping
fid3 = fopen(rankFilenameOut, 'w');% binary file used for mapping
curr = 1; % current taxid to consider

%=== create a map between taxid and rank
%=== create a map between taxid and parent taxid
while(~feof(fid1))
         
    data = textscan(fid1, '%d%d%s%s%d%d%d%d%d%d%d%d%s', blockSize, 'delimiter', '|'); 
    node = data{1};
    parent = data{2};
    classif = strtrim(data{3});
    classifInt = ones(numel(classif),1) * 22; % no rank
    for i = 1:numel(rank)
        K = strmatch(rank{i}, classif); % entries that are classified with the specific rank
        classifInt(K) = code{i}; % assign rank code
    end
    gap = node(1) - curr; 
    
    %=== missing node numbers between blocks are assigned a rank = -1
    if gap
        P = -1 * ones(gap, 1); % parent
        R = -1 * ones(gap, 1); % rank
        fwrite(fid2, P, 'int32');
        fwrite(fid3, R, 'int32');
    end
        
    %=== populate array P such that P(node) = parent of node
    %=== populate array R such that R(node) = rank of node (integer)
    
    curr = node(end) + 1;   % current gi position in the final list
    offset = min(node) - 1; % starting gi in the current block
    N = max(node) - offset; % number of gi's  to consider
    
    %=== if not found, the parent is the root and the rank is 1
    P = ones(N,1);   
    P(node - offset) = parent;
    R = ones(N,1);
    R(node - offset) = classifInt;
          
    %=== write array D into binary file
    fwrite(fid2, P, 'int32');
    fwrite(fid3, R, 'int32');
end
fclose all;
    

