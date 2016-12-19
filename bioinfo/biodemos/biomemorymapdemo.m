%% Working with Whole Genome Data
% This example shows how to create a memory mapped file for sequence data
% and work with it without loading all the genomic sequence into memory. Whole
% genomes are available for human, mouse, rat, fugu, and several other
% model organisms. For many of these organisms one chromosome can be
% several hundred million base pairs long. Working with such large data
% sets can be challenging as you may run into limitations of the hardware
% and software that you are using. This example shows one way to work
% around these limitations in MATLAB(R).

%   Copyright 2005-2012 The MathWorks, Inc.

%% Large Data Set Handling Issues
% Solving technical computing problems that require processing and
% analyzing large amounts of data puts a high demand on your computer
% system. Large data sets take up significant memory during processing and
% can require many operations to compute a solution. It can also take a
% long time to access information from large data files.
% 
% Computer systems, however, have limited memory and finite CPU speed.
% Available resources vary by processor and operating system, the latter of
% which also consumes resources. For example:
% 
% 32-bit processors and operating systems can address up to 2^32 =
% 4,294,967,296 = 4 GB of memory (also known as virtual address space).
% Windows(R) XP and Windows(R) 2000 allocate only 2 GB of this virtual memory to
% each process (such as MATLAB). On UNIX(R), the virtual memory allocated to a
% process is system-configurable and is typically around 3 GB. The
% application carrying out the calculation, such as MATLAB, can require
% storage in addition to the user task. The main problem when handling
% large amounts of data is that the memory requirements of the program can
% exceed that available on the platform. For example, MATLAB generates an
% "out of memory" error when data requirements exceed approximately 1.7 GB
% on Windows XP.
%
% More details on memory management and large data sets can be found in this
% <http://www.mathworks.com/company/newsletters/digest/nov04/newfeatures.html 
% MATLAB Digest article> or in the
% <http://www.mathworks.com/support/tech-notes/1100/1106.html Memory  
% Management Guide> in our support site. 
% 
% On a typical 32-bit machine, the maximum size of a single data set that
% you can work with in MATLAB is a few hundred MB, or about the size of a
% large chromosome. Memory mapping of files allows MATLAB to work around
% this limitation and enables you to work with very large data sets in an
% intuitive way.

%% Whole Genome Data Sets
% The latest whole genome data sets can be downloaded from the
% <http://www.ensembl.org/info/data/ftp/index.html Ensembl Website>. The
% data are provided in several formats. These are updated regularly as new
% sequence information becomes available. This example will use human DNA
% data stored in FASTA format. Chromosome 1 is (in the GRCh37.56 Release of
% September 2009) a 65.6 MB compressed file. After uncompressing the file it is
% about 250MB. MATLAB uses 2 bytes per character, so if you read the file
% into MATLAB, it will require about 500MB of memory. 
%
% This example assumes that you have already downloaded and uncompressed
% the FASTA file into your local directory. Change the name of the variable
% |FASTAfilename| if appropriate.
                
FASTAfilename = 'Homo_sapiens.GRCh37.56.dna.chromosome.1.fa';
fileInfo = dir(which(FASTAfilename))

%% Memory Mapped Files
% Memory mapping allows MATLAB to access data in a file as though it is in
% memory. You can use standard MATLAB indexing operations to access data.
% See the documentation for |memmapfile| for more details.

%%
% You could just map the FASTA file and access the data directly from
% there. However the FASTA format file includes new line characters. The
% |memmapfile| function treats these characters in the same way as all
% other characters. Removing these before memory mapping the file will make
% indexing operations simpler. Also, memory mapping does not work directly
% with character data so you will have to treat the data as 8-bit integers
% (uint8 class). The function |nt2int| in the Bioinformatics Toolbox(TM) can be
% used to convert character information into integer values. |int2nt| is
% used to convert back to characters.

%%
% First open the FASTA file and extract the header.

fidIn = fopen(FASTAfilename,'r');
header = fgetl(fidIn)

%%
% Open the file to be memory mapped.
[fullPath, filename, extension] = fileparts(FASTAfilename);
mmFilename = [filename '.mm']
fidOut = fopen(mmFilename,'w');

%%
% Read the FASTA file in blocks of 1MB, remove new line characters, convert
% to uint8, and write to the MM file. 

newLine = sprintf('\n');
blockSize = 2^20;
while ~feof(fidIn)
    % Read in the data
    charData = fread(fidIn,blockSize,'*char')';
    % Remove new lines
    charData = strrep(charData,newLine,'');
    % Convert to integers
    intData = nt2int(charData);
    % Write to the new file
    fwrite(fidOut,intData,'uint8');
end

%% 
% Close the files. 
fclose(fidIn);
fclose(fidOut);

%%
% The new file is about the same size as the old file but does not contain
% new lines or the header information.
mmfileInfo = dir(mmFilename)

%% Accessing the Data in the Memory Mapped File
% The command |memmapfile| constructs a memmapfile object that maps the new
% file to memory. In order to do this, it needs to know the format of the
% file. The format of this file is simple, though much more complicated
% formats can be mapped.

chr1 = memmapfile(mmFilename, 'format', 'uint8')

%% The MEMMAPFILE Object
% The memmapfile object has various properties. |Filename| stores the full
% path to the file. |Writable| indicates whether or not the data can be
% modified. Note that if you do modify the data, this will also modify the
% original file. |Offset| allows you to specify the space used by any
% header information. |Format| indicates the data format. |Repeat| is used to
% specify how many blocks (as defined by |Format|) to map. This can be useful
% for limiting how much memory is used to create the memory map. These
% properties can be accessed in the same way as other MATLAB data. For more
% details see type |help memmapfile| or |doc memmapfile|.

chr1.Data(1:10)

%%
% You can access any region of the data using indexing operations.

chr1.Data(10000000:10000010)'

%%
% Remember that the nucleotide information was converted to integers. You
% can use |int2nt| to get the sequence information back.

int2nt(chr1.Data(10000000:10000010)')

%%
% Or use |seqdisp| to display the sequence.

seqdisp(chr1.Data(10000000:10001000)')

%% Analysis of the Whole Chromosome
% Now that you can easily access the whole chromosome, you can analyze the
% data. This example shows one way to look at the GC content along the chromosome.

%%
% You extract blocks of 500000bp and calculate the GC content. 

%%
% Calculate how many blocks to use.

numNT = numel(chr1.Data);
blockSize = 500000;
numBlocks = floor(numNT/blockSize);

%%
% One way to help MATLAB performance when working with large data sets is
% to "preallocate" space for data. This allows MATLAB to allocate enough
% space for all of the data rather than having to grow the array in small
% chunks. This will speed things up and also protect you from problems of
% the data getting too large to store. For more details on pre-allocating
% arrays, see:
% http://www.mathworks.com/support/solutions/data/1-18150.html?solution=1-18150   

%%
% An easy way to preallocate an array is to use the |zeros| function.

ratio = zeros(numBlocks+1,1);

%%
% Loop over the data looking for C or G and then divide this number by the
% total number of A, T, C, and G. This will take about a minute to run.

A = nt2int('A'); C = nt2int('C'); G = nt2int('G'); T = nt2int('T'); 

for count = 1:numBlocks
    % calculate the indices for the block
    start = 1 + blockSize*(count-1);
    stop = blockSize*count;
    % extract the block
    block = chr1.Data(start:stop);
    % find the GC and AT content
    gc = (sum(block == G | block == C));
    at = (sum(block == A | block == T));
    % calculate the ratio of GC to the total known nucleotides
    ratio(count) = gc/(gc+at);
end

%%
% The final block is smaller so treat this as a special case.
block = chr1.Data(stop+1:end);
gc = (sum(block == G | block == C));
at = (sum(block == A | block == T));
ratio(end) = gc/(gc+at);

%% Plot of the GC Content for the Homo Sapiens Chromosome 1

xAxis = [1:blockSize:numBlocks*blockSize, numNT];
plot(xAxis,ratio)
xlabel('Base pairs');
ylabel('Relative GC content');
title('Relative GC content of Homo Sapiens Chromosome 1')

%%
% The region in the center of the plot around 140Mbp is a large region of
% Ns.

seqdisp(chr1.Data(140000000:140001000))

%% Finding Regions of High GC Content 
% You can use |find| to identify regions of high GC content.

indices = find(ratio>0.5);
ranges = [(1 + blockSize*(indices-1)), blockSize*indices];
fprintf('Region %d:%d has GC content %f\n',[ranges ,ratio(indices)]')

%% 
% If you want to remove the temporary file, you must first clear the
% |memmapfile| object.

clear chr1
delete(mmFilename)

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Enhancement%20for%20BIOMEMORYMAPDEMO%20in%20Bioinformatics%20Toolbox%204.2 Request enhancement for this example.>*

displayEndOfDemoMessage(mfilename)
