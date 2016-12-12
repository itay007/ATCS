%% Working with SFF Files from the 454 Genome Sequencer FLX System
% This example shows how to read SFF (Standard Flowgram Format) files,
% which store the sequencing trace data produced by the 454 GS FLX System.

% Copyright 2009-2012 The MathWorks, Inc.


%%
% *Note*: This example assumes the SFF file SRR000001.sff is in your
% current folder. You can download the SRA file for accession |SRR000001|
% from this
% <ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX000/SRX000007/SRR000001/
% NCBI FTP site> and convert it to SFF formatted file |SRR000001.sff| using
% the <http://www.ncbi.nlm.nih.gov/books/NBK47540/ NCBI SRA Toolkit>.
% Alternatively, you can change the name of the variable |filename| to run
% this example on a different file.

filename = 'SRR000001.sff';

%% Displaying Summary Information About the SFF File
% Every SFF file includes a header section with general information
% regarding the file content. To display this information, you can use the
% |sffinfo| function as shown below:

info = sffinfo(filename)

%%
% The summary includes the name of the file, the size of the file in bytes,
% the version of the SFF format, the number of reads included in the file,
% the key sequence, etc. 

%% Selecting SFF File Features
% SFF files are generally large files (in the order of hundreds of
% megabytes), and reading their entire content into memory can require
% ample computational resources. One possible solution to this problem is
% to read the file in blocks, and process relevant information in
% chunks. Another possible solution is to limit the amount of information 
% extracted from the file, by selecting only those features needed for the
% analysis at hand.
%
% The |sffread| function lets you employ both the
% strategies described above through the |feature| and |block| options. For
% example, you can read the file in blocks of size M entries:

N = info.NumberOfReads; % total number of reads in the file
M = 100; % number of reads to consider

sffread(filename, 'block', [1 M])
sffread(filename, 'block', [201 200+M])

%%
% Each read data section included in the SFF file consists of the reads'
% universal accession numbers (|h|), sequence information (|s|), quality
% scores of basecalls (|q|), clipping positions (|c|), flowgram values
% (|f|), and flowgram indices (|i|). By using the |feature| option you can
% select the combination of features that are relevant to your analysis and
% ignore the data that is not.  
%
% You can use the |sffread| function to read an SFF file into an array of
% structures. Each structure contains the default combination of features,
% including the reads's universal accession numbers (in the |Header|
% field), the nucleotide bases (in the |Sequence| field) and the quality
% scores (in the |Quality| field).

sffread(filename, 'block', [1 M])

%%
% Other feature combinations are shown below:

%=== extract only sequences
s = sffread(filename, 'block', [1 100], 'feature', 's')

%=== extract only accession numbers and sequences
hs = sffread(filename, 'block', [1 100], 'feature', 'hs')

%=== extract only sequences and quality scores
sq = sffread(filename, 'block', [1 100], 'feature', 'sq')

%=== extract only sequences and flowgram values
sf = sffread(filename, 'block', [1 100], 'feature', 'sf')

%=== extract only sequences and clipping info
sc = sffread(filename, 'block', [1 100], 'feature', 'sc')

%=== extract all info available
all = sffread(filename, 'block', [1 100], 'feature', 'hscqfi')

whos s hs sq sf sc all

clear s hs sq sf sc

%%
% *Note*: the nucleotide sequence is represented by 8-bit unsigned integers
% to minimize the memory space needed to store the sequence information
% and to comply with the original representation in the SFF file. To
% convert the sequence representation from uint8 to char, consider the
% following transformation:

s1 = all(1).Sequence(1:10)

s2 = char(double(s1))

%% Clipping the Sequences
% The four integers stored in the |Clipping| field of the output structure
% represent the clipping points for the given read. The four integers represent:
%
% * int1: position of the first base at the beginning of the sequence after
% the clipping point for quality;
% * int2: position of the last base at the end of the sequence before the 
% clipping point for quality;
% * int3: position of the first base at the beginning of the sequence after
% the clipping point for an adapter;
% * int4: position of the last base at the end of the sequence before the 
% clipping point for an adapter.
%
% Let us consider one specific read and apply the clipping to the termini
% as specified by the clipping points:

targetName = 'EM7LVYS01DSWC4'; 

targetIndex = find(strcmp({all.Header}, targetName));
t = all(targetIndex);

%=== transform sequence into character representation
seq = char(double(t.Sequence));

%=== first base after clipping, last base before clipping
clip = t.Clipping(1:2)

%=== sequence before clipping
seqdisp(seq) 

%=== sequence after clipping
seqdisp(seq(clip(1):clip(2)))

%% Plotting Quality Scores
% Each base in a given read has a quality score, Q, associated with it. Q =
% -10 * log10(p), where p is the probability error. The quality scores
% range from 0 (worst) to 40 (best). You can visualize the quality trend in
% a given sequence read as follows:

q = t.Quality;

%=== determine nucleotide positions by type
xa = find(seq == 'A');
xc = find(seq == 'C');
xg = find(seq == 'G');
xt = find(seq == 'T');

%=== plot in different colors
stem(xa, q(xa), 'linewidth', 2, 'MarkerSize', 0, 'color', 'g'); hold on;
stem(xc, q(xc), 'linewidth', 2, 'MarkerSize', 0, 'color', 'b');
stem(xg, q(xg), 'linewidth', 2, 'MarkerSize', 0, 'color', 'k');
stem(xt, q(xt), 'linewidth', 2, 'MarkerSize', 0, 'color', 'r'); hold off;
xlabel('Sequence position');
ylabel('Quality score');
title(['Base quality scores for sequence ' targetName]);

%% Representing Quality and Sequence
% You can display the nucleotides in the sequence with the corresponding
% quality scores as shown below:

sq = '';
for i = 1:floor(length(q)/10) % rows
	for j = 1:10
		x = (i-1)*10 + j; % position within sequence
		sq = [sq sprintf('%s:%-2d\t', seq(x), q(x))];
	end
	sq = [sq sprintf('\n')]; 
end

for j = 1:mod(length(q),10)
	sq = [sq sprintf('%s:%-2d\t', seq(x+j), q(x+j))];
end
sq = [sq sprintf('\n')]


%%
% You can append the sequence position before each base, as follows:

sq = '';
for i = 1:floor(length(q)/10) % rows
	for j = 1:10
		x = (i-1)*10 + j; % position within sequence
		sq = [sq sprintf('%4d_%s:%-2d\t', x,seq(x), q(x))];
	end
	sq = [sq sprintf('\n')]; 
end

for j = 1:mod(length(q),10)
	sq = [sq sprintf('%4d_%s:%-2d\t', x+j, seq(x+j), q(x+j))];
end
sq = [sq sprintf('\n')]

%% Performing Metrics Analysis on Reads' Length and Quality
% To have a sense of the quality of the reads produced by a given run, you
% can perform standard statistical analysis on the reads' length and
% quality scores.

%=== consider a block of size M
M = 10000;
all = sffread(filename, 'block', [1 M], 'feature', 'sqc');

%=== metrics on the length
len = cellfun(@length, {all(:).Sequence});
meanLength = mean(len)
medianLength = median(len)
stdLength = std(len)

figure(); hist(len); 
xlabel('Sequence length');
ylabel('Number of reads');
title('Length distribution');

%=== metrics on the quality
meanQuality = cellfun(@mean, {all.Quality});
medianQuality = cellfun(@(x) median(double(x)), {all.Quality});

figure(); hist(medianQuality);
xlabel('Median Quality Score');
ylabel('Number of reads');
title('Quality score distribution')

%% Plotting Flowgram Intensity
% The sequencing method behind the reads produced with the 454 GS FLX
% System evaluates the intensity of the luminiscent reaction by luciferase
% every time a nucleotide is incorporated in the elongating chain during
% the sequencing by synthesis process. This signal intensity is reported in
% the |FlowgramValue| field of the output structure and can be plotted as
% shown below:

%=== determine the flow index where each base is called
fi = double(t.FlowgramIndex);
fv = double(t.FlowgramValue);
x = cumsum(fi);

N = info.NumberOfFlowsPerRead;
f = info.FlowChars(1:4);
z = 1:4:N;

%=== determine flow index by nucleotide type
za = strfind(f,'A')-1;
zc = strfind(f,'C')-1;
zg = strfind(f,'G')-1;
zt = strfind(f,'T')-1;

%=== determine flow index with base called by nt type
xa = intersect(z+za,x);
xc = intersect(z+zc, x);
xg = intersect(z+zg, x);
xt = intersect(z+zt, x);

%=== plot flowgram intensity with nt color scheme
figure(); hold on;
h = stem(fv, 'linewidth', 2, 'MarkerSize', 0, 'color', [.8 .8 .8]); 
ha = stem(xa, fv(xa), 'linewidth', 2, 'MarkerSize', 0, 'color', 'g');
hc = stem(xc, fv(xc), 'linewidth', 2, 'MarkerSize', 0, 'color', 'b');
hg = stem(xg, fv(xg), 'linewidth', 2, 'MarkerSize', 0, 'color', 'k');
ht = stem(xt, fv(xt), 'linewidth', 2, 'MarkerSize', 0, 'color', 'r');
ylabel('Flow Intensity');
xlabel('Flow Sequence');
hold off;

%% Displaying Homopolymers in the Flowgram Plot
% Highest peaks in the plot of the flowgram intensity correspond to the
% presence of homopolymers, or uninterrupted stretches of a single
% nucleotide, in the sequence. You can observe the presence of homopolymers
% in the plot below: 

values = accumarray(x',1);
labels = regexp(info.FlowChars(1:length(values)), '[ACGTN]', 'match')';

figure(); hold on; 
stem(xa, values(xa), 'linewidth', 2, 'markersize', 0, 'color', 'g');
stem(xc, values(xc), 'linewidth', 2, 'markersize', 0, 'color', 'b');
stem(xg, values(xg), 'linewidth', 2, 'markersize', 0, 'color', 'k');
stem(xt, values(xt), 'linewidth', 2, 'markersize', 0, 'color', 'r');
yl = ylim; ylim([0 yl(2)+1]);

hpindex = values > 1;
hpvalues = values(hpindex);
hplabels = labels(hpindex);

text(find(hpindex), hpvalues+.1, hplabels, 'fontsize', 6, ...
	'horizontalAlignment', 'center')
xlabel('Sequence position');
ylabel('Number of bases');
hold off;

%%
% The intensity of the luciferase activity is proportional to the number of
% identical nucleotides (forming homopolymers) incorporated in the
% synthesized chain. When this number is too high, the intensity can reach
% saturation, thus making it difficult to quantify the exact number of
% contiguous identical nucleotides that have been incorporated. This leads
% to miscalled bases in the presence of homopolymers. You can survey the number
% of bases in a homopolymer stretch and the various types of homopolymers as
% follows:
 
%=== consider numeric representation of the sequence
seqNum = [0 double(t.Sequence)];
n = length(seqNum);

%=== identify bases in homopolymers (hp)
x = sum(diff(seqNum) == 0); % bases in hp (less starting base)
numHomopolymers = sum(diff((diff(seqNum) == 0)) == 1) % initial bases of hp

homopolBases = x + numHomopolymers

%=== consider homopolymers up to K = 9
K = 9; 
C = zeros(4,K);
pos = cell(4,K);

for j = 1:K
	k = j+1; % size of homopolymer
	
	%=== position of homopolymers
	pos{1,j} = find(values(xa) == k);
	pos{2,j} = find(values(xc) == k);
	pos{3,j} = find(values(xg) == k);
	pos{4,j} = find(values(xt) == k);
	
	%=== counts of homopolymers
	C(1,j) = sum(values(xa) == k);
	C(2,j) = sum(values(xc) == k);
	C(3,j) = sum(values(xg) == k);
	C(4,j) = sum(values(xt) == k);
end

figure(); bar(C');
set(gca, 'Xtick', [1:K-1], 'XtickLabel', [2:K])
xlabel('Homopolymer length');
ylabel('Number of homopolymer');
title(['Homopolymers in sequence ' targetName]);
legend('A', 'C', 'G', 'T')

%%
% The element in row |i| and column |j| of the matrix |C| holds the number
% of homopolymers of size k (k = 2,...,9, k = j+1) for nucleotide type i,
% where i = 1,2,3,4 correspond to A,C,G,T respectively. The element in row
% |i| and column |j| of the cell array |pos| holds the sequence flow
% position, where a homopolymer of size k (k = 2,...,9, k = j+1) for
% nucleotide i (i = 1,2,3,4) is found.

%% Saving Information into FASTQ Files
% Once you have parsed and processed the information contained in the SFF
% file, it can be convenient to save some (or all) sequences and
% quality scores into flat text files for future consideration and
% analysis. For this purpose, you can use the |fastqwrite| function.


%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20SFFDEMO%20in%20Bioinformatics%20Toolbox%204.2 Provide feedback for this example.>*

displayEndOfDemoMessage(mfilename)
