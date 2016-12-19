function plotcngenomeprofile(X,Y, CID, chromIDs, sample, type)
%PLOTCNGENOMEPROFILE plots copy number genome profile.
% 
%   PLOTCNGENOMEPROFILE(GP,CNRATIO,CHROMID, ALLCHROMIDS, SAMPLENAME) plots
%   copy number genome profile.  CP is a vector contains genomic positions
%   of the copy number target. CNRATIO is a vector contains the log2-based
%   intensity ratio for each target. CHROMID specifies the chromosome
%   numbers to be plotted. ALLCHROMIDS specifies the total chromosomes for
%   the whole genome.  SAMPLENAME is the name of the sample to be plot.
% 
%   PLOTCNGENOMEPROFILE(..., TYPE) specifies the type of the plot. The
%   choices are 'W' (default) to plot all chromosomes in one plot, and 'S'
%   to plot each chromosome in a subplot.

%   Copyright 2008 The MathWorks, Inc.


if nargin < 6
    type = 'W';
end

if strcmpi(type, 'S')
    plotprofileSub(X,Y, CID, chromIDs, sample)
else
    plotprofile(X,Y, CID, chromIDs, sample)
end

end

function plotprofile(X,Y, CID, chromIDs, sample)
% Plot the entile genome
n = numel(chromIDs);
chr_endIdx = zeros(1, n);
chr_data_len = zeros(1,n);
for i = 1:n
    tmp = CID == chromIDs(i);
    chr_endIdx(i) = find(tmp, 1, 'last');
    chr_data_len(i) = length(find(tmp));
end

x_lims = [0 chr_endIdx(n)];
y_lims = [fix(min(min(Y), -2)), ceil(max(max(Y), 2))];

% Draw a vertical bar at the end of a chromosome to indicate the border
x_vbar = repmat(chr_endIdx, 2, 1);
y_vbar = repmat(y_lims', 1, n);
offset = double(y_lims(1))-0.2;
% Label the autosome with their chromosome numbers
x_label = chr_endIdx - ceil(chr_data_len/2);
% % y_label = repmat(offset, 1, length(x_label));
y_label = zeros(1, length(x_label))+ offset;
chr_labels = cellstr(num2str(chromIDs'));
chr_labels = strrep(chr_labels, '23', 'X');
% A gray zero line
x_zero = [0; length(Y)];
y_zero = [0;0];
%== Plot
figure; hold on
h_ratio = plot(Y, '.');

line(x_vbar, y_vbar, 'color', [0.8 0.8 0.8]);
line(x_zero, y_zero, 'Color', [1 0.3 0.3], 'Linewidth', 1.5);
text(x_label, y_label, chr_labels,...
    'Fontsize', 8, 'HorizontalAlignment', 'Center');
h_axis = get(h_ratio, 'parent');
set(h_axis, 'xtick', [], 'ygrid', 'on', 'box', 'on',...
            'xlim', x_lims, 'ylim', y_lims)

title(sample, 'Interpreter', 'none')
xlabel({'', 'Chromosome'})
ylabel('Log2 Ratio')
hold off
end % end of function

function plotprofileSub(X,Y, CID, chromIDs, sample)
% X - genomic positions
% Y - log2ratio
% CID - chromosome ids
% sample - Sample name

n = numel(chromIDs);

if n <= 2
	spn = [1 n];
else 
    spn = [ceil(n/4), 4]; 
end

y_lims = [min(min(Y), -1.5), max(max(Y), 1.5)];

figure;
for c =1:n
	chrom = chromIDs(c);
	idx = CID == chrom;
	x= X(idx);
	y =Y(idx);

    x_lims = [0 max(x)];
	subplot(spn(1), spn(2), c)
	hp = plot(x, y, '.','Color', [0.3 0.3 1]);
    line(x_lims, [0 0], 'Color', [1 0.3 0.3], 'Linewidth', 1.2); %Zero line
    
	ha = get(hp, 'Parent');
	set(ha, 'xtick', [],...
			'box', 'on',...
			'xlim', x_lims,...
			'ylim', y_lims);
	title(sprintf('chr %d', chrom), 'Interpreter', 'none');
end
sth = suptitle(sprintf('%s',sample));
set(sth, 'Interpreter', 'none');
end 
