function hax = plotPerSequenceGC(qc,hax)
%PLOTPERSEQUENCEGC Line graph showing read count versus perctange GC content
%   PLOTPERSEQUENCEGC(qs) displays a line graph showing the percentage of
%   short reads versus GC content.
%
%   h = PLOTPERSEQUENCEGC(qs) returns the handle to the axes object
%   containing the plot.



if nargin == 1
    hax = gca;
end
pos = 2.5:5:97.5;

nqs = numel(qc);
totals = zeros(length(pos),nqs);
for i = 1:nqs
    totals(:,i) = qc(i).PerSeqGCDist/qc(i).NumberOfReads;
end 
totals = totals*100;
bar(hax, pos,totals);
xlabel(hax, '% GC-Content');
ylabel(hax, '% Read-Count');
title(hax, 'GC Distribution');
set(hax,'Xtick',pos);
labels = cellstr(num2str((0:5:100)'));
labels(2:2:end) = {''};
set(hax,'XTickLabel',labels);