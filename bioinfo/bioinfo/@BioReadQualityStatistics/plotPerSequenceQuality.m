function hax = plotPerSequenceQuality(qc,hax)
%PLOTPERSEQUENCEQUALITY Bar graph showing read count versus average quality
%   PLOTPERSEQUENCEQUALITY(qs) displays a bar  graph showing the percentage
%   of short reads versus GC content.
%
%   h = PLOTPERSEQUENCEQUALITY(qs) returns the handle to the axes object
%   containing the plot.

if nargin == 1
    hax = gca;
end


nqs = numel(qc);

mn = min(arrayfun(@(q)(q.MinEncodingPhred) , qc ));
mx = max(arrayfun(@(q)(q.MaxEncodingPhred) , qc ));
x = mn:mx;
y = zeros(length(x),nqs);
lim = [0 40];

for i = 1:numel(qc)
    col = qc(i).PerSeqAverageQualityDist/sum(qc(i).PerSeqAverageQualityDist)*100;
    ind  = ( qc(i).MinEncodingPhred:(qc(i).MaxEncodingPhred-1) )-mn+1;
    y(ind,i) = col;
    lim = [min(lim(1),qc(i).MinScore) max(lim(2),qc(i).MaxScore)]; 
end

bar(x,y);
xlim(lim);
set(hax,'xtick',lim(1):2:lim(2));

labels = cellstr(get(hax,'XTickLabel'));
labels(2:2:end) = {''};
set(gca,'XTickLabel',labels);

xlabel(hax,'Average Score');
ylabel(hax,'% Read-Count');
title(hax,'Quality Distribution');




