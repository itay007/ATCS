function hax = plotTotalGC(qc,hax)
%PLOTTOTALGC Bar graph displaying distribution of nucleotides
%   PLOTTOTALGC(qs) displays a bar graph showing the  nucleotide counts in
%   the short reads described by the BioReadQualityStatistics object qs. 
%
%   h = PLOTTOTALGC(qs) returns the handle to the axes object containing the
%   plot.



if nargin == 1
    hax = gca;
end


nqs = numel(qc);
totals = zeros(5,nqs);
for i =1:nqs
    col = sum(qc(i).PerPosBaseDist,2);
    col = col/sum(col)*100;
    totals(:,i) = col;
end
bar(hax,1:5,totals);
axis(hax,[0 6 0 1.2*max(totals(:))]);
set(hax,'XTick',1:5);
set(hax,'XTickLabel',{'A' 'C' 'G' 'T' 'Other'})
set(hax,'TickLength',[0 0]);
xlabel(hax,'Nucleotide');
ylabel(hax,'% Nucleotides');
title(hax,'Nucleotide Distribution');