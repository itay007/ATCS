function hax = plotPerPositionCountByQuality(qcs,hax)
%PLOTPERPOSITIONCOUNTBYQUALITY Line plot displaying the fraction of all
%reads with a Phred scores in various ranges at each position.
% 
%   PLOTPERPOSITIONCOUNTBYQUALITY(qs) displays a  plot showing fraction of
%   all reads that have a Phred score in one of four ranges at each
%   position.
%
%   Note that the scores at a given position need not sum to 100% as some
%   reads may not have a Nucleotide at all positions.
%
%   h = PLOTPERPOSITIONGC(qs) returns the handle to the axes object
%   containing the plot.

if nargin == 1
    hax = gca;
end

mx = 10;
hold(hax,'all');
for i = 1:length(qcs)
    qc = qcs(i);
   ind = qc.MinEncodingPhred:qc.MaxEncodingPhred;
    reads = [sum(qc.PerPosQualities(ind>=0 & ind <=10,:))
         sum(qc.PerPosQualities(ind>=11 & ind <=20,:))
         sum(qc.PerPosQualities(ind>=21 & ind <=30,:))
         sum(qc.PerPosQualities(ind>=31 & ind <=40,:))];

   h = plot(hax,1:qc.MaxReadLength,100*reads/qc.NumberOfReads,'-o','MarkerFaceColor','Auto');
   set(h,{'MarkerFaceColor'},get(h,'Color'))
   set(h,{'Marker'},{'o';'v';'s';'d'});
   mx = max(mx,qc.MaxReadLength);
end
hold(hax,'off');

xlabel(hax,'Base Position');
ylabel(hax,'% Read Count');
hleg = legend(hax,'0-10','11-20','21-30','31-40');

if ~isscalar(qcs)
    hl = findobj(hleg,'Type','line');
    set(hl,'Color','k');
    set(hl,'MarkerFaceColor','auto');
end

axis(hax,[1 mx 0 100]); 
