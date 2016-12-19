 function hax = plotPerPositionGC(qc,hax)
%PLOTPERPOSITIONGC Line plot displaying the percentage of nucliotides that
%are either G or C versus position
% 
%   PLOTPERPOSITIONGC(qs) displays a plot showing the percentage of
%   nucliotides at each position that are either G or C.
%
%   h = PLOTPERPOSITIONGC(qs) returns the handle to the axes object
%   containing the plot.
if nargin == 1
    hax = gca;
end


mx = 10;
hold(hax,'all');
for i = 1:numel(qc);
    gc = 100 * (qc(i).PerPosBaseDist(nt2int('G'),:) + qc(i).PerPosBaseDist(nt2int('C'),:))./sum(qc(i).PerPosBaseDist);
    hl = plot(hax,1:qc(i).MaxReadLength,gc,'Marker','o');
    set(hl,'MarkerFaceColor',get(hl,'Color'));
    mx = max(mx,qc(i).MaxReadLength); 
end
hold(hax,'off');
axis(hax,[0 mx 0 110]);
xlabel(hax,'Base Position');
ylabel(hax,'% GC-Content');
title(hax,'Percentage GC vs. Position');