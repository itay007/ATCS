function hax = plotPerPositionQuality(qcs,hax)
%PLOTPERPOSITIONQUALITY Box plot displaying the distribution of Phred
%scores at a given position versus position.
% 
%   PLOTPERPOSITIONQUALITY(qs) displays a box plot showing the median Phred
%   score, the 25 and 75th percentiles, and the most extreme scores that
%   are not considered outliers. 
%
%   A point is not an outlier if it is within 1.5 times the interquartile
%   distance from the median.
%
%   h = PLOTPERPOSITIONQUALITY(qs) returns the handle to the axes object
%   containing the plot.
if nargin == 1
    hax = gca;
end

switch numel(qcs)
    case 1
        colors = {'b' 'r'};
    otherwise
        co = get(hax,'ColorOrder');
        co_dark = co/2;
        colors = mat2cell([co co_dark],ones(size(co,1),1),[3 3]);
end

%average = ((0:qc.MaxPhredScore)*qc.PerPosQualities )./ sum(qc.PerPosBaseDist);%
%plot(hax,average,'g-o','MarkerFaceColor','g');
lim = [0 40];
hold(hax,'on');
for i = 1:numel(qcs)
    qc = qcs(i);
    c = mod(i-1,size(colors,1))+1;
    for j = 1:qc.MaxReadLength
        boxplotFromFreq(hax,qc.PerPosQualities(:,j),j,1/2,qc.MinEncodingPhred,colors{c,1},colors{c,2});
    end
     lim = [min(lim(1),qc.MinScore) max(lim(2),qc.MaxScore)]; 
end

mxread = max(arrayfun(@(q)q.MaxReadLength,qcs));
hold(hax,'off');
axis(hax,[1 mxread lim]);
xlabel(hax,'Base Position');
ylabel(hax,'Average Score');
title(hax,'Quality Average vs. Base Position');




function boxplotFromFreq(hax,freq,x,width,offset,lc,mc)
ind = quantileScoreFromFreq(freq,[.25 .5 .75]);
shift = offset - 1;
q = ind + shift;
height = q(3)-q(1);
if height > 0
rectangle('Position',[x-width/2 q(1) width height],...
          'Parent',hax,...
          'EdgeColor',lc,...
          'FaceColor','none'); 
end
line([x-width/2 x+width/2],[q(2) q(2)],'Parent',hax,'Color',mc);

%%%
range = (1:length(freq))';
upper_ind = ind(3)+1.5*(height);
lower_ind = ind(1)-1.5*(height);


whisker_upper = find( range >=lower_ind & freq ,1,'first')+shift;

whisker_lower = find( range <=upper_ind & freq ,1,'last' )+shift;
line([x x],[whisker_upper whisker_lower],'Parent',hax,'Color',mc);
line( [x-width/3 x+width/3],[whisker_lower whisker_lower],'Parent',hax,'Color',mc);
line( [x-width/3 x+width/3],[whisker_upper whisker_upper],'Parent',hax,'Color',mc);




 
function score = quantileScoreFromFreq(pmf,q)
pmf = pmf(:);
cdf = cumsum(pmf);
score = zeros(length(q),1);
for i = 1:length(q)
    eff_index = q(i)*cdf(end); 
    pot =(cdf >= eff_index) & ((cdf - pmf) <= eff_index);
    score(i) = (find(pot,1,'first') + find(pot,1,'last'))/2;         
end

