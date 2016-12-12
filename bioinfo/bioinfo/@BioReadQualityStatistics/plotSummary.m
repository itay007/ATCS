function h = plotSummary(qc)
% PLOTSUMMARY Summary plots of a BioReadQualityStatistics object
%  PLOTSUMMARY(qs) generates a figure containing six plots that present
%  summary statistics of the data stored in BioReadQualityStatisitcs object
%  qs.    



[w, l] = screenSizePixels();

reduce = 0.8;
ratio = [8.5 11]; %[Width Height]
if w >= l
    p = l/ratio(2);
else
    p = w/ratio(1);
end 
p = floor(reduce*p*ratio);



hf = figure('Position',[0 0 p],'PaperPositionMode','auto','Visible','off');
movegui(hf,'center');
set(hf,'Visible','on');

m = 4; n = 3;
pos = {'plotPerPositionQuality', [1,3]
       'plotPerPositionCountByQuality' [4,6]
       'plotPerPositionGC'  ,[7 9]
       'plotPerSequenceQuality', 10
       'plotPerSequenceGC' , 11
       'plotTotalGC', 12};
h = zeros(size(pos,1),1);
for i = 1:length(h)
     h(i) = subplot(m,n,pos{i,2},'Parent',hf);
     setUpPlot(qc,pos{i,1},h(i));
end
linkaxes(h(1:3),'x');
    
function setUpPlot(qc,plotName,h)
qc.(plotName)(h);
hf = get(h,'parent');
set(h,'ButtonDownFcn',@(~,~)(createNewFigure(hf,qc,plotName)));

function createNewFigure(hf,qc,plotName)
if strcmp(get(hf,'SelectionType'),'normal')
    hax = axes('Parent',figure);
    qc.(plotName)(hax);
end
    
    
function [w, l] = screenSizePixels()
u = get(0,'Units');
set(0,'Units','pixels')
p = get(0,'ScreenSize');
set(0,'Units',u);
w = p(3);
l = p(4);












