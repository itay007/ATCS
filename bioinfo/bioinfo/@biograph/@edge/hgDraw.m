function hgDraw(h,UIcontextMenu)
%HGDRAW  


% Copyright 2003-2010 The MathWorks, Inc.

hgAxes = h(1).up.hgAxes;

if isempty(hgAxes)
    error(message('bioinfo:edge:hgDraw:nohgaxes'));
end

for i = 1:numel(h)
    hl = line('Visible','off','UIContextMenu',UIcontextMenu,'Parent',hgAxes);
    setappdata(hl,'edge',h(i));
    h(i).hgLine = hl;
end

for i = 1:numel(h)
    ha = patch('LineStyle','none','Visible','off','UIContextMenu',UIcontextMenu,'Parent',hgAxes);
    setappdata(ha,'edge',h(i));
    h(i).hgArrow = ha;
end
                         
% update rendering with info from node
h.hgUpdate
