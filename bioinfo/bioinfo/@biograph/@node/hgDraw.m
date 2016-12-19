function hgDraw(h,UIcontextMenu)
%HGDRAW  


% Copyright 2003-2010 The MathWorks, Inc.

hgAxes = h(1).up.hgAxes;

if isempty(hgAxes)
    error(message('bioinfo:node:hgDraw:nohgaxes'));
end

for i = 1:numel(h)
    hp = patch('Visible','off','UIContextMenu',UIcontextMenu,'Parent',hgAxes);
    setappdata(hp,'node',h(i));
    h(i).hgPatch = hp;

    ht = text('Visible','off','horizontalalignment','center',...
        'Clipping','on','Interpreter','none',...
        'UIContextMenu',UIcontextMenu,'Parent',hgAxes);
    setappdata(ht,'node',h(i));
    h(i).hgText = ht;
end

% update rendering with info from node (vectorized)
h.hgUpdate

        
