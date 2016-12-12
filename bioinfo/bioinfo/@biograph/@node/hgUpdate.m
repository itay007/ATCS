function hgUpdate(hin)
%HGUPDATE

% Copyright 2003-2012 The MathWorks, Inc.

minAllowableFontSize = 1.5;

if isempty(hin(1).hgPatch)
    error(message('bioinfo:node:hgUpdate:nohgpatch'));
end

bg = hin(1).up;
biographScale = bg.Scale;
axesScale = min(1,getappdata(bg.hgAxes,'Scale'));
doCustomDraw = ~isempty(bg.CustomNodeDrawFcn);
if doCustomDraw
    cmenu = getappdata(bg.up.hgFigure,'biographContextMenus');
end


% Update Patches
for i = 1:numel(hin)
    h = hin(i);
    
    rx = h.Size(1);
    ry = h.Size(2);

    x = h.Position(1)*biographScale;
    y = h.Position(2)*biographScale;

    switch h.Shape
        case {'box','rect','rectangle'}
            px = x + rx/2 * [-1 1 1 -1 -1];
            py = y + ry/2 * [-1 -1 1 1 -1];
        case {'circle','ellipse'}
            t = 0:pi/30:2*pi;
            px = rx/2*cos(t)+x;
            py = ry/2*sin(t)+y;
        case 'house'
            px = x + rx/2 * [-1 -1  0 1   1  -1];
            py = y + ry/2 * [-.833 1/3 1 1/3 -.833 -.833];
        case 'invhouse'
            px = x + rx/2 * [-1 -1  0 1   1  -1];
            py = y + ry/2 * [.833 -1/3 -1 -1/3 .833 .833];
        case 'trapezium'
            px = x + rx/2 * [-1 -.6 .6 1 -1];
            py = y + ry/2 * [-1 1 1 -1 -1];
        case 'invtrapezium'
            px = x + rx/2 * [-1 -.6 .6 1 -1];
            py = y + ry/2 * [1 -1 -1 1 1];
        case 'parallelogram'
            px = x + rx/2 * [-1 -.6 1 .6 -1];
            py = y + ry/2 * [-1 1 1 -1 -1];
        case {'diamond'}
            px = x + rx/2 * [0 -1 0 1 0];
            py = y + ry/2 * [-1 0 1 0 -1];
        otherwise
            error(message('bioinfo:node:hgUpdate:InvalidShape', h.Shape));
    end

    if h.Visible && ~doCustomDraw
        vis = 'on';
    else
        vis = 'off';
    end

    if h.isSelected
        LineColor = [1 0.2 0.2];
        LineWidth = max(2,h.LineWidth);
    else
        LineColor = h.LineColor;
        LineWidth = h.LineWidth;
    end

    set(h.hgPatch,...
        'XData',px,...
        'YData',py,...
        'FaceColor',h.Color,...
        'EdgeColor',LineColor,...
        'LineWidth',LineWidth,...
        'Visible',vis);

    if ~isempty(get(h.hgPatch,'UserData'))
        delete(get(h.hgPatch,'UserData'))
        set(h.hgPatch,'UserData',[])
    end
    if doCustomDraw
         set(h.hgPatch,'UserData',bg.CustomNodeDrawFcn(h))
         if h.isSelected
             set(findobj(get(h.hgPatch,'UserData'),'Type','Patch'),'EdgeColor',[1 0 0])
             set(findobj(get(h.hgPatch,'UserData'),'Type','Line'),'Color',[1 0 0])
         end
         ud = get(h.hgPatch,'UserData');
         set(ud,'UIContextMenu',cmenu(1))
         for i = 1:numel(ud)
            setappdata(ud(i),'node',h)
         end
    end
end

if ~doCustomDraw

    % Update Text
    for i = 1:numel(hin)
        h = hin(i);

        x = h.Position(1)*biographScale;
        y = h.Position(2)*biographScale;
        FontSize =  h.FontSize * axesScale;

        if h.Visible
            vis = 'on';
        else
            vis = 'off';
        end

        switch bg.ShowTextInNodes
            case 'label'
                label = h.label;
            case 'id'
                label = h.id;
            case 'none'
                label = '';
        end

        if minAllowableFontSize>FontSize
            set(h.hgText,...
                'Color',h.TextColor,...
                'Position',[x,y,0],...
                'String',label,...
                'FontSize',minAllowableFontSize,...
                'Visible','off');
        else
            set(h.hgText,...
                'Color',h.TextColor,...
                'Position',[x,y,0],...
                'String',label,...
                'FontSize',FontSize,...
                'Visible',vis);
        end
    end

end
