function hgUpdate(hin)
%HGUPDATE

% Copyright 2003-2012 The MathWorks, Inc.

minAllowableFontSize = 1.5;

if isempty(hin(1).hgLine)
    error(message('bioinfo:edge:hgUpdate:nohgline'));
end

bg = hin(1).up;
lineType = bg.EdgeType;
sc = bg.Scale;
showArrows = isequal(bg.showArrows,'on');
showWeights = isequal(bg.showWeights,'on');
arrowSize = bg.ArrowSize;

xa = zeros(numel(hin),3);
ya = zeros(numel(hin),3);
axesScale = min(1,getappdata(bg.hgAxes,'Scale'));
edgeFontSize =  max(bg.EdgeFontSize * axesScale,minAllowableFontSize);

for i = 1:numel(hin)
    h = hin(i);

    xy = h.ControlPoints;

    xys = xy(:,1)+(sc-1)*h.FromNode.Position';
    xye = xy(:,end)+(sc-1)*h.ToNode.Position';
    switch lineType
        case 'straight'
            xy = [xys xye];
        otherwise % 'curved' &&  'segmented'
            useStraight = false;
            np = size(xy,2); % number of points in the vector
            nc = (np+1)/3;   % number of control points
            dxy = xy(:,np)-xy(:,1); % x & y dist between start end end of edge
            if all(abs(dxy)>1) % then we can guess the rescaling for both
                % independently
                lsc = (xye-xys)./dxy;
            else % one of the distances is close to zero, we rely on the scale
                % derived from the other dimension
                co = diff(abs(dxy));
                if co<0
                    lsc = repmat((xye(1)-xys(1))/dxy(1),2,1);
                elseif co>0
                    lsc = repmat((xye(2)-xys(2))/dxy(2),2,1);
                elseif ~any(any(diff(xy,[],2))) % both seem to be 0 !
                    lsc = [1;1];
                else
                    useStraight = true;
                    xy = [xys xye];
                end
            end
            if ~useStraight
                if isequal(lineType,'segmented')
                    cpi = [1:3:np np];
                    % now re-scale the control points to the scale we have
                    xy = repmat(xys,1,nc+1) + ...
                        (xy(:,cpi)-repmat(xy(:,1),1,nc+1)).*repmat(lsc,1,nc+1);
                else % 'Bezier splines'
                    t = 0:0.05:1; nt=21;
                    % now re-scale the control points to the scale we have
                    scp = repmat(xys,1,np) + ...
                        (xy-repmat(xy(:,1),1,np)).*repmat(lsc,1,np);
                    xy = zeros(2,numel(t)*(nc-1)+1);
                    j = 1;
                    % do Bezier Spline for every segment
                    for l = 1:nc-1
                        xy(:,(l-1)*nt+1:l*nt) =  scp(:,j)*((1-t).^3) + ...
                            3*scp(:,j+1)*(t.*(1-t).^2) + ...
                            3*scp(:,j+2)*((1-t).*t.^2) + ...
                            scp(:,j+3)*t.^3;
                        j = j+3;
                    end
                    % add arrow (or last) point
                    xy(:,end) = scp(:,end);
                end
            end
    end
    
    xa(i,:) = xy(1,[1 end-1 end]);  % use later to figure out arrows
    ya(i,:) = xy(2,[1 end-1 end]);

    if h.Visible
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
    
    set(h.hgLine,...
        'XData',xy(1,:),...
        'YData',xy(2,:),...
        'Color',LineColor,...
        'LineWidth',LineWidth,...
        'Visible',vis);

    if showWeights
        sxy = size(xy,2);
        if rem(sxy,2)
            txy = xy(:,(sxy+1)/2);
            qxy = diff(xy(:,(sxy+1)/2+[-1 1]),[],2);
        else
            txy = sum(xy(:,(sxy+[0 2])/2),2)/2;
            qxy = diff(xy(:,(sxy+[0 2])/2),[],2);
        end
        if isempty(get(h.hgline,'UserData'))
            set(h.hgline,'UserData',text(txy(1),txy(2),[' ' num2str(h.Weight) ' '] ,...
                'FontSize',edgeFontSize,'Visible',vis,...
                'Color',bg.edgeTextColor,'Parent',bg.hgAxes));
        else
            set(get(h.hgline,'UserData'),...
                'Position',[txy' 0],...
                'String',[' ' num2str(h.Weight) ' '],...
                'FontSize',edgeFontSize,'Visible',vis,'Color',bg.edgeTextColor);
        end
        if diff(abs(qxy))>0
            if prod(qxy)>0
                set(get(h.hgline,'UserData'),... 
                    'VerticalAlignment','middle','HorizontalAlignment','right')
            else
                set(get(h.hgline,'UserData'),... 
                    'VerticalAlignment','middle','HorizontalAlignment','left')
            end
        else
            if prod(qxy)>0
                set(get(h.hgline,'UserData'),... 
                    'VerticalAlignment','top','HorizontalAlignment','center')
            else
                set(get(h.hgline,'UserData'),... 
                    'VerticalAlignment','bottom','HorizontalAlignment','center')
            end
        end      
    elseif ~isempty(get(h.hgline,'UserData'))
        set(get(h.hgline,'UserData'),'Visible','off')
    end
end

for i = 1:numel(hin)
    h = hin(i);

    xy = [xa(i,:);ya(i,:)];
        
    % draw arrow heads when necessary
    if showArrows
        if all(abs(diff(xy(:,[end-1 end]),[],2))<sqrt(eps)) % both seem to be 0 !
            xy = xy(:,[1 end]);
        end
        if all(abs(diff(xy,[],2))<sqrt(eps)) % still zero ?
            xya = xy(:,end);
        else
            s = arrowSize; % arrow size in pts
            tip = xy(:,end);
            uv = (xy(:,end-1)-tip)/sqrt(sum((xy(:,end-1)-tip).^2));
            xya =[tip tip+[s s/4;-s/4 s]*uv tip+s*5/6*uv tip+[s -s/4;s/4 s]*uv];
        end
        if h.Visible
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
        set(h.hgArrow,...
            'XData',xya(1,:),...
            'YData',xya(2,:),...
            'FaceColor',LineColor,...
            'EdgeColor',LineColor,...
            'LineStyle','-',...
            'Visible',vis,...
            'LineWidth',LineWidth); 
    else
        set(h.hgArrow,'Visible','off')
    end
end

