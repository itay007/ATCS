function [] =car_plot(ax, parsed_route_osm)
% plot simulated driving car entity(over the map)
%
% usage
%   car_plot(ax, parsed_route_osm)
%
% input
%
% 2017.04.22 (c) Itay Levitan, itay007@gmail.com


map_img_filename = [];

[bounds, node, way, ~] = assign_from_parsed(parsed_route_osm);
show_ways(ax, bounds, node, way, map_img_filename);

function [] = show_ways(hax, bounds, node, way, map_img_filename)
show_map(hax, bounds, map_img_filename)

key_catalog = {};
for i=1:size(way.id, 2)
    % plot highway
    way_nd_ids = way.nd{1, i};
    num_nd = size(way_nd_ids, 2);
    nd_coor = zeros(2, num_nd);
    nd_ids = node.id;
    for j=1:num_nd
        cur_nd_id = way_nd_ids(1, j);
        if ~isempty(node.xy(:, cur_nd_id == nd_ids))
             nd_coor(:, j) = node.xy(:, cur_nd_id == nd_ids);
        end
    end
    
    % remove zeros
    nd_coor(any(nd_coor==0,2),:)=[];
    
    if ~isempty(nd_coor)
        plot(hax, nd_coor(1,:), nd_coor(2,:), 'g-', 'LineWidth',2)
    end
end
disp(key_catalog.')

%h1 = plot(nd_coor(1,1), nd_coor(2,1),'^','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10); %// plot initial position


% n=size(nd_coor, 2)-1;
% for i=1:n
%     m_nd_coor(1,1) = ((nd_coor(2,2)-nd_coor(2,1))/(nd_coor(1,2)-nd_coor(1,1)));
%     b = -m_nd_coor(1,1)*nd_coor(1,1)+nd_coor(2,1);
%     x = linspace(nd_coor(1,1),nd_coor(1,2),20);
%     %x = nd_coor(1,1):0.0001:nd_coor(1,2);
%     y=m_nd_coor(1,1)*x+b;
%     %plot(x,y);
%     %plot(x(2),y(2),'^','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
% 
%     for i=1:20
%         h2 = plot(x(i),y(i),'*','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
%         pause(0.3);
%         set(h2,'Visible','off')
%     end
% end




%%% one way %%%
m_nd_coor(1,1) = ((nd_coor(2,2)-nd_coor(2,1))/(nd_coor(1,2)-nd_coor(1,1)));
b = -m_nd_coor(1,1)*nd_coor(1,1)+nd_coor(2,1);
x = linspace(nd_coor(1,1),nd_coor(1,2),20);
%x = nd_coor(1,1):0.0001:nd_coor(1,2);
y=m_nd_coor(1,1)*x+b;
%plot(x,y);
%plot(x(2),y(2),'^','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);

for i=1:20
    h2 = plot(x(i),y(i),'*','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
    pause(0.3);
    set(h2,'Visible','off')
end

%%%%%%%%%%%%%%%

% for n = 1:num_nd
%     set(h1, nd_coor(1,n),nd_coor(2,n)); %// update position of object 2
%     drawnow %// refresh figure
% end
%  figure
%  hold on
%  h1 = plot(1:10,'r');
%  h2 = plot(2:11,'g');
%  if <your condition doesn't hold>
%      set(h2,'Visible','off')
%  end



function [] = disp_info(bounds, Nnode, Nway)
disp( ['Bounds: xmin = ' num2str(bounds(1,1)),...
    ', xmax = ', num2str(bounds(1,2)),...
    ', ymin = ', num2str(bounds(2,1)),...
    ', ymax = ', num2str(bounds(2,2)) ] )
disp( ['Number of nodes: ' num2str(Nnode)] )
disp( ['Number of ways: ' num2str(Nway)] )



% t = linspace(0,2,1000); %// time parameter
% x1 = 10*cos(2*pi*t+1);
% y1 = 5*sin(2*pi*t+1); %// trajectory of object 1
% x2 = 2*cos(6*pi*t-2);
% y2 = 3*sin(6*pi*t-2); %// trajectory of object 1
% plot(x1,y1,'color',[.5 .5 .5]); %// plot trajectory of object 1
% hold on
% plot(x2,y2,'color',[.5 .5 .5]); %// plot trajectory of object 2
% h1 = plot(x1(1),y1(1),'ro'); %// plot initial position of object 1
% h2 = plot(x2(1),y2(1),'b*'); %// plot initial position of object 2
% axis([-12 12 -12 12]) %// freeze axis size
% grid on
% for n = 1:numel(t)
%     set(h1, 'XData', x1(n), 'YData', y1(n)); %// update position of object 2
%     set(h2, 'XData', x2(n), 'YData', y2(n)); %// update position of object 2
%     drawnow %// refresh figure
% end