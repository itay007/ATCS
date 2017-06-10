function [] =car_plot(ax, parsed_route_osm1, parsed_route_osm2)
% plot simulated driving car entity(over the map)
%
% usage
%   car_plot(ax, parsed_route_osm)
%
% input
%
% 2017.04.22 (c) Itay Levitan, itay007@gmail.com


map_img_filename = [];
%% route1
[bounds, node, way, ~] = assign_from_parsed(parsed_route_osm1);
%%show_ways(ax, bounds, node, way, map_img_filename);
hax=ax;
%%function [] = show_ways(hax, bounds, node, way, map_img_filename)
show_map(hax, bounds, map_img_filename)

key_catalog = {};
for i=1:size(way.id, 2)
    % plot highway
    way_nd_ids = way.nd{1, i};
    num_nd = size(way_nd_ids, 2);
    nd_coor1 = zeros(2, num_nd);
    nd_ids = node.id;
    for j=1:num_nd
        cur_nd_id = way_nd_ids(1, j);
        if ~isempty(node.xy(:, cur_nd_id == nd_ids))
             nd_coor1(:, j) = node.xy(:, cur_nd_id == nd_ids);
        end
    end
    
    % remove zeros
    nd_coor1(any(nd_coor1==0,2),:)=[];
    
    if ~isempty(nd_coor1)
        plot(hax, nd_coor1(1,:), nd_coor1(2,:), 'g-', 'LineWidth',2)
    end
end
disp(key_catalog.')

%% route2
[bounds2, node2, way2, ~] = assign_from_parsed(parsed_route_osm2);
%show_ways(ax, bounds, node, way, map_img_filename);

%function [] = show_ways(hax, bounds, node, way, map_img_filename)
%show_map(hax, bounds, map_img_filename)

key_catalog2 = {};
for i=1:size(way2.id, 2)
    % plot highway
    way_nd_ids2 = way2.nd{1, i};
    num_nd2 = size(way_nd_ids2, 2);
    nd_coor2 = zeros(2, num_nd2);
    nd_ids2 = node2.id;
    for j=1:num_nd2
        cur_nd_id2 = way_nd_ids2(1, j);
        if ~isempty(node2.xy(:, cur_nd_id2 == nd_ids2))
             nd_coor2(:, j) = node2.xy(:, cur_nd_id2 == nd_ids2);
        end
    end
    
    % remove zeros
    nd_coor2(any(nd_coor2==0,2),:)=[];
    
    if ~isempty(nd_coor2)
        plot(hax, nd_coor2(1,:), nd_coor2(2,:), 'g-', 'LineWidth',2)
    end
end
disp(key_catalog2.')
%%
%h1 = plot(nd_coor(1,1), nd_coor(2,1),'^','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10); %// plot initial position
i1=1;
i2=1;
% %% route1+ route2
% for k=1: max(size(way.id, 2),size(way2.id, 2))
% if (i1 <= size(way.id, 2))
% %for i=1:size(way.id, 2)
%     % plot highway
%     way_nd_ids = way.nd{1, i1};
%     num_nd = size(way_nd_ids, 2);
%     nd_coor1 = zeros(2, num_nd);
%     nd_ids = node.id;
%     for j1=1:num_nd
%         cur_nd_id = way_nd_ids(1, j1);
%         if ~isempty(node.xy(:, cur_nd_id == nd_ids))
%              nd_coor1(:, j1) = node.xy(:, cur_nd_id == nd_ids);
%         end
%     end
% end
% if (i2 <= size(way2.id, 2))
% %for i=1:size(way.id, 2)
%     % plot highway
%     way_nd_ids2 = way2.nd{1, i2};
%     num_nd2 = size(way_nd_ids2, 2);
%     nd_coor2 = zeros(2, num_nd2);
%     nd_ids2 = node.id;
%     for j2=1:num_nd2
%         cur_nd_id2 = way_nd_ids2(1, j2);
%         if ~isempty(node.xy(:, cur_nd_id2 == nd_ids2))
%              nd_coor2(:, j2) = node.xy(:, cur_nd_id2 == nd_ids2);
%         end
%     end
% end
%     % remove zeros
%     nd_coor1(any(nd_coor1==0,2),:)=[];
%     n1=size(nd_coor1, 2)-1;
%     nd_coor2(any(nd_coor2==0,2),:)=[];
%     n2=size(nd_coor2, 2)-1;
%     t1=1;
%     t2=1;
%     for t=1:max(n1,n2)
%         if(t1 <=n1)
%             m_nd_coor1(1,t1) = ((nd_coor1(2,t1+1)-nd_coor1(2,t1))/(nd_coor1(1,t1+1)-nd_coor1(1,t1)));
%             b = -m_nd_coor1(1,t1)*nd_coor1(1,t1)+nd_coor1(2,t1);
%             x1(1,:) = linspace(nd_coor1(1,t1),nd_coor1(1,t1+1),5);
%             %x = nd_coor(1,1):0.0001:nd_coor(1,2);
%             y1=m_nd_coor1(1,t1)*x1+b;
%             %plot(x,y);
%             %plot(x(2),y(2),'^','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
%             
%         end
%         if(t2 <=n2)
%             m_nd_coor2(1,t2) = ((nd_coor2(2,t2+1)-nd_coor2(2,t2))/(nd_coor2(1,t2+1)-nd_coor2(1,t2)));
%             b = -m_nd_coor2(1,t2)*nd_coor2(1,t2)+nd_coor2(2,t2);
%             x2(1,:) = linspace(nd_coor2(1,t2),nd_coor2(1,t2+1),5);
%             %x = nd_coor(1,1):0.0001:nd_coor(1,2);
%             y2=m_nd_coor2(1,t2)*x2+b;
%             
%         end
%         for z=1:5
%             if(t1 <=n1)
%                 h1 = plot(x1(z),y1(z),'*','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
%                 pause(0.3);
%                 set(h1,'Visible','off')
%             end
%             if(t2 <=n2)
%                 h2 = plot(x2(z),y2(z),'*','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
%                 pause(0.3);
%                 set(h2,'Visible','off')
%             end 
%             h1 = plot(x1(z),y1(z),'*','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
%             h2 = plot(x2(z),y2(z),'*','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
%         end
%         t1=t1+1;
%         t2=t2+1;
%     end
%     i1=i1+1;
%     i2=i2+1;
% end

%% route1
for i=1:size(way.id, 2)
    % plot highway
    way_nd_ids = way.nd{1, i};
    num_nd = size(way_nd_ids, 2);
    nd_coor1 = zeros(2, num_nd);
    nd_ids = node.id;
    for j=1:num_nd
        cur_nd_id = way_nd_ids(1, j);
        if ~isempty(node.xy(:, cur_nd_id == nd_ids))
             nd_coor1(:, j) = node.xy(:, cur_nd_id == nd_ids);
        end
    end
    
    % remove zeros
    nd_coor1(any(nd_coor1==0,2),:)=[];
    n=size(nd_coor1, 2)-1;

    for i1=1:n
        m_nd_coor1(1,i1) = ((nd_coor1(2,i1+1)-nd_coor1(2,i1))/(nd_coor1(1,i1+1)-nd_coor1(1,i1)));
        b = -m_nd_coor1(1,i1)*nd_coor1(1,i1)+nd_coor1(2,i1);
        x(1,:) = linspace(nd_coor1(1,i1),nd_coor1(1,i1+1),5);
        %x = nd_coor(1,1):0.0001:nd_coor(1,2);
        y=m_nd_coor1(1,i1)*x+b;
        %plot(x,y);
        %plot(x(2),y(2),'^','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);

        for z=1:5
            h2 = plot(x(z),y(z),'*','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
            pause(0.3);
            set(h2,'Visible','off')
        end
    end
end

%% itay
% n=size(nd_coor, 2);
% 
% for i=1:n
%     m_nd_coor(1,i) = ((nd_coor(2,i+1)-nd_coor(2,i))/(nd_coor(1,i+1)-nd_coor(1,i)));
%     b = -m_nd_coor(1,i)*nd_coor(1,i)+nd_coor(2,i);
%     x(1,:) = linspace(nd_coor(1,i),nd_coor(1,i+1),5);
%     %x = nd_coor(1,1):0.0001:nd_coor(1,2);
%     y=m_nd_coor(1,i)*x+b;
%     %plot(x,y);
%     %plot(x(2),y(2),'^','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
% 
%     for z=1:5
%         h2 = plot(x(z),y(z),'*','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
%         pause(0.3);
%         set(h2,'Visible','off')
%     end
% end


function [] = disp_info(bounds, Nnode, Nway)
disp( ['Bounds: xmin = ' num2str(bounds(1,1)),...
    ', xmax = ', num2str(bounds(1,2)),...
    ', ymin = ', num2str(bounds(2,1)),...
    ', ymax = ', num2str(bounds(2,2)) ] )
disp( ['Number of nodes: ' num2str(Nnode)] )
disp( ['Number of ways: ' num2str(Nway)] )


