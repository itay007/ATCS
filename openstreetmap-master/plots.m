function [  ] = Plots( waynd,route,counting,n,n1,parsed_osm, intersection_node_indices,dist)
% 2017.6.23 (c) Ofer Keren, ofer293@gmail.com; Itay Levitan, itay007@gmail.com


%% XY of routes for n cars
for i=1:n
    [routeXY1] = XY_Route(waynd,route{1,i},counting,i,n,n1);
    routeXY{1,i}(1,:) = routeXY1(1,:);
    routeXY{1,i}(2,:) = routeXY1(2,:);
end


%% plot
fig = figure(1);
ax = axes('Parent', fig);
hold(ax, 'on')

%% plot the network, optionally a raster image can also be provided for the
% map under the vector graphics of the network
plot_way(ax, parsed_osm)
plot_nodes(ax, parsed_osm, intersection_node_indices,route)
plot_route(routeXY,n)  %plot for n car

only_nodes = 1:10:1000; % not all nodes, to reduce graphics memory & clutter
plot_nodes(ax, parsed_osm, only_nodes)

%show intersection nodes (unseful, but may result into a cluttered plot)
plot_nodes(ax, parsed_osm, intersection_node_indices,route)

hold(ax, 'off')

end

