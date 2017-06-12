function [  ] = Plots( waynd,route,counting,n,n1,parsed_osm, intersection_node_indices)
%PLOTS Summary of this function goes here
%   Detailed explanation goes here

%% XY of routes for one car
%[routeXY] = XY_Route(waynd,route1,counting);
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

plot_nodes(ax, parsed_osm, intersection_node_indices)

%plot_route_new(routeXY)  %plot for one car
plot_route_new_n(routeXY,n)  %plot for n car

only_nodes = 1:10:1000; % not all nodes, to reduce graphics memory & clutter
plot_nodes(ax, parsed_osm, only_nodes)

%show intersection nodes (unseful, but may result into a cluttered plot)
plot_nodes(ax, parsed_osm, intersection_node_indices)

hold(ax, 'off')

end

