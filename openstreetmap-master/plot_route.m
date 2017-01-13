function [] = plot_route(ax, route, parsed_osm)
% plot (over map) the route found by route planner
%
% usage
%   PLOT_ROUTE(ax, route, parsed_osm)
%
% input
%   ax = axes object handle.
%   route = matrix of nodes traversed along route, as returned by the
%           route_planner function.
%   parsed_osm = parsed OpenStreetMap XML data, as returned by the
%                parse_openstreetmap function.
%
% 2017.01.12 (c) Itay Levitan, itay007@gmail.com
%
% See also ROUTE_PLANNER, PARSE_OPENSTREETMAP.

% empty path ?
if isempty(route)
    warning('path is empty. This means that the planner found no path.')
    return
end

%get osm data

nodes = parsed_osm.node;
parsed_route_osm.node = parsed_osm.node;
node_ids = nodes.id;
node_xys = nodes.xy;

ways = parsed_osm.way;
node_ways_ids = ways.id;
node_ways_tag = ways.tag;
waynd = ways.nd;



held = takehold(ax);
nodelist = [];

%get passed nodes from selected route
% nodexy = parsed_osm.node.xy;
% start_xy = nodexy(:, route(1, 1) );
% path_xy = nodexy(:, route);
% path_end = nodexy(:, route(1, end) );
%held = takehold(ax);

%%%%%%%%%%%%%%%%%%% algo try 2 %%%%%%%%%%%%%%%%%&%%%

route_ids = node_ids(1,route);
n=1;

%%find route_ids in waynd
%optional_route_ids = [];
for i=1:size(waynd,2)%loop for each way node struct
    for j=1: size(waynd{1,i},2) %loop for each node struct element
     way_size= size(waynd{1,i},2);
    optional_route_ids_cell(1,j) = {find(waynd{1,i}(1,j) == route_ids(1,:))}; %find if there are 2 routed node in way's
    end
    optional_route_ids = cell2mat(optional_route_ids_cell);
    if(size(optional_route_ids,2)>1)
        routed_ways(n) = {i}; % routed ways
        n=n+1;
     end
end

%routed_ways = cell2mat(routed_ways);


%%%%%%%%%%%%%%% route_parsed_osm %%%%%%%%%%%%%%%%%%%%%%

parsed_route_osm.bounds = parsed_osm.bounds;
parsed_route_osm.node = parsed_osm.node;
    
for i=1:size(routed_ways , 2)
    j=cell2mat(routed_ways(i));
    [parsed_route_osm.way.id(1,i)] = node_ways_ids(1 , j);
    [parsed_route_osm.way.nd(1,i)] = waynd(1, j);
    [parsed_route_osm.way.tag(1,i)] = node_ways_tag(1 , j);
    
end

plot_route_way(ax, parsed_route_osm);
