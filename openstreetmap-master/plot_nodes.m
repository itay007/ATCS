function [] = plot_nodes(ax, parsed_osm, only_node_indices,route)
% plot (selected) nodes and label each with index and id
%
% usage
%   PLOT_NODES(ax, parsed_osm, only_node_indices, show_id)
%
% input
%   ax = axes object handle where to plot the nodes.
%   parsed_osm = MATLAB structure containing the OpenStreetMap XML data
%                after parsing, as returned by function
%                parse_openstreetmap.
%   only_node_indices = selected node indices in the global node matrix.
%   show_id = select whether to show or not the ID numbers of nodes as text
%             labes within the plot
%           = 0 (do not show labels) | 1 (show labels)
%
% 2012.04.24 (c) Ioannis Filippidis, jfilippidis@gmail.com
%
% See also PARSE_OPENSTREETMAP, ROUTE_PLANNER.

% do not show node id (default)
if nargin < 5
    show_id = 0;
end

nodes = parsed_osm.node;
node_ids = nodes.id;
node_xys = nodes.xy;

% which nodes to plot ?
n = size(node_xys, 2);
if nargin < 3
    only_node_indices = 1:n;
end

%% plot
held = takehold(ax);

% nodes selected exist ?
if n < max(only_node_indices)
    warning('only_node_indices contains node indices which are too large.')
    return
end

% plot nodes
xy = node_xys(:, only_node_indices);
plotmd(ax, xy, 'yo')

% label plots
for i=only_node_indices
    node_id_txt = num2str(node_ids(1, i) );
    if show_id
        curtxt = {['index=', num2str(i) ], ['id=', node_id_txt] }.';
    else
        curtxt = ['index=', num2str(i) ];
    end
    colors={'Color','black','black','black','black','black'};
    road_color=2;
    sizer=size(route,2);
    if(sizer>0)
        if(i==route{1,1}(1,1) || i==route{1,1}(1,size(route{1,1},2)))
            road_color=3;
        end
    end
    if(sizer>1)
        if(i==route{1,2}(1,1) || i==route{1,2}(1,size(route{1,2},2)))
            road_color=4;
        end
    end
    if(sizer>2)
        if(i==route{1,3}(1,1) || i==route{1,3}(1,size(route{1,3},2)))
            road_color=5;
        end
    end
    if(sizer>3)
        if(i==route{1,4}(1,1) || i==route{1,4}(1,size(route{1,4},2)))
            road_color=6;
        end
    end
    textmd(node_xys(:, i), curtxt, 'Parent', ax,colors(1),colors(road_color))
end

restorehold(ax, held)
