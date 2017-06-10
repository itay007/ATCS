function [ parsed_osm1,counting,connectivity_matrix ] = arrange_map( intersection_node_indices,parsed_osm )
%ARRANGE_MAP Summary of this function goes here
%   Detailed explanation goes here
parsed_osm1=parsed_osm;
[counting,intersections_id,parsed_osm1]=split_ways(parsed_osm1,intersection_node_indices);
%        parsed_osm=parsed_osm1;
% % 
   [connectivity_matrix, intersection_node_indices] = extract_connectivity(parsed_osm1);
      intersection_nodes = get_unique_node_xy(parsed_osm1, intersection_node_indices);

end

