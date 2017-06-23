function [ parsed_osm1,counting,connectivity_matrix ] = arrange_map( intersection_node_indices,parsed_osm )
% 2017.6.23 (c) Ofer Keren, ofer293@gmail.com; Itay Levitan, itay007@gmail.com

[counting,intersections_id,parsed_osm1]=split_ways(parsed_osm,intersection_node_indices);
   [connectivity_matrix, intersection_node_indices] = extract_connectivity(parsed_osm1);
      intersection_nodes = get_unique_node_xy(parsed_osm1, intersection_node_indices);

end

