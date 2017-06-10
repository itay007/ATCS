function [ intersection_nodes,connectivity_matrix,intersection_node_indices,parsed_osm ] = getmap(  )
%GETMAP Summary of this function goes here
%   Detailed explanation goes here
%% name file
openstreetmap_filename = 'map_MP1.osm';
%map_img_filename = 'map.png'; % image file saved from online, if available

%% convert XML -> MATLAB struct
% convert the OpenStreetMap XML Data file donwloaded as map.osm
% to a MATLAB structure containing part of the information describing the
% transportation network
[parsed_osm, osm_xml] = parse_openstreetmap(openstreetmap_filename);

%% find connectivity
[connectivity_matrix, intersection_node_indices] = extract_connectivity(parsed_osm);
intersection_nodes = get_unique_node_xy(parsed_osm, intersection_node_indices);


end

