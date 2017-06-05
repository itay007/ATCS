function [dg,arr_dg] = calc_weight_dist(dg,dists,parsed_osm1,counting)
%% arrange dg as parsed_osm1

ways = parsed_osm1.way;
waynd = ways.nd;










for i=1:size(waynd,2)%loop for each way node struct
%     counter=0;
    for k=1:size(waynd{1,i},2) %loop for each node struct element
        for j=1:size(intersections_id,2)
            if (waynd{1,i}(1,k)==intersections_id(j))
               waynd{1,i}(2,k)=intersection_node_indices(j);
               counter=counter+1;
            end
        end
    end
    if (counter == 0)
        waynd{1,i}(2,1:size(waynd{1,i},2)) = 0;
    end
end
