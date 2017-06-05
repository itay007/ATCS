function [dg] = w_dist(dists,counting,parsed_osm1,connectivity_matrix)
%% arrange dg as parsed_osm1

ways = parsed_osm1.way;
waynd = ways.nd;

dg = connectivity_matrix;




for i=1:size(counting,2)
    if(counting(2,i)~=0)
        first = counting(1,i);
        last = counting(2,i);
         dists(2,i) = waynd{1,i}(2,first);
         dists(3,i) = waynd{1,i}(2,last);
        dg(dists(2,i),dists(3,i)) = dists(1,i,1);
    end
end






% for i=1:size(waynd,2)%loop for each way node struct
% %     counter=0;
%     for k=1:size(waynd{1,i},2) %loop for each node struct element
%         for j=1:size(intersections_id,2)
%             if (waynd{1,i}(1,k)==intersections_id(j))
%                waynd{1,i}(2,k)=intersection_node_indices(j);
%                counter=counter+1;
%             end
%         end
%     end
%     if (counter == 0)
%         waynd{1,i}(2,1:size(waynd{1,i},2)) = 0;
%     end
% end
