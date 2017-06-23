function [dists,waynd] = calc_dist_mat(parsed_osm1,counting)
% 2017.6.23 (c) Ofer Keren, ofer293@gmail.com; Itay Levitan, itay007@gmail.com


nodes = parsed_osm1.node;
node_ids = nodes.id;
node_xys = nodes.xy;
ways = parsed_osm1.way;
node_ways_ids = ways.id;
waynd = ways.nd;
waytag= ways.tag;

dists = zeros(1,54);

for i=1:size(waynd,2)
    mikumid=0;
    for j=1:size(waynd{1,i},2)
        if(waynd{1,i}(1,j) ~= 0)
           for k=1:size(node_ids,2)
               if (node_ids(k) == waynd{1,i}(1,j))
                   mikum = k;
               end
           end
           waynd{1,i}(3,j) = node_xys(1,mikum);      %X cord
           waynd{1,i}(4,j) = node_xys(2,mikum);      %Y cord
        end
    end
end

       for i=1:size(waynd,2) 
            d=0;
            for j=1:size(waynd{1,i},2)-1
                if (counting(2,i) ~= 0)
                    p1=waynd{1,i}(3:4,j);
                    p2=waynd{1,i}(3:4,j+1);
                    if (p2(1) ~=0 && p2(2)~=0)
                         d=d+norm(p1-p2);
                    end
                end
            end
            dists(i)=d;
       end