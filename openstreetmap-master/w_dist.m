function [dg,dists] = w_dist(dists,counting,parsed_osm1,connectivity_matrix)
% 2017.6.23 (c) Ofer Keren, ofer293@gmail.com; Itay Levitan, itay007@gmail.com

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
        dg(dists(2,i),dists(3,i)) = dists(1,i);
    end
end

