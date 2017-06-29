% 2017.6.23 (c) Ofer Keren, ofer293@gmail.com; Itay Levitan, itay007@gmail.com

% clear all
[ intersection_nodes,connectivity_matrix ,intersection_node_indices,parsed_osm ] = getmap(  );
%% split ways
[ parsed_osm1,counting,connectivity_matrix ] = arrange_map( intersection_node_indices,parsed_osm);
                                               
%% calc distance
[dists,waynd] = calc_dist_mat(parsed_osm1,counting);
[dg,dists] = w_dist(dists,counting,parsed_osm1,connectivity_matrix);
dg= dg.'+dg;


%% self calculate

totalSelf=0;
for i=1:45
    for j=1: (size(routeSelf{1,i}, 2)-1)
         first = routeSelf{1,i}(1,j);
         last = routeSelf{1,i}(1,j+1);
         totalSelf=totalSelf+dg(first,last);
    end
end

totalPc=0;
for i=1:45
    start=ways_to_route(i,1);
    target=ways_to_route(i,2);
    [route1, dist] = route_planner(dg, start, target);
    routePc{1,i}(1,:)=route1;
    routePc{1,i}(2,:)=dist;
    totalPc=dist+totalPc;
end

totalSelf=totalSelf/45;
totalPc=totalPc/45;

disp(totalSelf);
disp(totalPc);