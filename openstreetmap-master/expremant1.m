% 2017.6.23 (c) Ofer Keren, ofer293@gmail.com; Itay Levitan, itay007@gmail.com

% clear all
[ intersection_nodes,connectivity_matrix ,intersection_node_indices,parsed_osm ] = getmap(  );
%% split ways
[ parsed_osm1,counting,connectivity_matrix ] = arrange_map( intersection_node_indices,parsed_osm);
                                               
%% calc distance
[dists,waynd] = calc_dist_mat(parsed_osm1,counting);
[dg,dists] = w_dist(dists,counting,parsed_osm1,connectivity_matrix);
dg= dg.'+dg;


dg(860,495)=0;
dg(865,866)=0;
dg(485,453)=0;


start=860;
target=484;

weight=1.01;
for i=1:100
     [route1, dist] = route_planner(dg, start, target);
     route{1,i}(1,:)=route1;
    
    for j=1: (size(route{1,i}, 2)-1)
        first = route{1,i}(1,j);
        last = route{1,i}(1,j+1);
        dg(first,last) = dg(first,last)*weight;
        A{1,i}(1,:)=dg(860,484);
        B{1,i}(1,:)=dg(860,865)+dg(865,485)+dg(485,484);
        if(dg(first,last) >1)
            dg(first,last)=0.99;
        end
    end
end

A1 = cell2mat(A);
B1 = cell2mat(B);

plot(A1,'gs',...
    'LineWidth',1,...
    'MarkerSize',2,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5])


title('dg of A(red) vs B(green)')
hold on;
plot(B1,'gs',...
    'LineWidth',1,...
    'MarkerSize',2,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',[0.5,0.5,0.5])





