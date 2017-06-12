% example script for using the OpenStreetMap functions
%
% use the example map.osm file in the release on github
%
% or
%
% download an OpenStreetMap XML Data file (extension .osm) from the
% OpenStreetMap website:
%   http://www.openstreetmap.org/
% after zooming in the area of interest and using the "Export" option to
% save it as an OpenStreetMap XML Data file, selecting this from the
% "Format to Export" options. The OSM XML is specified in:
%   http://wiki.openstreetmap.org/wiki/.osm
%
% See also PARSE_OPENSTREETMAP, PLOT_WAY, EXTRACT_CONNECTIVITY,
%          GET_UNIQUE_NODE_XY, ROUTE_PLANNER, PLOT_ROUTE, PLOT_NODES.
%
% 2010.11.25 (c) Ioannis Filippidis, jfilippidis@gmail.com
% clear all
[ intersection_nodes,connectivity_matrix,intersection_node_indices,parsed_osm ] = getmap(  );
%% split ways
[ parsed_osm1,counting,connectivity_matrix ] = arrange_map( intersection_node_indices,parsed_osm);
                                               
%% calc distance
[dists,waynd] = calc_dist_mat(parsed_osm1,counting);
[dg,dists] = w_dist(dists,counting,parsed_osm1,connectivity_matrix);

%% choose the mode you want to plot
choose=2;
for x=1:2
    switch choose
        case 1 % block a way
                    %%
            if (x==1)
              n=1; %number of cars in total
              n1=0; % number of cars for second route
              weight=1;  %this weight control the traffic for the roads!
            elseif (x == 2)
              n=1; %number of cars in total
              n1=0; % number of cars for second route
              weight=2;  %this weight control the traffic for the roads! big number= more traffic!
                           %input 0 for total block for all the roads on this path!
              [dg] = calc_weight_dist(dg,n,n1,route,weight);
            end
            
            start_tar(1,1:3)=495;
            start_tar(2,1:3)=484;
            [route,dist] = plan_first_route(dg,n,n1,start_tar);


        case 2 % route cars to a new road
            %explain -  for n=4 and n1=3 with weight of 1.4 we see all the cars
            %going in the same way from  495 to 484
            %after changing the n1=1 (the number of cars for the route after calc
            %the weights) we can see that the 3 cars for the first road
            %changed the weight of the road, so the last car goes from another
            %road.  --showing that the number of cars for road make traffic


            %%
            if (x==1)  %no traffic - do not change!!!
              n=4; %number of cars in total
              n1=3; % number of cars for second route
              weight=1;  %this weight control the traffic for the roads!
              
              start1 = 495; % node global index
              target1 = 484;
              start2=69;
              
              start_tar(1,1)=start1;
              start_tar(2,1)=target1;
               start_tar(1,2:4)=start2;
              start_tar(2,2:4)=target1;
            elseif(x==2)  %traffic - do not change!!!
              n=4; %number of cars in total
              n1=1; % number of cars for second route
              weight=1.2;  %this weight control the traffic for the roads!
              start1 = 495; % node global index
              target1 = 484;
              start2=69;
              
              start_tar(1,1:3)=start1;
              start_tar(2,1:3)=target1;
              start_tar(1,4)=start2;
              start_tar(2,4)=target1;
            end
            %% plan a route1 for one car
            % start = 69; % node global index
            % target =484;
            % [route1, dist] = route_planner(dg, start, target);
            %% plan a route2 for one car
            % start = 102; % node global index
            % target =505;
            % [route2, dist] = route_planner(dg, start, target);
            %% plan a route for N cars
            [route,dist] = plan_first_route(dg,n,n1,start_tar);

            % dg*weights - plan a route for N cars after weights of roads
            [dg] = calc_weight_dist(dg,n,n1,route,weight);
            [ route ] = plan_second_route( route,dg,n,n1,start_tar);



        case 3 %route cars random --- what do we do here? ideas!
    end


    %%plots
    plots( waynd,route,counting,n,n1,parsed_osm, intersection_node_indices)
end