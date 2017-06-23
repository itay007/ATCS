% 2017.6.23 (c) Ofer Keren, ofer293@gmail.com; Itay Levitan, itay007@gmail.com

% clear all
[ intersection_nodes,connectivity_matrix ,intersection_node_indices,parsed_osm ] = getmap(  );
%% split ways
[ parsed_osm1,counting,connectivity_matrix ] = arrange_map( intersection_node_indices,parsed_osm);
                                               
%% calc distance
[dists,waynd] = calc_dist_mat(parsed_osm1,counting);
[dg,dists] = w_dist(dists,counting,parsed_osm1,connectivity_matrix);
dg= dg.'+dg;
%% choose the mode you want to plot
disp('1.Block a Path');
disp('2.Traffic Accumulation');
disp('3.Random Mode');
choose=input('Choose Simulation Mode:');

if (choose==1)
    loops=1;
elseif (choose==2)
    loops=2;
elseif(choose==3)
    loops=1;
end
for x=1:loops
    switch choose
        case 1 % Block a Path Mode
%%
                n=1; %number of cars in total
                n1=0; % number of cars for second route
                start_tar(1,1)=input('start :');
                start_tar(2,1)=input('stop :');
                %block
                b_start=input('way to block start :');
                b_stop=input('way to block stop :');
                dg(b_start,b_stop)=0;
                dg(b_stop,b_start)=0;
                %end block

                [route,dist] = plan_first_route(dg,n,n1,start_tar);

        case 2 %% Traffic Accumulation Mode
            %explain -  for n=4 and n1=3 with weight of 1.4 we see all the cars
            %going in the same way from  495 to 484
            %after changing the n1=1 (the number of cars for the route after calc
            %the weights) we can see that the 3 cars for the first road
            %changed the weight of the road, so the last car goes from another
            %road.  --showing that the number of cars for road make traffic
            
            if (x==1)  %no traffic - do not change!!!
              n=2; %number of cars in total
              n1=1; % number of cars for second route
              weight=1.2;  %this weight control the traffic for the roads!
              
              start1 = 495; % node global index
              target1 = 484;
              start2=69;
              
              start_tar(1,1)=start1;
              start_tar(2,1)=target1;
              start_tar(1,2)=start2;
              start_tar(2,2)=target1;
                
              % plan a route for N cars
              [route,dist] = plan_first_route(dg,4,3,start_tar);

              % dg*weights - plan a route for N cars after weights of roads
              [dg] = calc_weight_dist(dg,4,3,route,weight,1);
              [ route ] = plan_second_route( route,dg,2,1,start_tar);

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
 
              % plan a route for N cars
              [route,dist] = plan_first_route(dg,n,n1,start_tar);

              % dg*weights - plan a route for N cars after weights of roads
              [dg] = calc_weight_dist(dg,n,n1,route,weight,1);
              [ route ] = plan_second_route( route,dg,n,n1,start_tar);
            end



        case 3 %Random Mode
            for k1=1:4
                num1=0;
                num2=0;
                while(num1==num2)
                    ran=randi(38);
                    num1=intersection_node_indices(ran);
                    start_tar(1,k1)=num1;
                    ran=randi(38);
                    num2=intersection_node_indices(ran);
                    start_tar(2,k1)=num2;
                end
            end 

            n=4; %number of cars in total
            n1=3; % number of cars for second route
            weight=1.2;  %this weight control the traffic for the roads!
            [route,dist] = plan_first_route(dg,n,n1,start_tar);
             [dg] = calc_weight_dist(dg,n,n1,route,weight,1);
             
            n=2; %number of cars in total
            n1=1; % number of cars for second route
            weight=1.2;  %this weight control the traffic for the roads!
            [ route ] = plan_second_route( route,dg,n,n1,start_tar);
            [dg] = calc_weight_dist(dg,4,2,route,weight,2);

            n=3; %number of cars in total
            n1=1; % number of cars for second route
            weight=1.2;  %this weight control the traffic for the roads!
            [ route ] = plan_second_route( route,dg,n,n1,start_tar);
            [dg] = calc_weight_dist(dg,4,1,route,weight,3);

            n=4; %number of cars in total
            n1=1; % number of cars for second route
            weight=1.2;  %this weight control the traffic for the roads!
            [ route ] = plan_second_route( route,dg,n,n1,start_tar);
            [dg] = calc_weight_dist(dg,4,0,route,weight,4);
    end


    %%plots
    plots( waynd,route,counting,n,n1,parsed_osm, intersection_node_indices,dist)
end