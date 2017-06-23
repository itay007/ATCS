function [route,dist] = plan_first_route(dg,n,n1,start_tar )
% 2017.6.23 (c) Ofer Keren, ofer293@gmail.com; Itay Levitan, itay007@gmail.com


for i=1:(n-n1)
    if i==1
        start = start_tar(1,1); % node global index
        target =  start_tar(2,1);
    elseif i==2
        start = start_tar(1,2); % node global index
        target =  start_tar(2,2);
    elseif i==3
        start = start_tar(1,3); % node global index
        target =  start_tar(2,3);
    elseif i==4
        start = start_tar(1,4); % node global index
        target =  start_tar(2,4);
    end
     [route1, dist] = route_planner(dg, start, target);
     route{1,i}(1,:)=route1;
end

end

