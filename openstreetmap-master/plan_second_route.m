function [ route ] = plan_second_route( route,dg,n,n1 )
%PLAN_SECOND_ROUTE Summary of this function goes here
%   Detailed explanation goes here

for i=(n-n1+1):n
     if i==1
        start = 495; % node global index
        target = 484;
    elseif i==2
        start = 495; % node global index
        target = 484;
    elseif i==3
        start = 495; % node global index
        target = 484;
    elseif i==4
        start = 495; % node global index
        target = 484;
    end
     [route1, dist] = route_planner(dg, start, target);
     route{1,i}(1,:)=route1;
end
end

