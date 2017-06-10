function [ routeXY ] = XY_Route( waynd,route1,counting )
%XY_ROUTE Summary of this function goes here
%   Detailed explanation goes here

z=1;
for i=1:(size(route1,2)-1)
    for j=1:size(waynd,2)
        if (counting(1,j)~=0 && counting(2,j)~=0)
            if route1(i) == waynd{1,j}(2,counting(1,j)) 
               if route1(i+1) == waynd{1,j}(2,counting(2,j))

                routeXY(1,z:z+counting(2,j)-1)= waynd{1,j}(3,1:counting(2,j));  %% X cord
                routeXY(2,z:z+counting(2,j)-1)= waynd{1,j}(4,1:counting(2,j)); %% Y cord
                z=z+counting(2,j)-1;
               end
            end
        end
    end
end
end

