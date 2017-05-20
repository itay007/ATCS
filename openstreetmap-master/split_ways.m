function [intersections_id,parsed_osm1] = split_ways( parsed_osm1,intersection_node_indices)
% function [ parsed_osm1] = split_ways( parsed_osm1,intersection_node_indices)
%SPLIT_WAYS Summary of this function goes here
%   Detailed explanation goes here
%%
nodes = parsed_osm1.node;
node_ids = nodes.id;

ways = parsed_osm1.way;
node_ways_ids = ways.id;
waynd = ways.nd;
waytag= ways.tag;
%% Id To Location number
intersections_id=zeros(1,38);
for j=1:size(intersection_node_indices,2)
            intersections_id(j)= node_ids(1,intersection_node_indices(j));
end
%% Location number to Waynd
for i=1:size(waynd,2)%loop for each way node struct
    counter=0;
    for k=1:size(waynd{1,i},2) %loop for each node struct element
        for j=1:size(intersections_id,2)
            if (waynd{1,i}(1,k)==intersections_id(j))
               waynd{1,i}(2,k)=intersection_node_indices(j);
               counter=counter+1;
            end
        end
    end
    if (counter == 0)
        waynd{1,i}(2,1:size(waynd{1,i},2)) = 0;
    end
end

%% Try 2
counting=zeros(1,28);
for i=1:size(waynd,2)
   mikum=1;
   for k=1:size(waynd{1,i},2)
       if (waynd{1,i}(2,k)~=0)
          counting(mikum,i)=k;
          mikum=mikum+1;
      end
   end   
end    
%% Split function -New
Waylength=size(waynd,2);
for i=1:size(waynd,2)
    num1=0;
    num2=0;
    mikum=2;
    delflag=0;
    mikumflag=0;
    if(sum(counting(:,i)~=0)-1 >1)
        for k=1:(sum(counting(:,i)~=0)-2)
            num1=counting(mikum,i);
            num2=counting(mikum+1,i);
            waynd{1,Waylength+1}(1,1:num2-num1+1)=waynd{1,i}(1,num1:num2);
            waynd{1,Waylength+1}(2,1:num2-num1+1)=waynd{1,i}(2,num1:num2);
            if(k==1  && delflag==0)
                mikumflag = num1;
                delflag=1;
            end;
            waytag(Waylength+1)=waytag(i);
            node_ways_ids(Waylength+1)= randi([1e3 9e3],1,1);
            
            Waylength=Waylength+1;
            mikum=mikum+1;
        end
    elseif (sum(counting(:,i)~=0)-1 == 1)
        delflag=1;
        mikumflag=counting(2,i);
    end
    if (delflag==1)
        waynd{1,i}(1,mikumflag+1:size(waynd{1,i},2))=0;
    end
end

%%


parsed_osm1.way.tag=waytag;
parsed_osm1.way.id=node_ways_ids;
parsed_osm1.way.nd=waynd;