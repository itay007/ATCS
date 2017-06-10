function [dg] = calc_weight_dist(dg,n,n1,route,weight)
for i=1:(n-n1)
%     for i=1:size(counting,2)
%         if(counting(2,i)~=0)
%             first = counting(1,i);
%             last = counting(2,i);
%              dists(2,i) = waynd{1,i}(2,first);
%              dists(3,i) = waynd{1,i}(2,last);
% 
%         end
%     end
    for j=1: (size(route{1,i}, 2)-1)
        first = route{1,i}(1,j);
        last = route{1,i}(1,j+1);
        dg(first,last) = dg(first,last)*weight;
        if(dg(first,last) >1)
            dg(first,last)=0.99;
        end
    end
end
