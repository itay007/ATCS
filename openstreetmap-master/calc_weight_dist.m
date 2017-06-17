function [dg] = calc_weight_dist(dg,n,n1,route,weight)
for i=1:(n-n1)
    for j=1: (size(route{1,i}, 2)-1)
        first = route{1,i}(1,j);
        last = route{1,i}(1,j+1);
        dg(first,last) = dg(first,last)*weight;
        if(dg(first,last) >1)
            dg(first,last)=0.99;
        end
    end
end
