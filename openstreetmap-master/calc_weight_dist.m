function [dg] = calc_weight_dist(dg,n,n1,route,weight,i0)
% 2017.6.23 (c) Ofer Keren, ofer293@gmail.com; Itay Levitan, itay007@gmail.com

for i=i0:(n-n1)
    for j=1: (size(route{1,i}, 2)-1)
        first = route{1,i}(1,j);
        last = route{1,i}(1,j+1);
        dg(first,last) = dg(first,last)*weight;
        if(dg(first,last) >1)
            dg(first,last)=0.99;
        end
    end
end
