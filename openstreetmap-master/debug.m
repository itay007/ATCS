n=size(parsed_osm.way.nd);
n1=n(1:2);
s1=0;
s2=0;
for i=1:n1
s=size(parsed_osm.way.nd{n(i)});
s1=s(1:2);
s2=s2+s1;
end