function [] = plot_route_new_n(routeXY,n)
% 2017.6.23 (c) Ofer Keren, ofer293@gmail.com; Itay Levitan, itay007@gmail.com


for i=1:n
     n1(i)=size(routeXY{1,i}, 2)-1;
end   
t= max(n1) 
t1=ones(1,n);

for i=1:t
    for j=1:n
        if (t1(j) <= size(routeXY{1,j}, 2)-1)
            m_nd_coor{1,j}(1,i) = ((routeXY{1,j}(2,i+1)-routeXY{1,j}(2,i))/(routeXY{1,j}(1,i+1)-routeXY{1,j}(1,i)));
            b(j) = -m_nd_coor{1,j}(1,i)*routeXY{1,j}(1,i)+routeXY{1,j}(2,i);
            x{1,j}(1,:) = linspace(routeXY{1,j}(1,i),routeXY{1,j}(1,i+1),5);
            y{1,j}(1,:)=m_nd_coor{1,j}(1,i)*x{1,j}+b(j);
        end
    end
    for z=1:5
        for j=1:n
             if (t1(j) <= size(routeXY{1,j}, 2)-1) 
                 if(j==1)
                    h(j) = plot(x{1,j}(z),y{1,j}(z),'*','LineWidth',2,'MarkerEdgeColor','red','MarkerFaceColor','g','MarkerSize',10);
                 elseif(j==2)
                     h(j) = plot(x{1,j}(z),y{1,j}(z),'*','LineWidth',2,'MarkerEdgeColor','blue','MarkerFaceColor','g','MarkerSize',10);
                 elseif(j==3)
                     h(j) = plot(x{1,j}(z),y{1,j}(z),'*','LineWidth',2,'MarkerEdgeColor','green','MarkerFaceColor','g','MarkerSize',10);
                 elseif(j==4)
                     h(j) = plot(x{1,j}(z),y{1,j}(z),'*','LineWidth',2,'MarkerEdgeColor','Magenta','MarkerFaceColor','g','MarkerSize',10);
                 end
             end
        end
        for j=1:n
            if (t1(j) <= size(routeXY{1,j}, 2)-1)
                 pause(0.3);
                 set(h(j),'Visible','off')
                 
            end
        end
    end
    for j=1:n
        if (t1(j) <= size(routeXY{1,j}, 2)-1)
            t1(j)=t1(j)+1;
        end
    end
end