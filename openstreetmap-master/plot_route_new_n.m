function [] = plot_route_new_n(routeXY,n)
% routeXY{1,i}(1,:)  x cords
% routeXY{1,i}(2,:)  y cords
for i=1:n
     n1(i)=size(routeXY{1,i}, 2)-1;
end   
    t= max(n1) 
     
%t=size(routeXY, 2)-1;
t1=ones(1,n);

for i=1:t
   % t1=ones(1,n);
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
                h(j) = plot(x{1,j}(z),y{1,j}(z),'*','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
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