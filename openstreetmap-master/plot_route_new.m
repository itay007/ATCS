function [] = plot_route_new(routeXY)

%% for debug only
n=size(routeXY, 2)-1;

for i=1:n
    m_nd_coor(1,i) = ((routeXY(2,i+1)-routeXY(2,i))/(routeXY(1,i+1)-routeXY(1,i)));
    b = -m_nd_coor(1,i)*routeXY(1,i)+routeXY(2,i);
    x(1,:) = linspace(routeXY(1,i),routeXY(1,i+1),5);
    y=m_nd_coor(1,i)*x+b;
    
    for z=1:5
        h2 = plot(x(z),y(z),'*','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
        pause(0.3);
        set(h2,'Visible','off')
    end
end