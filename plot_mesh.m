function plot_mesh(IEN,n_c)

nelem = size(IEN,1);
nnode = size(n_c,1);

shift = 0;  %0.05;
figure
for k=1:nelem
    node = IEN(k,:);
    coord_elem = n_c(node,:);
    x = coord_elem(:,1);
    y = coord_elem(:,2);
    line([x;x(1)],[y;y(1)],'color','black','linewidth',0.5)
    text(mean(x)-shift,mean(y),num2str(k),...
        'color','black','fontsize',11,...
        'BackgroundColor',[0.7 0.9 0.7],...
        'HorizontalAlignment','center',... 
	    'Margin',1)
    hold on 
end
for i=1:nnode
    text(n_c(i,1),n_c(i,2),num2str(i),...
        'color','red','fontsize',10)
    hold on
end
xlabel('x (m)') 
ylabel('y (m)')
axis equal
hold off
