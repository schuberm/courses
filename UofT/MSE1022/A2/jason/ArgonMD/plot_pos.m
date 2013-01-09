function plot_pos3D(x,L)

plot3(x(:,1),x(:,2),x(:,3),'.')
axis([0 L(1,1) 0 L(1,2) 0 L(1,3)])
