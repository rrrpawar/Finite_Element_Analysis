
%Example mesh of triangle elements
% nnodes = 4
% nelements = 2

ien = [1 3 2; 3 4 2]   %index of element nodes

coordinates = [1.0 1.0; 1.25 1.5; 1.8 1.0; 1.9 1.7]  %coordinates of all nodes

 plot_mesh(ien,coordinates)

u = [1 2 3 4]'  %nodal values to plot contour

plot_contour(ien,coordinates,u)
axis equal