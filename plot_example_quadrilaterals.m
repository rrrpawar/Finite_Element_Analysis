%Example mesh of quadrilateral elements
% nnodes = 6
% nelements = 2

ien = [1 3  4 2; 4 3 5 6]   %index of element nodes

coordinates = [1.0 1.0; 1.1 1.5; 1.7 1.0; 1.5 1.6; 2.3 1.3;  2.0 1.8]  %coordinates of all nodes

 plot_mesh(ien,coordinates)

u = [1 2 3 4 5 6]'  %example nodal values to plot contour

plot_contour(ien,coordinates,u)
axis equal