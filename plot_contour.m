function  plot_contour(IEN,n_c,D)
%Input nodal vector C to plot; 
%C must be a column vector with same number of rows as coord
figure 
patch('Faces', IEN, 'Vertices', n_c, 'FaceVertexCData',D,'Facecolor','interp');
colorbar;
colormap jet
cmin = min(D);
cmax = max(D);
caxis([cmin cmax]);
end

