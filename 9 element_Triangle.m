%% 9 element-triangular domain %%
close all
clear all
n = 10;         % number of nodes
ne = 9;
l = 1;
j = l/3;
k = 1;
ab_co = [0 0; 0.333 0; 0.666 0; 1 0; 0.833 0.289; 0.666 0.577; 
         0.5 0.866; 0.333 0.577; 0.166 0.289; 0.5 0.289];
ien = [1 2 9; 2 3 10; 3 4 5; 2 10 9; 3 5 10;
       9 10 8; 10 5 6; 10 6 8; 8 6 7];
ien_t = [4 5; 5 6; 6 7; 7 8; 8 9; 9 1];
u_essential = [1 2 3 4];
u_nat = [5 6 7 8 9 10];
% To find K Omega
 K = zeros(n,n);
for e = 1 : ne
    node = ien(e,:);
    coord = ab_co(node,:);
    A = [ones(3,1), coord];
    Ainv = inv (A);
    B = Ainv(2:3,:);
    Area = det(A)/2 ;
    ke = k * Area * B' * B;
    K(node,node) = K(node,node) + ke;
end
% To find K Tau
K_t=zeros(n,n);
for i=1:n/3
    node= ien_t(i,:);
    ke_t=(2*j/6)*[2 1; 1 2];
    K_t(node,node)=K_t(node,node)+ke_t;
    node= ien_t(length(ien_t)/2 +i,:);
    ke_t=(3*j/6)*[2 1; 1 2];
    K_t(node,node)=K_t(node,node)+ke_t;
end
Ke = K + K_t; 
% To find F Tau
F=zeros(n,1);
for i=1:n/3
    node=ien_t(i,:);
    fe_t=(2*j/2)*[1; 1];
    F(node,1)=F(node,1)+fe_t;
    node=ien_t(length(ien_t)/2+i,:);
    fe_t=(3*j/2)*[1; 1];
    F(node,1)=F(node,1)+fe_t;
end
d(u_essential) = 2;
d(5:10) = 0;
D(u_nat) = Ke(u_nat,u_nat)\(F(u_nat) - Ke(u_nat,u_essential)*d(u_essential)');
f_D = d' + D';
plot_mesh(ien,ab_co);
% To find the displacement (Post-Processing)
del_D = zeros(ne,2);
for i = 1:ne
    node = ien(i,:);
    coord = ab_co(node,:);
    A = [ones(3,1),coord];
    Ainv = inv(A);
    B = Ainv(2:3,:);
    del_D(i,:) = (B * f_D(node))';
    Del_D(i,1) = sqrt(del_D(i,1)^2 + del_D(i,2)^2);
end
plot_contour(ien,ab_co,f_D);
disp(['temperature at the node = ',num2str(D(1,7))]);
fprintf(['gradient vector temperature at the centroid =\n '])
disp(del_D);
disp(['Magniturde of gradient vector temperature at the centroid = '])
disp(Del_D);