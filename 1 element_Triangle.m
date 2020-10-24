%% 1 element-triangular domain %%
close all
clear all
n = 3;       % number of nodes
ne = 1;
l = 1;
j = l/2;
k = 1;
ab_co = [0 0; 1 0; 0.5 0.866];
ien = [1 2 3];
ien_t = [2 3; 3 1];
u_essential=[1 2];
u_nat=[3];
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
    ke_t=(2*l/6)*[2 1; 1 2];
    K_t(node,node)=K_t(node,node)+ke_t;
    node= ien_t(length(ien_t)/2 +i,:);
    ke_t=(3*l/6)*[2 1; 1 2];
    K_t(node,node)=K_t(node,node)+ke_t;
end
Ke = K + K_t;
% To find F Tau
F=zeros(n,1);
for i=1:n/3
    node=ien_t(i,:);
    fe_t=(2*l/2)*[1; 1];
    F(node,1)=F(node,1)+fe_t;
    node=ien_t(length(ien_t)/2+i,:);
    fe_t=(3*l/2)*[1; 1];
    F(node,1)=F(node,1)+fe_t;
end
d(u_essential) = 2;
d(3) = 0;
D(u_nat) = Ke(u_nat,u_nat)\(F(u_nat) - Ke(u_nat,u_essential)*d(u_essential)');
f_D = d'+ D';
plot_mesh(ien,ab_co);
plot_contour(ien,ab_co,f_D);
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
disp(['temperature at the node = ',num2str(D(1,3))]);
disp(['gradient vector temperature at the centroid = ',num2str(del_D)]);
disp(['Magniturde of gradient vector temperature at the centroid = ',num2str(Del_D)]);