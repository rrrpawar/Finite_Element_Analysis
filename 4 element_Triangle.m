%% 4 element-triangular domain %%
close all
clear all
syms u
n = 6;      % number of nodes
ne = 4;
u_inf = 1;
h2 = 2;
h3 = 3;
l = 1;
j = l/2;
k = 1;
n_c = [0 0; 0.5 0; 0.25 0.433; 0.75 0.433; 1 0; 0.5 0.866];
IEN = [1 2 3; 2 4 3; 2 5 4; 4 6 3];
IEN_T = [5 4; 4 6; 6 3; 3 1];
u_essential=[1 2 5];
u_nat=[4 6 3];
% To find K Omega
 K = zeros(n,n);
for e = 1 : ne
    node = IEN(e,:);
    coord = n_c(node,:);
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
    node= IEN_T(i,:);
    ke_t=(2*j/6)*[2 1; 1 2];
    K_t(node,node)=K_t(node,node)+ke_t;
    
    node= IEN_T(length(IEN_T)/2 +i,:);
    ke_t=(3*j/6)*[2 1; 1 2];
    K_t(node,node)=K_t(node,node)+ke_t;
end
Ke = K + K_t;
% To find F Tau
F=zeros(n,1);
for i=1:n/3
    node=IEN_T(i,:);
    fe_t=(2*j/2)*[1; 1];
    F(node,1)=F(node,1)+fe_t;
    node=IEN_T(length(IEN_T)/2+i,:);
    fe_t=(3*j/2)*[1; 1];
    F(node,1)=F(node,1)+fe_t;
end
d(u_essential) = 2;
d(6) = 0;
D(u_nat) = Ke(u_nat,u_nat)\(F(u_nat) - Ke(u_nat,u_essential)*d(u_essential)');
f_D = d'+ D';
plot_mesh(IEN,n_c);
plot_contour(IEN,n_c,f_D);
% To find the displacement (Post-Processing)
del_D = zeros(ne,2);
for i = 1:ne
    node = IEN(i,:);
    coord = n_c(node,:);
    A = [ones(3,1),coord];
    Ainv = inv(A);
    B = Ainv(2:3,:);
    del_D(i,:) = (B * f_D(node))';
    Del_D(i,1) = sqrt(del_D(i,1)^2 + del_D(i,2)^2); 
end
disp(['temperature at the node = ',num2str(D(1,6))]);
fprintf(['gradient vector temperature at the centroid =\n '])
disp(del_D);
disp(['Magniturde of gradient vector temperature at the centroid = '])
disp(Del_D);