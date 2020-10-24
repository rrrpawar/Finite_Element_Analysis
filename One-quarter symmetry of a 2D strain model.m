%% Problem: One-quarter symmetry of two-dimensional strain model
% Created by Raghavendra Rao
%% Main program
clear all 
close all
E = 30e6;               % Modulus of ELasticity
nu = 0.3;               % Poisson's ratio
p = 5000;               % Internal pressure
D = (E/((1+nu)*(1-2*nu)))*[1-nu nu 0; nu 1-nu 0; 0 0 0.5*(1-2*nu)];
IEN = [1 2 5 4; 2 3 6 5; 4 5 8 7; 5 6 9 8; 7 8 11 10; 8 9 12 11;10 11 14 13; 11 12 15 14];
IEN_T = [1 4; 4 7; 7 10; 10 13]; % node indexing on elements
nodes = 15;             % Number of nodes
n_dof = 2;              % Degrees of freedom
N = nodes*n_dof;        % Defining total dof  
n_ele = 8;              % Number of elements
R = [60,65,70];         % Radius
alpha = 22.5;           % Angle to define polar coordinates
d_alpha = [0,alpha,2*alpha,3*alpha,4*alpha];
co_ord = zeros(nodes,2);% Defining coordinates
i = 1;                  % Defining counter
for j = 1:N/6
    for k = 1:3
         co_ord(i,1) = R(k)*cosd(d_alpha(j));
         co_ord(i,2) = R(k)*sind(d_alpha(j)); 
         i = i + 1;
    end
end

% To define Stiffnes Matrix
x_co = [-1; 1]/sqrt(3);  
weight = [1;1];
K = zeros(N,N);                   % To define element matrix
for e = 1 : n_ele
    K_ele = zeros(n_ele,n_ele);   % Re-define element matrix
    n = IEN(e,:);                 % Nodes for the element
    co_elem = co_ord(n,:);
    lm = [2*n(1)-1, 2*n(1), 2*n(2)-1, 2*n(2), 2*n(3)-1, 2*n(3), 2*n(4)-1, 2*n(4)];
    for i = 1:2
        for j = 1:2
            psi = x_co(j); eta = x_co(i);
            B_hat = 0.25*[eta-1, 1-eta, 1+eta, -eta-1; psi-1, -psi-1, 1+psi, 1-psi];
            JT = B_hat * co_elem;
            J = det(JT);
            G = JT \ B_hat;       % inv(JT)*Bhat
            B = [];
            for c = 1:4
                B_node = [G(1,c) 0; 0 G(2,c); G(2,c) G(1,c)];
                B = [B , B_node]; % Stack nodal partitions by columns
            end
            K_ele = K_ele + B'*D*B*J*weight(i)*weight(j);
        end
    end    
    K(lm,lm) = K(lm,lm) + K_ele;  % Equivalent K matrix
end

% To define force matrix
edge_ele = 4;
edge_nodes = 5;                   
ndof_node = 2;
n_eq = edge_nodes * ndof_node; 
F = zeros(N,1);
for k = 1:edge_ele
    n = IEN_T(k,:);               % nodes for the boundary element
    coord_edge = co_ord(n,:);     % to extract coordinates for boundary element
    x = coord_edge(:,1);
    y = coord_edge(:,2);
    edge_length = sqrt( (x(2)-x(1))^2 + (y(2)-y(1))^2 );
    pos_vector = [x(2)-x(1); y(2)-y(1)];
    unit_tangent = pos_vector / edge_length;
    unit_normal = [unit_tangent(2); -unit_tangent(1)]; % using inward pointing normal
    traction_vector = p * unit_normal;  % Workequivalent nodal force components using analytical integration 
    f_node = (0.5*edge_length) * (traction_vector);
    f_gamma = [f_node; f_node];
    lm = [2*n(1)-1, 2*n(1), 2*n(2)-1, 2*n(2)];% vector components odd_even
    F(lm,1) = F(lm,1) + f_gamma;
end

% To define displacement
u = zeros(N,1);   
n_l = ([1,3,5,7: 24,26,28,30]);    % Node listing
K_K = K(n_l,n_l);    
F_F = F(n_l);     
u(n_l) = K_K\F_F;     

% To plot the deformed element shape
u_node = reshape(u,[2,15])';
new_point = co_ord + 20*u_node;    % scaling the deformation
plot_mesh(IEN,co_ord,new_point);   % Call mesh plot function
grid on
title(['2D plane strain model for infinite cylinder with internal pressure'])
xlabel('x') 
ylabel('y')
axis equal

% To define the location of largest magnitude displacement vector
for i=1:15
    u_n(i) = sqrt(u_node(i,1)^2+u_node(i,2)^2); 
end
u_max = max(u_n);                            % largest magnitude
for l=1:15
    if u_n(l) == u_max
        coord_u_max = co_ord(l,:);
    end
end

% Post processing to find principal stresses

% To find principal normal Stress at the center of the element
St_c = [];
for e = 1 : n_ele
    n = IEN(e,:);              % nodes for the element
    co_elem = co_ord(n,:);
    u_e = u_node(n,:);
    u_elem = reshape((u_e)',[1,8])';
    lm = [2*n(1)-1, 2*n(1), 2*n(2)-1, 2*n(2), 2*n(3)-1, 2*n(3), 2*n(4)-1, 2*n(4)];
    psi_c = 0; eta_c = 0;
    Bhat_c = 0.25*[eta_c-1, 1-eta_c, 1+eta_c, -eta_c-1; psi_c-1, -psi_c-1, 1+psi_c, 1-psi_c];
    JT_c = Bhat_c * co_elem;
    J_c = det(JT_c);
    G_c = JT_c \ Bhat_c;       % inv(JT)*B_hat
    B_c = [];
    for c = 1:4
        Bnode_c = [G_c(1,c) 0; 0 G_c(2,c); G_c(2,c) G_c(1,c)];
        B_c = [B_c , Bnode_c];  % Stack nodal partitions by columns
    end
    St_c = [St_c, D*B_c * u_elem];
end

% To find maximum and minimum principal normal stresses
for i=1:8
    St_avg = (St_c(1,i)+St_c(2,i))/2;
    R = sqrt(((St_c(1,i)-St_avg)^2)+St_c(3,i)^2);
    St_p(1,i) = St_avg + R;
    St_p(2,i) = St_avg - R;
end
% To deifne shape function at the centroid
St_pmax = max(St_p(1,:));
St_pmin = min(St_p(2,:));
N_c = 0.25*ones(1,4);
coord_c = zeros(n_ele,2); 
for e=1:n_ele
    n = IEN(e,:);
    lm = [2*n(1)-1, 2*n(1), 2*n(2)-1, 2*n(2), 2*n(3)-1, 2*n(3), 2*n(4)-1, 2*n(4)];
    co_elem = co_ord(n,:);
    for j=1:2
        for i = 1:4
         coord_c(e,j) = coord_c(e,j) + N_c(i)*[co_elem(i,j)];
        end
    end
end

for l=1:8
    if St_p(1,l) == St_pmax
        coord_St_max = coord_c(l,:);
    end
    if St_p(2,l) == St_pmin
        coord_St_min = coord_c(l,:);
    end
end

% Tabulating the results
disp('The displacement at the nodes'); disp(u_node);
disp(['Magnitude of the highest displacement vector = ', num2str(u_max),' in']);
disp(['location ','(',num2str(coord_u_max(1)),',',num2str(coord_u_max(2)),')']);
fprintf('\n');
disp('The value of principle normal stresses at the center');
disp(St_c); fprintf('\n');
disp(['The value of the minimum principal stress = ', num2str(St_pmin),' psi']);
disp(['location ','(',num2str(coord_St_min(1)),',',num2str(coord_St_min(2)),')']);
disp(['The value of the maximum principal stress = ', num2str(St_pmax),' psi']);
disp(['location ','(',num2str(coord_St_max(1)),',',num2str(coord_St_max(2)),')']);
