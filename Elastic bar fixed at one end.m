%% FEA of an elastic bar fixed at one end vertically
clear all
close all
Ne=4;                       % number of elements
Nn=Ne+1;                    % number of nodes
d1=ones(1,Ne/2)*0.25/100;   % brass dian1eter
d2=ones(1,Ne/2)*0.5/100;    % aluminum diameter
d=[d1 d2];                  % diameter
L=ones(1,Ne)*10/(100*Ne);   % element length
rhol=ones(1,Ne/2)*8470;        
rho2=ones(1,Ne/2)*2700;
r=[rhol rho2];              % density
E1=ones(1,Ne/2)*97e9;
E2=ones(1,Ne/2)*69e9;
E=[E1 E2];                  % Young's modulus
A=(d.^2)*pi/4;              % Area matrix of Al & Br
g=9.81;                     % gravity

% Calculate Body force %
for e=1:Ne
    b(e)=-r(e)*g*A(e);          
end

% Index of Element Nodes %
for e=1:Ne
    IEN(e,1)=e;
    IEN(e,2)=e+1;
end

% Stiffness matrix %
K=zeros(Nn,Nn);
for e=1:Ne
    k=E(e)*A(e)/L(e);
    ke=[k -k;-k k];
    for i=1:2
        i_node=IEN(e,i);
        for j=1: 2
            j_node=IEN(e,j);
            K(i_node,j_node)=K(i_node,j_node)+ke(i,j);
        end
    end
end

% Force matrix %
for e=1:Ne
    F(1,e)=b(e)*L(e)/2;
end
for e=2:Ne
    F(1,e)=b(e-1)*L(e-1)/2+b(e)*L(e)/2;
end
F(1,Nn)=b(Ne)*L(Ne)/2;

% Reduced matrix %
F_FF=F([1:Ne]);           % Reduced force matrix
Disp=zeros(1,Ne);
K_FF=K( [1:Ne], [1:Ne]);  % Reduced stiffness matrix
D=K_FF\F_FF';             % displacement
D(Nn,1)=0;

% Stress matrix %
stress=zeros(Ne,1);
for e=1:Ne
    stress(e)=E(e)*(D(e+1)-D(e))/L(e); 
end
FL=0;
for e=1:Nn
    FL=FL+K(Nn,e)*D(e);
end

% Plot %
figure(1)
plot([0:sum(L)/Ne:sum(L)],D,'-*b','LineWidth',2)
xlabel('Length(m)','Interpreter','LaTex','fontsize',10)
ylabel('Displacement(m)','Interpreter','LaTex','fontsize',10)
grid on
set(gca,'GridLineStyle',': ','LineWidth',1.5)
title(['Displacement vs Posit.on for' ' num2str(N) '' elements'])

figure(2)
stairs([0:sum(L)/Ne:sum(L)],[stress(:,1);stress(Ne,1)],'LineWidth',2)
xlabel('Length (m)', 'Interpreter', 'LaTex', 'fontsize',10)
ylabel('Stress (N/mA2)', 'Interpreter', 'LaTex', 'fontsize',10)
grid on
set(gca,'GridLineStyle', ': ', 'LineWidth',1.5)
title(['Stress vs Positon for ' num2str(Ne) ' elements'])