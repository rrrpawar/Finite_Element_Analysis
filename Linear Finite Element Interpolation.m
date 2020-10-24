%% Linear Finite Element Interpolation

close all
N=2;                     % Number of elements
L=10;                    % Length
E=30*10^6;               % Young's modulus
A=1.5;                   % Area
x=0:L/N:L;               % Defining nodes
arb_x=0:0.5:L;           % Arbitrary position
 
%% To compute displacement at nodes on the bar
for j=1:N+1
Ux(j) = (1000/E*A)*((L*x(j))-((x(j)^2)/2)+((log(x(j)+1)-1)*(x(j)+1))-(x(j)*log(L+1))+1);
end
 
%% To compute exact displacement at nodes
for j=1:(L/0.5)+1
ex_u(j) = (1000/E*A)*((L*arb_x(j))-((arb_x(j)^2)/2)+((log(arb_x(j)+1)-1)*(arb_x(j)+1))-(arb_x(j)*log(L+1))+1);
end
 
%% Initializing the variables
u=zeros(1,(L/0.5)+1);   % Initialising the displacement matrix 
le =L/N;                % Element length
Ap_Strain =[0];         % Initialising approx strain
Ex_Strain =[0];         % Initialising exact strain
 
%% To compute linear finite element interpolation
for j = 1:(L/0.5)+1
    for k = 1:N
        if arb_x(j)>x(k) && arb_x(j)<=x(k+1)
            n1=((arb_x(j)-x(k+1))/(x(k)-x(k+1)));
            n2=((arb_x(j)-x(k))/(x(k+1)-x(k)));
            u(1,j)=n1*Ux(k)+n2*Ux(k+1);
            Be=[-1, 1]*(1/le);
            St=(Ux(k+1)-Ux(k))/le;
            Ex_Strain =[Ex_Strain,St];
            Eps=Be*[Ux(k);Ux(k+1)];
            Ap_Strain = [Ap_Strain,Eps];
        end
    end
end
 
%% Plot position vs displacement
figure(1)
plot(arb_x,ex_u)
hold on
plot(arb_x,u)
title('Position vs Displacement')
xlabel('Position (in)');
ylabel('Displacement (in)');
legend('Exact Solution','FEM Solution')
 
%% Plot strain vs position
figure(2)
plot(arb_x,Ex_Strain)
hold on
stairs(arb_x,Ap_Strain)
title('Strain vs Position')
xlabel('Position (in)')
ylabel('Strain')
legend('Exact Solution','FEA Solution')
