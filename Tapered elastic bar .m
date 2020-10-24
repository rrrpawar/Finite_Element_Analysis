%% FEA of tapered elastic bar with mesh of 4 and 8 elements
clear all;
close all;
Ne=8; N=Ne+1; 
TL=10; 
L=(TL/Ne); 
A_1(1)=1;
A_1(N)=2; 
% Area at nodes
for x=2:Ne
A_1(x)=A_1(x-1)+(TL/(Ne*10)); 
end
% Matrix reversal
At=flip(A_1); 
% finding average
for x=1:Ne
A(x)=((At(x)+At(x+1))/2); 
end
E=30e3*ones(1,Ne); 
% To calculate the body force
y=0;
for x=0:L:TL
    y=y+1;
    b(y)=(x/(1+x));
end
for x=1:Ne
    y=0;
    y=y+1;
    be(x,y)=((2*b(x))+b(x+1));
    be(x,y+1)=(b(x)+(2*b(x+1)));
end

bv(1,1)=(L/6)*be(1,1);

for x=2:Ne
    bv(x,1)=(L/6)*(be(x-1,2)+be(x,1));
end

bv(N,1)=(L/6)*be(Ne,2);

%Index of element nodes 
for x=1:Ne
    for i=1:2
        IEN(x,i)=x+(i-1);
    end
end

K = zeros(N,N);
for x=1:Ne
    k=E(x)*A(x)/L;
    ke=k*[1, -1; -1, 1];
    node = IEN(x,:);
    K(node,node) = K(node,node) + ke(1:2,1:2);
end

Be=bv([2:N]);        % Reduced force matrix
Ke=K([2:N],[2:N]);  % Reduced Stiffness matrix
de=Ke\Be;           % To compute displacement
d(1,1)=0;           % Displacement at the fixed end 

for x=1:Ne
d(x+1,1)=de(x,1);
end

figure(1)
plot(0:TL/Ne:TL,d)
xlabel('Position (in)');
ylabel('Displacement (in)');
grid on
title(['For ' num2str(Ne) ' elements']);

%To compute stress
for x=1:Ne
S(x)=(E(x)*(d(x+1)-d(x)))/L;
end
figure(2)
stairs([0:TL/Ne:TL-1],S)
xlabel('Position (in)');
ylabel('Stress (ksi)');
grid on
title(['For ' num2str(Ne) ' elements']);
