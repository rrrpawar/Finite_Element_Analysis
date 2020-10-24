%% Gaussian Quad method
close all
E = 30*10^6;   % Young's modulus
L = 10;        % Length
Ao = 2;        % Area at fixed end
AL = 1;        % Area at free end

for nG=[1:4,8]         % Condition for seelcting the number of Gaussian points
[Xg,Wg] = GLTable(nG); % To select points and weights from the GLTable
e=[0,L];               % Interval
Je=(xe(2)-xe(1))/2;    % Jacobian at the interval

%% To compute the integral using Gaussian-quadrature method
I=0;
for i=1:nG
    X = 0.5*(1-Xg(i))*xe(1) + 0.5*(1+Xg(i))*xe(2);
    A(i)=Ao+(((AL-Ao)/L)*X);
    Eps(i)=(1/E)*(1000*(log(1+X)-log(1+L)-(X-L)))/A(i);
    I = I + Eps(i)*Je*Wg(i);
    e = abs(I- 8.7271e-04);
end
fprintf(['\n Finite element by ' num2str(nG) '-point Gaussian quadrature method, u(L) : %g \n'], I);
fprintf(['\n Difference between ' num2str(nG) '-point and 8 point , e(L) : %g \n'], e);
end