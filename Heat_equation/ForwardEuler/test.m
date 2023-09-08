clc
close all
clear all
format long
%-------------------------------------------------------------------------%
% Solve PDE: u_t = kappa*u_xx
% For all 0<=x<=1
%-----------------------%
a = 0; b = 1;           %Interval of space
t0 = 0; T = 1;          %Interval of time
%-----------------------%
beta = 1/16;           %kappa

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%------------------------------------------%
%----------------- Error ------------------%
%------------------------------------------%
% Consider error on spaceline %
%--- Error on spaceline ---%
%--------------------------%
% Fixed timestep
k = 0.0001;
Nt = (T - t0)/k;
%------------------------------------------%
n=5;
h = zeros(n,1);
error1 = zeros(n,1);
%-----------------------%
Nx = 10;
%------------------------------------------%

for q = 1:n
    h(q) = (b-a)/Nx;
    x = a:h(q):b;
    temp1 = 0;
    u1 = u0(x)';
    for jj=1:Nt+1
        t_k = jj*k;
        %---------- Dicrete solution -----------%
        u_dis = Dirichlet(beta,Nx,k,x,u1);
        u_dis(1) = uex(0,t_k);
        u_dis(Nx+1) = uex(1,t_k);
        
        u1 = u_dis;
        %----------- Exact solution ------------%
        u_ex = zeros(Nx+1,1);
        for ii=1:Nx+1
            u_ex(ii) = uex(x(ii),t_k);
        end
        %--------- 
        %-------------  Error  ----------------%
        for ii=1:Nx+1
            temp1 = temp1 + abs(u_ex(ii) - u_dis(ii))^2*h(q);
        end
    end
%     hold off
    error1(q) = sqrt(temp1);
    Nx = Nx*2;
end
figure
plot(-log(h),-2*log(h)-2,-log(h),-log(h),-log(h),-log(error1)+3)
legend('2x','1x','Forward Euler')
title('Order')
%------------------------------------------%