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
kappa = 1/16;           %kappa

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%------------------------------------------%
%----------------- Error ------------------%
%------------------------------------------%
% Consider error on spaceline %
%--- Error on spaceline ---%
%--------------------------%
% Fixed timestep
Nt = 1200;
k = (T-t0)/Nt;
t = t0:k:T;
%------------------------------------------%
n=5;
h = zeros(n,1);
error1 = zeros(n,1);
error2 = zeros(n,1);
error3 = zeros(n,1);
%-----------------------%
Nx = 10;
%------------------------------------------%

for q = 1:n
    h(q) = (b-a)/Nx;
    x = a:h(q):b;
    temp1 = zeros(Nt+1,1);
    temp2 = zeros(Nt+1,1);
    temp3 = zeros(Nt+1,1);
    u1 = u0(x);
    u2 = u0(x);
    u3 = u0(x);
    for j=1:Nt
        u1 = FDM_for_HeatEquation(h(q),k,Nt,Nx,kappa,0,u1);     %Forward Euler
        u2 = FDM_for_HeatEquation(h(q),k,Nt,Nx,kappa,1/2,u2);   %Crank-Nilcoson
        u3 = FDM_for_HeatEquation(h(q),k,Nt,Nx,kappa,1,u3);     %Backward Euler
        for i=1:Nx+1
            temp1(j) = temp1(j)+abs(u1(i)-uex(t(j+1),x(i)))^2*h(q)*k;
            temp2(j) = temp2(j)+abs(u2(i)-uex(t(j+1),x(i)))^2*h(q)*k;
            temp3(j) = temp3(j)+abs(u3(i)-uex(t(j+1),x(i)))^2*h(q)*k;
        end
    end
    error1(q) = sqrt(sum(temp1));
    error2(q) = sqrt(sum(temp2));
    error3(q) = sqrt(sum(temp3));
    Nx = Nx*2;
end
plot(-log(h),-2*log(h),-log(h),-log(h)+2,-log(h),-log(error1),'x-',-log(h),-log(error2),'x-',-log(h),-log(error3),'x-')
legend('2x','1x','Forward','Crank-Nilcoson','Backward')
title('Order of Error (consider spaceline) for each method')
%------------------------------------------%