% Finite Volume method for Diffusion Equation on 1D
%-------------------------------------------------------------------------%
%-------------------- Solve 1D Heat Equation: ----------------------------%
%----         u_t = beta*u_xx      ,in Omega              ----%
%-------------------------------------------------------------------------%
% Consider on Domain Omega = (0 1), t>0
% with initial condition u0(x)
% and Dirichlet boudary condition u(0,t)=g1(t); u(1,t)=g2(t)
%-------------------------------------------------------------------------%
clc
clear all
close all
%----------------%
% Domain of x
ax = 0.0;
bx = 1.0;
%Coefficience
beta = 1/16;
%----------------%
Nx = 10 ;
%-------------------------------------------------------------------------%
m = 5;
h = zeros(m,1);
error1 = zeros(m,1);
% error2 = zeros(m,1);
% error3 = zeros(m,1);
%------------------------------------------%
t0 = 0;
T = 0.2; % number of step on time
iter = 1000;
k = (T-t0)/iter;
t = t0:k:T;
for q = 1:m
    h(q) = (bx-ax)/Nx;
    % Create the mesh point
    x = zeros(Nx+1,1);
    for ii=1:Nx+1
        x(ii) = ax+(ii-1)*h(q);                           % uniform_Mesh
    end
    
    temp1 = zeros(iter+1,1);
%     temp2 = zeros(iter+1,1);
%     temp3 = zeros(iter+1,1);
    u1 = u0(x);
    u2 = u0(x);
    u3 = u0(x);
    for jj=1:iter
        t_k = jj*k;
        u1 = Dirichlet(beta,Nx,k,x,u1);                     %Forward Euler
        u1(1) = uex(0,t_k);
        u1(Nx+1) = uex(1,t_k);
%         u2 = Dirichlet(alpha,beta,1/2,Nx,k,x,x_cp,u2);      %Crank-Nilcoson
%         u3 = Dirichlet(alpha,beta,0,Nx,k,x,x_cp,u3);        %Backward Euler
        for i=1:Nx+1
             temp1(jj) = temp1(jj)+abs(u1(i)-uex(t_k+k,x(i)))^2*h(q)*k;
%             temp2(jj) = temp2(jj)+abs(u2(i)-uex(t_k,x_cp(i),alpha))^2*h(q);
%             temp3(jj) = temp3(jj)+abs(u3(i)-uex(t_k,x_cp(i),alpha))^2*h(q);
        end
    end
    error1(q) = sqrt((temp1(iter)));
%     error2(q) = sqrt((temp2(iter)));
%     error3(q) = sqrt((temp3(iter)));
    Nx = Nx*2;
end
plot(-log(h),-2*log(h),-log(h),-log(h),-log(h),-log(error1))
legend('Oder 2','Oder 1','Error Forward Euler Method')
title('Order of Error (consider spaceline) for each method')
%------------------------------------------%