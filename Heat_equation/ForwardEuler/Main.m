% Finite Volume method for Diffusion Equation on 1D
% Consider on Domain Omega = (0 1), t>0
% with initial condition u0(x)
% and Dirichlet boudary condition u(0,t)=g1(t); u(1,t)=g2(t)
%-------------------------------------------------------------------------%
%-------------------- Solve 1D Heat Equation: ----------------------------%
%----         u_t = beta*u_xx      ,in [a,b]              ----%
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
Nx = 50 ;
%-------------------------------------------------------------------------%

dx = (bx-ax)/Nx;
dt = 1/2*1/(2*beta/dx^2);
% Create the mesh point
x = zeros(Nx+1,1);
for ii=1:Nx+1
    x(ii) = ax+(ii-1)*dx;                           % uniform_Mesh
    %     x(ii) = bx-cos(pi/2*(i_iter-1)/N);        % non-uniform_Mesh
end
%-----------------%
%Initial condition
u1 = (u0(x));
% theta = 1/2;
%------------------------------------------------------------------------%
for jj=1:300
    clf
    t_k = jj*dt;
    %---------- Dicrete solution -----------%
    u_dis = Dirichlet(beta,Nx,dt,x,u1);
    u_dis(1) = uex(0,t_k);
    u_dis(Nx+1) = uex(1,t_k);
    
    u1 = u_dis;
    %----------- Exact solution ------------%
    u_ex = zeros(Nx+1,1);
    for ii=1:Nx+1
        u_ex(ii) = uex(x(ii),t_k);
    end
    %---------------------------------------%
    %-------------- Drawing ----------------%
    %---------------------------------------%
    hold on
    plot(x,u_dis,'o r','LineWidth',2,'MarkerSize',4);
    plot(x,u_ex,'b');
    
    %     legend('discrete solution','exact solution')
    axis([0 1 -1 1])
    title(['Particle at t = ', num2str(jj),'*dt seconds'])
    pause(10^(-2))
end
hold off