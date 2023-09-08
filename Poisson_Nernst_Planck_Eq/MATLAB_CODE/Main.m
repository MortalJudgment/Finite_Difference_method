 clc
close all
clear all
format long
%-------------------------------------------------------------------------%
%----------- Coefficient ----------%
d1  = 1;    d2  = d1;
v1  = 0.3;  v2  = -4;
y1  = 1/2;  y2  = -1;
g11 = 7;    g12 = 8;    g21 = 9;    g22 = 33;
%---------------------------%
a = -1; b = 1;               %Interval of space
t0 = 0; T = 0.1;             %Interval of time
%---------------------------%
dt = 0.000001;
Nt = (T-t0)/dt;
Nx = 100;
dx = (b-a)/Nx;
t = t0:dt:T;
x = a:dx:b;
%---------------------------%
u1 = zeros(Nx+1,1);
u2 = zeros(Nx+1,1);
for i = 1:Nx+1
    u1(i) = u0(x(i));
    u2(i) = v0(x(i));
end
phi1 = Poisson_Eq(Nx,dx,y1,y2,u1,u2);
% figh = figure;
for j=2:Nt+1
    clf                             
    u11 = u1;
    u22 = u2;
    t_j = t(j);
    u1 = Nernst_Planck(Nx,dx,dt,d1,v1,g11,g12,u11,u22,phi1);
    u2 = Nernst_Planck(Nx,dx,dt,d2,v2,g22,g21,u22,u11,phi1);
    phi1 = Poisson_Eq(Nx,dx,y1,y2,u1,u2);
    
    %----- Draw all of them at t = t_k -----%
    plot(x,u1,x,u2,x,phi1)
    
    xlabel('x')
    title(['Particle at t = ', num2str(t_j*10^3),'.e-03 seconds'])
%     axis auto
    
    pause(0.005)
%     movieVector(j) = getframe(figh, [10 10 520 400]);
end
% myWriter = VideoWriter('Poisson-Nernst-Planck');
% myWriter.FrameRate = 20;
% open(myWriter);
% writeVideo(myWriter,movieVector);
% close(myWriter);