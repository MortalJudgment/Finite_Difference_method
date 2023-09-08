clc
close all
clear all
format long
%-------------------------------------------------------------------------%
% Solve PDE: u_t = kappa*u_xx
% For all 0<=x<=1
%---------------------------%
a = 0; b = 1;               %Interval of space
t0 = 0; T = 1;              %Interval of time
%---------------------------%
Nt = 100;
k = (T-t0)/Nt;
Nx = 20;
h = (b-a)/Nx;
t = t0:k:T;
x = a:h:b;
%---------------------------%
kappa = 1/16;               % Kappa
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%---------------------------%
%------ Forward Euler ------%
u1 = u0(x);
for i=1:Nt
    u1 = FDM_for_HeatEquation(h,k,Nt,Nx,kappa,0,u1);
end
plot(x,u1,'-. g','LineWidth',2,'MarkerSize',3)
%---------------------------%
hold on
%---------------------------%
%------ Crank-Nilcoson -----%
u2 = u0(x);
for i=1:Nt
    u2 = FDM_for_HeatEquation(h,k,Nt,Nx,kappa,1/2,u2);
end
plot(x,u2,'-. b','LineWidth',2,'MarkerSize',3)
%---------------------------%
%------ Backward Euler -----%
u3 = u0(x);
for i=1:Nt
    u3 = FDM_for_HeatEquation(h,k,Nt,Nx,kappa,1,u3);
end
plot(x,u3,'-. y','LineWidth',2,'MarkerSize',3)
%---------------------------%
%------ Exact Solution -----%
uexact = zeros(Nx+1,1);
for i=1:Nx+1
    uexact(i) = uex(T,x(i));
end
%------------------------------------------%
%---------------- Drawing -----------------%
%------------------------------------------%
plot(x,uexact,'-r','LineWidth',2,'MarkerSize',3)
legend('Forward Euler','Crank-Nicolson','Backward Euler','Exact Solution')
hold off