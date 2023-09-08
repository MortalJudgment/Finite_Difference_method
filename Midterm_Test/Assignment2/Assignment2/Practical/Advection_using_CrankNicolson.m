% Approximately for Advection equation
% Using Crank-Nicolson
clc
close all
clear all
format long
a = 0; b = 1;
T = 0.1;
%-----------------------------------%
k = 0.001;
Nt = T/k;
h = 0.0025;
Nx = (b-a)/h;
%-----------------------------------%
t = 0:k:T;
x = a:h:b;
%-----------------------------------%
speed = 2;
%-------------------------------------------------------------------------%
% type: type of method (e.g: Forward-Euler or Crank-Nicolson)
% nameFun: name of the function (e.g Advection or Diffusion)
% nameSche: name of the Schemes (e.g Lax-Friedrich, Lax-Wendroff or Upwind)
%----------------------- Lax Friedrichs ----------------------------------%
%------------- Using Forward Euler for advection equations ---------------%
%---- Checking stable conditions ----%
type = 'Crank-Nicolson';
nameFun = 'Advection';
kappa = 0.0;
r = Check_Stable(type,nameFun,'Lax-Wendroff',k,h,speed,kappa);
% 3 schemes for Advection equation have the same condition
% So check 1 for all
if r == 0
    disp('Not sastify stable condition !');
elseif r == 1
    u1 = (u0(x))';
    u2 = (u0(x))';
    u3 = (u0(x))';
    uexact = zeros(Nx+1,1);
    for i=1:Nt
        clf
        %-----------------------------------------------------------------%
        % Lax-Friedrich scheme for advection equation using Crank-Nicolson%
        
        u1 = Lax_Friedrich_advection(h,k,Nx,speed,u1,type);
        u1(1) = uex(x(1),t(i),speed,kappa,nameFun);
        u1(Nx+1) = uex(x(Nx+1),t(i),speed,kappa,nameFun);
        %-----------------------------------------------------------------%
        %  Lax-Wendroff scheme for advection equation using Crank-Nicolson%
        
        u2 = Lax_Wendroff_advection(h,k,Nx,speed,u2,type);
        u2(1) = uex(x(1),t(i),speed,kappa,nameFun);
        u2(Nx+1) = uex(x(Nx+1),t(i),speed,kappa,nameFun);
        %-----------------------------------------------------------------%
        %------ Upwind scheme for advection equation using Crank-Nicolson %
        
        u3 = Upwind_advection(h,k,Nx,speed,u3,type);
        u3(1) = uex(x(1),t(i),speed,kappa,nameFun);
        u3(Nx+1) = uex(x(Nx+1),t(i),speed,kappa,nameFun);
        %-----------------------------------------------------------------%
        %------ Exact Solution for advection equation --------------------%
        for j=1:Nx+1
            uexact(j) = uex(x(j),t(i),speed,kappa,nameFun);
        end
        hold on
        plot(x,uexact,' -r','LineWidth',2,'MarkerSize',3)
        plot(x,u1,' -g','LineWidth',2,'MarkerSize',3)
        plot(x,u2,' -b','LineWidth',2,'MarkerSize',3)
        plot(x,u3,' -c','LineWidth',2,'MarkerSize',3)
        legend('Exact-solution','Lax-Friedrich','Lax-Wendroff','Upwind')
        axis([0 1 -1 1])
%         drawnow
        pause(10^(-6))
    end
   hold off 
else
    disp('Not supported yet !');
end



