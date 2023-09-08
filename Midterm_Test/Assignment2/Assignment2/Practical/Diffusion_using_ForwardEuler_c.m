% Approximately for Diffusion equation
% Using Forward-Euler
% Consider kappa = 1/16
clc
close all
clear all
format long
a = 0; b = 1;
T = 0.5;
%-----------------------------------%
kappa = 1/16;
speed = 2;
h = 0.025;
k = (h^2/(4*kappa));
Nt = T/k;
Nx = (b-a)/h;
%-----------------------------------%
t = 0:k:T;
x = a:h:b;
%-----------------------------------%
%-------------------------------------------------------------------------%
% type: type of method (e.g: Forward-Euler or Crank-Nicolson)
% nameFun: name of the function (e.g Advection or Diffusion)
% nameSche: name of the Schemes (e.g Lax-Friedrich, Lax-Wendroff or Upwind)
%----------------------- Lax Friedrichs ----------------------------------%
%------------- Using Forward Euler for advection equations ---------------%
%---- Checking stable conditions ----%
type = 'Forward-Euler';
nameFun = 'Diffusion';
r1 = Check_Stable(type,nameFun,'Lax-Friedrich',k,h,speed,kappa);
r2 = Check_Stable(type,nameFun,'Lax-Wendroff',k,h,speed,kappa);
r3 = Check_Stable(type,nameFun,'Upwind',k,h,speed,kappa);
% 3 schemes for Advection equation have the same condition
% So check 1 for all
if (r1 == 0)||(r2 == 0)||(r3 == 0)
    if r1 == 0
        disp('Not sastify Lax-Friedrich stable condition !');
    elseif r2 == 0
        disp('Not sastify Lax-Wendroff stable condition !');
    else
        disp('Not sastify Upwind stable condition !');
    end
elseif (r1 == 1)&&(r2 == 1)&&(r3 == 1)
    u1 = (u0(x))';
    u2 = (u0(x))';
    u3 = (u0(x))';
    uexact = zeros(Nx+1,1);
    for i=1:Nt
        clf
        %-----------------------------------------------------------------%
        % Lax-Friedrich scheme for advection equation using Forward Euler %
        
        u1 = Lax_Friedrich_diffusion(h,k,Nx,speed,kappa,u1,type);
        u1(1) = uex(x(1),t(i),speed,kappa,nameFun);
        u1(Nx+1) = uex(x(Nx+1),t(i),speed,kappa,nameFun);
        %-----------------------------------------------------------------%
        %  Lax-Wendroff scheme for advection equation using Forward Euler %
        
        u2 = Lax_Wendroff_diffusion(h,k,Nx,speed,kappa,u2,type);
        u2(1) = uex(x(1),t(i),speed,kappa,nameFun);
        u2(Nx+1) = uex(x(Nx+1),t(i),speed,kappa,nameFun);
        %-----------------------------------------------------------------%
        %------ Upwind scheme for advection equation using Forward Euler  %
        
        u3 = Upwind_diffusion(h,k,Nx,speed,kappa,u3,type);
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
        axis([0 1 -1 1])
%         drawnow
        pause(10^(-6))
    end
   hold off 
else
    disp('Not supported yet !');
end



