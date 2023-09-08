% Approximately for Advection equation
% Using Forward Euler
clc
close all
clear all
format long
a = 0; b = 1;
T = 0.5;
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
type = 'Forward-Euler';
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
    % figh = figure;
    for i=1:Nt
        clf
        %-----------------------------------------------------------------%
        % Lax-Friedrich scheme for advection equation using Forward Euler %
        
        u1 = Lax_Friedrich_advection(h,k,Nx,speed,u1,type);
        u1(1) = uex(x(1),t(i),speed,kappa,nameFun);
        u1(Nx+1) = uex(x(Nx+1),t(i),speed,kappa,nameFun);
        %-----------------------------------------------------------------%
        %  Lax-Wendroff scheme for advection equation using Forward Euler %
        
        u2 = Lax_Wendroff_advection(h,k,Nx,speed,u2,type);
        u2(1) = uex(x(1),t(i),speed,kappa,nameFun);
        u2(Nx+1) = uex(x(Nx+1),t(i),speed,kappa,nameFun);
        %-----------------------------------------------------------------%
        %------ Upwind scheme for advection equation using Forward Euler  %
        
        u3 = Upwind_advection(h,k,Nx,speed,u3,type);
        u3(1) = uex(x(1),t(i),speed,kappa,nameFun);
        u3(Nx+1) = uex(x(Nx+1),t(i),speed,kappa,nameFun);
        %-----------------------------------------------------------------%
        %------ Exact Solution for advection equation --------------------%
        for j=1:Nx+1
            uexact(j) = uex(x(j),t(i),speed,kappa,nameFun);
        end
        %---------------------------- Drawing ----------------------------%
        hold on
        plot(x,uexact,' -r','LineWidth',2,'MarkerSize',3)
        plot(x,u1,' xg','LineWidth',2,'MarkerSize',1)
        %-----------------------------------------------------------------%
        plot(x,u2,' xb','LineWidth',2,'MarkerSize',1)
        %-----------------------------------------------------------------%
        plot(x,u3,' xc','LineWidth',2,'MarkerSize',1)
        axis([0 1 -1 1])
        legend('Exact-solution','Lax-Friedrich','Lax-Wendroff','Upwind')
        %-----------------------------------------------------------------%
%         drawnow
        pause(10^(-6)) 
    %     movieVector(k) = getframe(figh, [10 10 520 400]);
    end
hold off 
% myWriter = VideoWriter('Heat Equation');
% myWriter.FrameRate = 20;
% open(myWriter);
% writeVideo(myWriter,movieVector);
% close(myWriter);
else
    disp('Not supported yet !');
end



