clc
close all
clear all
format long
a = 0; b = 1;
T = 0.01;
%-----------------------------------%
Nt = 100;
k = T/Nt;
Nx = 30;
h = (b-a)/Nx;
%-----------------------------------%
t = 0:k:T;
x = (a:h:b);
%-----------------------------------%
speed = 1;
kappa = h^2/(2*k);
%--------- Lax Friedrichs ----------%
if ((kappa<=speed*h/2)&&(k<=2*kappa/speed^2))||((kappa>speed*h/2)&&(k<=h^2/(2*kappa)))
    u1 = (u0(x));
    for i=2:Nt+1
        u1 = Lax_Friedrichs(h,k,Nx,speed,kappa,(u1)');
        u1(1) = uexact(t(i),x(1),h,k,speed);
        u1(Nx+1) = uexact(t(i),x(Nx+1),h,k,speed);
    end
    hold on
    plot(x,u1,'-. g','LineWidth',2,'MarkerSize',2)
else
    disp('Error!');
end

%---------- Exact Solution ----------%
y = 0:0.01:1;
plot(x,uexact(T,x,h,k,speed),'- r','LineWidth',2,'MarkerSize',3)
hold off
