clc
close all
clear all
format long
%-------------------------------------------------------------------------%
% Input Data
a = 0; b = 1;
T = 1;
%-------------------------------------------------------------------------%
% Time and Space lines
% Here we consider 1D space depend only on x
Nt = 100;       % Time interval  
k = T/Nt;
Nx = 10;        % Space interval
h = (b-a)/Nx;
%------------------------%
% Let t and x is vector present for time and space.
t = 0:k:T;
x = a:h:b;
%------------------------%
kappa = 1/16;
%------------------------%
uexact = zeros(Nx+1,1);
%------------------------%
% figh = figure;

uk1 = u0(x);
uk2 = u0(x);
uk3 = u0(x);
for j=2:Nt+1
    clf                             
    
    t_j = t(j);
    uk1 = FDM_for_HeatEquation(h,k,Nt,Nx,kappa,0,uk1);      %Forward Euler    
    uk2 = FDM_for_HeatEquation(h,k,Nt,Nx,kappa,1/2,uk2);    %Crank-Nicolson
    uk3 = FDM_for_HeatEquation(h,k,Nt,Nx,kappa,1,uk3);      %Backward Euler
    for i=1:Nx+1
        uexact(i) = uex(t_j,x(i));      % Exactly solution
    end
    
    plot(x,uk1,x,uk2,x,uk3,x,uexact)    % Draw all of them at t = t_k
    xlabel('x')
    ylabel('u')
    title(['Particle at t = ', num2str(t_j),' seconds'])
    axis([0 1 -1 1])
    
%     drawnow
    pause(0.005)
%     movieVector(k) = getframe(figh, [10 10 520 400]);
end

% myWriter = VideoWriter('Heat Equation');
% myWriter.FrameRate = 20;
% open(myWriter);
% writeVideo(myWriter,movieVector);
% close(myWriter);