% Elliptic Equation on 1 demension
% Laplace Method
%========================================================================================%
%Ho va ten: Nguyen Tu Huy
clc
clear all
close all
format long

disp('===============---------------------------------------------===============')
disp('--------------- Solving differece equations -u_xx(x) = f(x) ---------------')
disp('===============---------------------------------------------===============')
disp('On (0,1) With Boundary condition')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
%Imput value
a=0;
b=1;

%-------------------------------------------------------------------------%
%------------------Dirichlet conditions on uniform mesh.------------------%
%-------------------------------------------------------------------------%
% Consider differece equations:                                           %
%              - u_xx  = -12.x^2 - 6.x + 4                                %
%   u(0) = 4, u(1) = 4.                                                   %
%-------------------------------------------------------------------------%
h=zeros(5,1);
errormax=zeros(5,1);
n=4;
type = 'Dirichlet';
figure
for j=1:5
    h(j)=(b-a)/n;
    x=a:h(j):b;
    u1 = solve_Dirichlet(x,h(j),n,4,4);
    yex = zeros(n+1,1);
    for i=1:n+1
        yex(i)=uexact(x(i),type);
    end
    subplot(2,3,j);
    plot(x,u1,x,yex)
    legend('Discrete solution','Exact solution');
    title('Dirichlet condition on Uniform Mesh');
    error = zeros(n+1,1);
    for i=1:n+1
        error(i)=abs(u1(i)-yex(i));
    end
    errormax(j)=max(error);
    n=n*2;
end
subplot(2,3,6)
plot(log(h),2*log(h)+2,log(h),log(errormax))
title('Bai toan hoi tu bac 2')

%-------------------------------------------------------------------------%
%--------------Neumann Dirichlet conditions on uniform mesh.--------------%
%-------------------------------------------------------------------------%
% Consider differece equations:                                           %
%              - u_xx  = 4*pi^2*cos(2*pi*x)                               %
%   left Neuman, right Dirichlet                                          %
%   u'(0) = 2, u(1) = 3.                                                  %
%-------------------------------------------------------------------------%
h=zeros(5,1);
errormax=zeros(5,1);
n=4;
type = 'NeumannDirichlet';
figure
for j=1:5
    h(j)=(b-a)/n;
    x=a:h(j):b;
    u2 = solve_NeumannDirichlet(x,h(j),n,2,3);
    yex = zeros(n+1,1);
    for i=1:n+1
        yex(i)=uexact(x(i),type);
    end
    subplot(2,3,j);
    plot(x,u2,x,yex)
    legend('Discrete solution','Exact solution');
    title('Neumann Dirichlet condition');
    error = zeros(n+1,1);
    for i=1:n+1
        error(i)=abs(u2(i)-yex(i));
    end
    errormax(j)=max(error);
    n=n*2;
end
subplot(2,3,6)
plot(log(h),2*log(h),log(h),log(errormax))
title('Bai toan hoi tu bac 2')

%-------------------------------------------------------------------------%
%-------------------Neumann conditions on uniform mesh.-------------------%
%-------------------------------------------------------------------------%
% Consider differece equations:                                           %
%               - u_xx  = 4*pi^2*cos(2*pi*x).                             %
%   u'(0) = 2*pi, u'(1) = 2*pi.                                           %
%   int(f,0,1) = u'(1) - u'(0);     int(u,0,1) = 0;                       %
%-------------------------------------------------------------------------%
h=zeros(5,1);
errormax=zeros(5,1);
n=4;
type = 'Neumann';
figure
for j=1:5
    h(j)=(b-a)/n;
    x=a:h(j):b;
    u3 = solve_Neumann(x,h(j),n,0,0);
    yex = zeros(n+1,1);
    for i=1:n+1
        yex(i)=uexact(x(i),type);
    end
    subplot(2,3,j);
    plot(x,u3,x,yex)
    legend('Discrete solution','Exact solution');
    title('Neumann condition');
    error = zeros(n+1,1);
    for i=1:n+1
        error(i)=abs(u3(i)-yex(i));
    end
    errormax(j)=max(error);
    n=n*6;
end
subplot(2,3,6)
plot(log(h),1.5*log(h)+1,log(h),log(errormax))
title('Bai toan hoi tu bac 1.5')



%-------------------------------------------------------------------------%
%----------------Dirichlet conditions on non-uniform mesh.----------------%
%-------------------------------------------------------------------------%
% Consider differece equations:                                           %
%               - u_xx  = 4.                                              %
%   u(0) = 0, u(1) = 0.                                                   %
%   x = 1-cos((pi*i)/(2*N))  forall i=0,...N                              %
%-------------------------------------------------------------------------%
h=zeros(5,1);
errormax=zeros(5,1);
n=4;
type = 'nonUniformMesh_Dirichlet';
figure
for j=1:5
    h(j)=(b-a)/n;
    x=zeros(n+1,1);
    for i=1:n+1
        x(i) = 1-cos((pi*i)/(2*(n+1)));
    end
    u4 = solve_nonUniformMesh(x,n,0,0);
    yex = zeros(n+1,1);
    for i=1:n+1
        yex(i)=uexact(x(i),type);
    end
    subplot(2,3,j);
    plot(x,u4,x,yex)
    legend('Discrete solution','Exact solution');
    title('Dirichlet condition on non-Uniform Mesh');
    error = zeros(n+1,1);
    for i=1:n+1
        error(i)=abs(u4(i)-yex(i));
    end
    errormax(j)=max(error);
    n=n*2;
end
subplot(2,3,6)
plot(log(h),2*log(h)+2,log(h),log(errormax))
title('Bai toan hoi tu bac gan nhu 2')
