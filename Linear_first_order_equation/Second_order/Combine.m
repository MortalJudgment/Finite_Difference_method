  
% Elliptic Equation on 1 demension
% Laplace Method
%========================================================================================%
%Ho va ten: Nguyen Tu Huy
clc
clear all
close all
format long

disp('===========================--------------------===========================')
disp('Giai bai toan -u_xx(x) = f(x)')
disp('Kem them dieu kien bien cho truoc (Dirichlet or Newmann)')
%Imput value
a=0;
b=1;
%Kiem tra viec nhap dieu kien bai toan
selection=Choose_Boundary_Condition();
if (selection==1) 
    fprintf('Moi nhap gia tri ham u tai %d:  ',a);
    alpha=input('');
    fprintf('Moi nhap gia tri ham u tai %d:  ',b);
    beta=input('');
elseif (selection==2) 
    fprintf('Moi nhap gia tri tai %d: ',a);
    alpha=input('');
    fprintf('Moi nhap gia tri dao ham tai %d:  ',b);
    duN=input('');
elseif (selection==3) 
    fprintf('Moi nhap gia tri dao ham tai %d: ',a);
    du0=input('');
    fprintf('Moi nhap gia tri dao ham tai %d: ',b);
    duN=input('');
else
    selection=Choose_Boundary_Condition();
end
h=zeros(5,1);
errormax=zeros(5,1);
n=4;
for j=1:5
    h(j)=(b-a)/n;
    x=a:h(j):b;
    if selection==1
        y_calcu = solve_Dirichlet(x,h(j),n,alpha,beta);
    elseif selection==2
        y_calcu = solve_Neumann(x,h(j),n,alpha,duN);
    else
        y_calcu = solve_twodifferent(x,h(j),n,du0,duN);
    end
    yex = zeros(n+1,1);
    for i=1:n+1 
        yex(i)=uexact(x(i));
    end
    subplot(2,3,j);
    plot(x,y_calcu,x,yex)
    legend('Discrete solution','Exact solution');
    error = zeros(n+1,1);
    for i=1:n+1
        error(i)=abs(y_calcu(i)-yex(i));
    end
    errormax(j)=max(error);
    n=n*2;
end
subplot(2,3,6)
plot(log(h),2*log(h)+2,log(h),log(errormax))
title('Bai toan hoi tu bac 2')