% Elliptic Equation on 1 demension
% Laplace Method
clc
clear all
close all
format long

%Imput value
a=0;
b=1;
u0=2;
uN=2;
%Set the step
h=zeros(6,1);
errormax=zeros(6,1);
n = 6;
for j=1:6
    n=2*n;
    %Calculate deta_x
    h(j)=(b-a)/n;
    x=a:h(j):b;
    y=zeros(n-1,1);
    A=zeros(n-1,n-1);
    B=zeros(n-1,1);
    for i=1:n-1
        if i==1
            A(i,i)=-20;
            A(i,i+1)=6;
            A(i,i+2)=4;
            A(i,i+3)=-1;
        elseif i==2
            A(i,i-1)=16;
            A(i,i)=-30;
            A(i,i+1)=16;
            A(i,i+2)=-1;
        elseif i==n-2
            A(i,i-2)=-1;
            A(i,i-1)=16;
            A(i,i)=-30;
            A(i,i+1)=16;
        elseif i==n-1
            A(i,i-3)=-1;
            A(i,i-2)=4;
            A(i,i-1)=6;
            A(i,i)=-20;
        else
            A(i,i-2)=-1;
            A(i,i-1)=16;
            A(i,i)=-30;
            A(i,i+1)=16;
            A(i,i+2)=-1;
        end
    end
    for i=1:n-1
        if i==1
            B(i)=F(x(i+1))+11*u0/(12*h(j)^2);
        elseif i==2
            B(i)=F(x(i+1))-u0/(12*h(j)^2);
        elseif i==n-2
            B(i)=F(x(i+1))-uN/(12*h(j)^2);
        elseif i==n-1
            B(i)=F(x(i+1))+11*uN/(12*h(j)^2);
        else
            B(i)=F(x(i+1));
        end
    end
    y = (12*h(j)^2)*((-A)^-1)*B;
    y_calcu = [u0;y;uN];
    yex = zeros(n+1,1);
    for i=1:n+1
        yex(i)=uexact(x(i));
    end
%     subplot(2,4,j);
%     plot(x,y_calcu,x,yex)
%     legend('Discrete solution','Exact solution');
    error = zeros(n+1,1);
    for i=1:n+1
        error(i)=abs(y_calcu(i)-yex(i));
    end
    errormax(j)=max(error);
end
figure
plot(log(h),4*log(h)+8,log(h),log(errormax))
            