%----------------------------------------------%
% Lax-Wendroff schemes for advection equations
%             u_t + speed*u_x = 0
%       0 < t < T , a < x < b.
% where speed > 0
% ---------------------------------------------%
function u = Lax_Wendroff_advection(h,k,Nx,speed,u_before,type)
% where type is considering for using Forward Euler(FE) or
% Crank-Nicolson(CN)
s = (speed*k)/h;
A = zeros(Nx+1,Nx+1);
I = eye(Nx+1,Nx+1);
for i=2:Nx
    A(i,i-1) = -1/2*s - 1/2*s^2;
    A(i,i) = s^2;
    A(i,i+1) = 1/2*s - 1/2*s^2;
end
switch type
    case 'Forward-Euler'
        u =(I-A)*u_before;
    case 'Crank-Nicolson'
        u = (I+1/2*A)^(-1)*(I-1/2*A)*u_before;
    otherwise
        disp('Non-available!!');
        u=0;
end