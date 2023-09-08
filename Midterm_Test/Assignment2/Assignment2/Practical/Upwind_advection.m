%----------------------------------------------%
% Upwind schemes for advection equations
%             u_t + speed*u_x = 0
%       0 < t < T , a < x < b.
% where speed > 0
%----------------------------------------------%
% where type is considering for using Forward Euler(FE) or
% Crank-Nicolson(CN)
function u = Upwind_advection(h,k,Nx,speed,u_before,type)
s = (speed*k)/h;
A = zeros(Nx+1,Nx+1);
I = eye(Nx+1,Nx+1);
for i=2:Nx
    A(i,i-1) = -s;
    A(i,i) = s;
end
switch type
    case 'Forward-Euler'
        u = (I-A)*u_before;
    case 'Crank-Nicolson'
        u = (I+1/2*A)^(-1)*(I-1/2*A)*u_before;
    otherwise
        disp('Non-available');
        u=0;
end