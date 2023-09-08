%----------------------------------------------%
% Lax-Wendroff schemes for diffusion equations
%             u_t + speed*u_x = kappa*u_xx
%       0 < t < T , a < x < b.
% where speed > 0
%----------------------------------------------%
function u = Lax_Wendroff_diffusion(h,k,Nx,speed,kappa,u_before,type)
% where type is considering for using Forward Euler(FE) or
% Crank-Nicolson(CN)
s = (speed*k)/h;
r = (kappa*k)/h^2;
A = zeros(Nx+1,Nx+1);
I = eye(Nx+1,Nx+1);
for i=2:Nx
    A(i,i-1) = 1/2*s + r;
    A(i,i) = - 2*r;
    A(i,i+1) = - 1/2*s + r;
end
switch type
    case 'Forward-Euler'
        u = (I+A)*u_before;
    case 'Crank-Nicolson'
        u = (I-1/2*A)^(-1)*(I+1/2*A)*u_before;
    otherwise
        disp('Non-available');
        u=0;
end