% Lax-Friedrich schemes for advection equations
%             u_t + speed*u_x = 0
%       0 < t < T , a < x < b.
% where speed > 0
%----------------------------------------------%
function u = Lax_Friedrich_advection(h,k,Nx,speed,u_before,type)
% where type is considering for using Forward Euler(FE) or
% Crank-Nicolson(CN)
s = (speed*k)/h;
I = eye(Nx+1,Nx+1);
switch type
    case 'Forward-Euler'
        for i=2:Nx
            u(i) = (u_before(i-1)+u_before(i+1))/2.0 - s/2.0*(-u_before(i-1)+u_before(i+1));
        end
    case 'Crank-Nicolson'
        A = zeros(Nx+1,Nx+1);
        for i=1:Nx+1
            if i==1
                A(i,i) = 1;
            elseif i==Nx+1
                A(i,i) = 1;
            else
                A(i,i-1) = 1/4*s;
                A(i,i) = 1;
                A(i,i+1) = -1/4*s;
            end
        end
        B = zeros(Nx+1,Nx+1);
        for i=1:Nx+1
            if i==1
                B(i,i) = 1;
            elseif i==Nx+1
                B(i,i) = 1;
            else
                B(i,i-1) = -1/4*s;
                B(i,i) = 1;
                B(i,i+1) = 1/4*s;
            end
        end
        u = B^(-1)*A*u_before;
    otherwise
        disp('Non-available');
        u=0;
end