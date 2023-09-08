 function u = Nernst_Planck(Nx,dx,dt,d,v,g1,g2,u1,v1,phi1)
u = zeros(Nx+1,1);
A = zeros(Nx+1,Nx+1);
for i=1:Nx+1
    if i==1
        A(i,i)      = 1 - (2*d*dt/(dx^2)) + (v*dt/(dx^2))*(2*phi1(i+1)-2*phi1(i))...
                        + (g1*dt/(dx^2))*(2*u1(i+1)-2*u1(i)) + (g2*dt/(dx^2))*(2*v1(i+1)-2*v1(i));
        A(i,i+1)    = (2*d*dt/(dx^2));
    elseif i==Nx+1
        A(i,i-1)    = (2*d*dt/(dx^2));
        A(i,i)      = 1 - (2*d*dt/(dx^2)) + (v*dt/(dx^2))*(2*phi1(i-1)-2*phi1(i))...
                        + (g1*dt/(dx^2))*(2*u1(i-1)-2*u1(i)) + (g2*dt/(dx^2))*(2*v1(i-1)-2*v1(i));
    else
        A(i,i-1)    = (d*dt/(dx^2)) - (v*dt/(dx^2))*1/2*1/2*(phi1(i+1)-phi1(i-1))...
                        - (g1*dt/(dx^2))*1/2*1/2*(u1(i+1)-u1(i-1)) - (g2*dt/(dx^2))*1/2*1/2*(v1(i+1)-v1(i-1));
        A(i,i)      = 1 - (2*d*dt/(dx^2)) + (v*dt/(dx^2))*(phi1(i+1)-2*phi1(i)+phi1(i-1))...
                        + (g1*dt/(dx^2))*(u1(i+1)-2*u1(i)+u1(i-1)) + (g2*dt/(dx^2))*(v1(i+1)-2*v1(i)+v1(i-1));
        A(i,i+1)    = (d*dt/(dx^2)) + (v*dt/(dx^2))*1/2*1/2*(phi1(i+1)-phi1(i-1))...
                        + (g1*dt/(dx^2))*1/2*1/2*(u1(i+1)-u1(i-1)) + (g2*dt/(dx^2))*1/2*1/2*(v1(i+1)-v1(i-1));
    end
end
u = A*u1;


