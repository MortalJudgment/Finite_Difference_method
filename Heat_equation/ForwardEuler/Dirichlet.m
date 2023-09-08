% Finite Volume method for Diffusion Equation on 1D
% Consider on Domain Omega = (0 1), t>0
% with initial condition u0(x)
% and Dirichlet boudary condition u(0,t)=g1(t); u(1,t)=g2(t)
function u = Dirichlet(beta,N,dt,x,u1)
%-------------------------------------------------------------------------%
%-----  Solve 1D Heat equation:  -----%
%                 u_t = beta*u_xx   ,in [a,b]
%-------------------------------------------------------------------------%
% Unit matrix
I = eye(N+1,N+1);
% Creare the Matrix
A=zeros(N+1,N+1);
for ii=2:N
    r = 1.0/(x(ii)-x(ii-1))^2;
    A(ii,ii-1)  = r;
    A(ii,ii)    = -2*r;
    A(ii,ii+1)  = r;
end
A = beta*dt.*A;
u = (I+A)*u1;