    %Finite Different Method for Heat Equation
    % u_t = 1/16*u_xx
    % Consider 0<=x<=1, 0<=t<=T.
    % u(0,x) = sin(2*pi*x);
    % u(t,0) = u(t,1) =0;
    % Stability condition: k/(h^2) <= 8.
function u = FDM_for_HeatEquation(h,k,Nt,Nx,kappa,theta,u_before)
    r = kappa*k/h^2;
    
    I = eye(Nx+1);
    A = zeros(Nx+1,Nx+1);
    u = zeros(Nx+1,1);
    
    for i=1:Nx+1
        u(i) = u_before(i);
    end
    for i=2:Nx
        A(i,i-1) = -r;
        A(i,i) = 2*r;
        A(i,i+1) = -r;
    end
    u = (I+theta*A)^(-1)*(I-(1-theta)*A)*u;
end
      