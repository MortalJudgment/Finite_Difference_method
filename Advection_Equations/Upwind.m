    %---------------------------------------------------------------------%
    %                       Upwind Method                                 %
    %---------------------------------------------------------------------%
    % Consider Convection Equations:                                      %
    %           u_t + a*u_x = kappa*u_xx                                  %
    %                                                                     %
    % S.t  a<=x<=b, 0<=t<=T.                                              %
    % u(x,0) = g(x);                                                      %
    % u(t,a) = h1(t);                                                     %
    % u(t,b) = h2(t);                                                     %
    % Stability conditions: k<= 1/(abs(a)/h+2*kappa/h^2)                  %
    %---------------------------------------------------------------------%
function u = Upwind(h,k,Nx,speed,u_before)

    r = speed*k/h;
    %-------------------------------%    
    A = zeros(Nx+1,Nx+1);
    u = zeros(Nx+1,1);
    for i=1:Nx+1
        u(i) = u_before(i);
    end
    %-------------------------------%
    for i=1:Nx+1
        if i==1
            A(i,i) = 1;
        elseif i == Nx+1
            A(i,i) = 1;
        else
            if (speed >=0)
                A(i,i) = 1-r;
                A(i,i-1) = r;
            else
                A(i,i) = 1+r;
                A(i,i+1) = -r;
            end
        end
    end
    u = A*u;   
end
      