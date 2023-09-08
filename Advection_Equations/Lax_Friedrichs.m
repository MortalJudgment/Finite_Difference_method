    %---------------------------------------------------------------------%
    %                       Lax Friedrichs Method                         %
    %---------------------------------------------------------------------%
    % Consider Convection Equations:                                      %
    %           u_t + a*u_x = kappa*u_xx                                  %
    %                                                                     %
    % S.t  a<=x<=b, 0<=t<=T.                                              %
    % u(x,0) = g(x);                                                      %
    % u(t,0) = h(t);                                                      %
    % Stability condition:                                                %
    % If (kappa<= abs(a)*h/2) then k<= 2*kappa/a^2                        %
    % If (kappa> abs(a)*h/2) then k<= h^2/(2*kappa)                       %
    %---------------------------------------------------------------------%
function u = Lax_Friedrichs(h,k,Nx,speed,kappa,u_before)
    s = speed*k/(2*h);
    r = kappa*k/h^2;
    %-------------------------------%    
    A = zeros(Nx+1,Nx+1);
    %-------------------------------%
    for i=1:Nx+1
        if i==1
            A(i,i) = 1;
        elseif i == Nx+1
            A(i,i) = 1;
        else
            A(i,i-1) = s+r;
            A(i,i) = 1-2*r;
            A(i,i+1) = -s+r;
        end
    end
    u1 = A*u_before; 
    u(2:Nx) = u1(2:Nx);
 
end
      