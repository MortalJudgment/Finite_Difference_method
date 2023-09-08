function phi = Poisson_Eq(Nx,dx,y1,y2,u1,v1)
phi = zeros(Nx+1,1);
A   = zeros(Nx+1,Nx+1);
b   = zeros(Nx+1,1);
for i = 1:Nx+1
    if i==1
        A(i,i)      = 1;
        A(i,i+1)    = -1;
    elseif i==Nx+1
        A(i,i-1)    = -1;
        A(i,i)      = 1;
    else
        A(i,i-1)    = -1;
        A(i,i)      = 2;
        A(i,i+1)    = -1;
    end
end
for i = 1:Nx+1
    if i==1
        b(i)    = 1/2*( y1*u1(i) + y2*v1(i) );
    elseif i==Nx+1
        b(i)    = 1/2*( y1*u1(i) + y2*v1(i) );
    else
        b(i)    = y1*u1(i) + y2*v1(i);
    end
end

epsilon11 = 0.0000001;

    phi = dx^2.*(A+epsilon11*eye(Nx+1,Nx+1))^(-1)*b;
