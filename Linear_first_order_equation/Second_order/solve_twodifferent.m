function y_discrete = solve_twodifferent(x,h,n,du0,duN)
    M=h/2*(F(0)+F(1))-(duN-du0);
    for i=1:n-1
        M=M+h*F(x(i+1));
    end
    f1=zeros(n+1,1);
    for i=1:n+1
        f1(i)=F(x(i))-M;
    end
    y=zeros(n+1,1);
    A=zeros(n+1,n+1);
    B=zeros(n+1,1);
    for i=1:n+1
        if i==1
            A(i,i)=1;
            A(i,i+1)=-1;
        elseif i==n+1
            A(i,i-1)=-1;
            A(i,i)=1;
        else
            A(i,i-1)=-1;
            A(i,i)=2;
            A(i,i+1)=-1;
        end
    end
    A=A/(h^2);
    for i=1:n+1
        if i==1
            B(i)=1/2*f1(i)+du0/h;
        elseif i==n+1
            B(i)=1/2*f1(i)+duN/h;
        else
            B(i)=f1(i);
        end
    end
    epsilon11=0.0000001/h^2;
    y_discrete = (A+epsilon11*eye(n+1,n+1))\B;
end