function y = solve_Dirichlet(x,h,n,alpha,beta)
    type = 'Dirichlet';
    y=zeros(n+1,1);
    A=zeros(n+1,n+1);
    B=zeros(n+1,1);
    for i=1:n+1
        if i==1
            A(i,i)=h^2;
        elseif i==n+1
            A(i,i)=h^2;
        else
            A(i,i-1)=-1;
            A(i,i)=2;
            A(i,i+1)=-1;
        end
    end
    for i=1:n+1
        if i==1
            B(i) = alpha;
        elseif i==n+1
            B(i) = beta;
        else
            B(i) = F(x(i),type);
        end
    end
    y = (A/h^2)\B;
end