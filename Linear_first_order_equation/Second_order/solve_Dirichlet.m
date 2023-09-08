function y_discrete = solve_Dirichlet(x,h,n,alpha,beta)
    y=zeros(n-1,1);
    A=zeros(n-1,n-1);
    B=zeros(n-1,1);
    for i=1:n-1
        if i==1
            A(i,i)=2;
            A(i,i+1)=-1;
        elseif i==n-1
            A(i,i-1)=-1;
            A(i,i)=2;
        else
            A(i,i-1)=-1;
            A(i,i)=2;
            A(i,i+1)=-1;
        end
    end
    for i=1:n-1
        if i==1
            B(i)=F(x(i+1))+alpha/(h^2);
        elseif i==n-1
            B(i)=F(x(i+1))+beta/(h^2);
        else
            B(i)=F(x(i+1));
        end
    end
    y = A\B*h^2;
    y_discrete = [alpha;y;beta];
end