function y_discrete = solve_NeumannDirichlet(x,h,n,alpha,duN)
    y=zeros(n,1);
    A=zeros(n,n);
    B=zeros(n,1);
    for i=1:n
        if i==1
            A(i,i)=2;
            A(i,i+1)=-1;
        elseif i==n
            A(i,i-1)=-1;
            A(i,i)=1;
        else
            A(i,i-1)=-1;
            A(i,i)=2;
            A(i,i+1)=-1;
        end
    end
    beta=duN;
    for i=1:n
        if i==1
            B(i)=1/2*F(x(i+1))+alpha/(h^2);
        elseif i==n
            B(i)=1/2*F(x(i+1))+beta/h;
        else
            B(i)=F(x(i+1));
        end
    end
    y = (A/h^2)\B;
    y_discrete = [alpha;y];
end