function y = solve_NeumannDirichlet(x,h,n,du0,beta)
    % Solving poisson problem on 1D with Boudary condition
    % Left Neuman, right Dirichlet
    type = 'NeumannDirichlet';
    y=zeros(n+1,1);
    A=zeros(n+1,n+1);
    B=zeros(n+1,1);
    for i=1:n+1
        if i==1
            A(i,i)=1;
            A(i,i+1)=-1;
        elseif i==n+1
            A(i,i)=h^2;
        else
            A(i,i-1)=-1;
            A(i,i)=2;
            A(i,i+1)=-1;
        end
    end
    alpha = du0;
    for i=1:n+1
        if i==1
            B(i) = 1/2*F(x(i),type)-alpha/h;
        elseif i==n+1
            B(i) = beta;
        else
            B(i) = F(x(i),type);
        end
    end
    y = (A/h^2)\B;
end