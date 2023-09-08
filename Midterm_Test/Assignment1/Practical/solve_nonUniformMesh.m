function y = solve_nonUniformMesh(x,n,alpha,beta)
    type = 'nonUniformMesh_Dirichlet';
    y=zeros(n+1,1);
    h=zeros(n,1);
    A=zeros(n+1,n+1);
    B=zeros(n+1,1);
    for i=1:n
        h(i)= x(i+1)-x(i);
    end
    for i=1:n+1
        if i==1
            A(i,i)=1;
        elseif i==n+1
            A(i,i)=1;
        else
            A(i,i-1)=-2/(h(i-1)*(h(i)+h(i-1)));
            A(i,i)=2/(h(i-1)*h(i));
            A(i,i+1)=-2/(h(i)*(h(i)+h(i-1)));
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
    y = A\B;
end