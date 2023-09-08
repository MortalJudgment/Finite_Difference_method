function f = F(x,type)
    switch type
        case 'Dirichlet'
            f = -12*x^2 - 6*x + 4;
        case 'NeumannDirichlet'
            f = (2*pi)^2*cos(2*pi*x);
        case 'Neumann'
            f = (2*pi)^2*cos(2*pi*x);
        case 'nonUniformMesh_Dirichlet'
            f = -4;
        otherwise
            disp(['Boudary Condition',type,' not yet supported !!'])
            f = 0;
    end
end