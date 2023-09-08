function uex=uexact(x,type)
    switch type
        case 'Dirichlet'
            uex = x^4 + x^3 - 2*x^2 + 4 ;
        case 'NeumannDirichlet'
            uex = cos(2*pi*x) + 2*x;
        case 'Neumann'
            uex = cos(2*pi*x) ;
        case 'nonUniformMesh_Dirichlet'
            uex = 2*x*(x-1); 
        otherwise
            disp(['Exactly Solution for Boudary Condition',type,' not yet supported !!'])
            uex = 0;
    end
end