function selection=Choose_Boundary_Condition()
selection=0;
while ((selection < 1)||(selection > 4))
    fprintf('Cac loai dieu kien bien: \n');
    fprintf('1. Dieu kien bien Dirichlet.\n');
    fprintf('u0 = alpha; uN = beta \n');
    fprintf('2. Dieu kien bien Neumann.\n');
    fprintf('u0 = alpha, diff_u(uN) = beta \n');
    fprintf('3.Dieu kien bien dao ham tai 2 bien.\n');
    fprintf('diff_u(u0) = du0, diff_u(uN) = duN \n');
    selection=input('Dang dieu kien bien cua bai toan la:  ');
end
end