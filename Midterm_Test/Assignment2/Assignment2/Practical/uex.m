% Exact solution 
function u = uex(x,t,speed,kappa,nameFun)
    switch nameFun
        case 'Advection'
            u = sin(4*pi*(x-speed*t));
        case 'Diffusion'
            u = exp(-16*pi^2*kappa*t)*sin(4*pi*(x-speed*t));
        otherwise
            disp('Not supported yet!!!')
            u = 0;
    end
end