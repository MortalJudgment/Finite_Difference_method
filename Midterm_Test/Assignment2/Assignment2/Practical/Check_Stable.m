% Check stable condition
function r = Check_Stable(type,nameFun,nameSche,k,h,speed,kappa)
% type: type of method (e.g: Forward-Euler or Crank-Nicolson)
% nameFun: name of the function (e.g Advection or Diffusion)
% nameSche: name of the Schemes (e.g Lax-Friedrich, Lax-Wendroff or Upwind)
% k: time step
% h: space step
% speed, kappa: const of function
% ------------------------------------------------------------------------%
% function return to 0(false), 1(true)
% Means 0 -> not sastify stable condition
%       1 -> sastify stable condition
if strcmp(type,'Crank-Nicolson')== 1
    disp('There is no stable condition at all')
    r = 1;
elseif strcmp(type,'Forward-Euler')== 1
    switch nameFun
        case 'Advection'
            if strcmp(nameSche,'Lax-Friedrich')== 1
                if (speed*k/h) > 1
                    r = 0;
                else r = 1;
                end
            elseif strcmp(nameSche,'Lax-Wendroff')==1
                if (speed*k/h) > 1
                    r = 0;
                else r = 1;
                end
            elseif strcmp(nameSche,'Upwind')==1
                if (speed*k/h) > 1
                    r = 0;
                else r = 1;
                end
            else
                disp('Not supported yet!!')
                r = 2;
            end
        case 'Diffusion'
            if strcmp(nameSche,'Lax-Friedrich')==1
                if kappa >= (speed*h)/2
                    if k <= (h^2/(2*kappa))
                        r = 1;
                    else
                        r = 0;
                    end
                else
                    if k <= ((2*kappa)/speed^2)
                        r = 1;
                    else
                        r = 0;
                    end
                end
            elseif strcmp(nameSche,'Lax-Wendroff')==1
                if kappa >= 1/2*(speed*h - speed^2*k)
                    if k <= (-kappa+ sqrt(kappa^2+speed^2*h^2))/speed^2
                        r = 1;
                    else
                        r = 0;
                    end
                else
                    if (kappa*k)/h^2 > 0
                        r = 1;
                    else
                        r = 0;
                    end
                end
            elseif strcmp(nameSche,'Upwind')==1
                if kappa > 0
                    if k < (h^2/(speed*h+2*kappa))
                        r = 1;
                    else
                        r = 0;
                    end
                else
                    if k <= (speed*h+2*kappa)/speed^2
                        r = 1;
                    else
                        r = 0;
                    end
                end
            else
                disp('Not supported yet!!')
                r = 2;
            end
        otherwise
            disp('Not supported yet!!')
            r = 2;
    end
else
    disp('Not supported yet!!')
    r = 2;
end
