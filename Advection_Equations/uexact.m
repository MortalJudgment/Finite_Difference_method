function u = uexact(t,x,h,k,speed)
    u = exp(-sqrt(h^2/(2*k))*pi^2*t)*sin(2*pi*(x-speed*t));
end