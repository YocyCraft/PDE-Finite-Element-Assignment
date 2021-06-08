function du=u0_grad(x0,y0)
    du=[(1-2*x0)*sin(pi*y0),x0*(1-x0)*pi*cos(pi*y0)];
end