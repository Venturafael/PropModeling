function [Re, Ma, Vo] = reynolds_mach_local(r_in, c_in, J, RPM, R_in, rho, mu, Vson)
% r_in, c_in, R_in em polegadas; rho [kg/m^3], mu [Pa·s], Vson [m/s].
    inch = 0.0254;
    r  = r_in*inch;  c = c_in*inch;  D = 2*R_in*inch;
    n  = RPM/60;     Omega = 2*pi*n;
    Vinf = J*n*D;    Vtan  = Omega*r;   Vo = hypot(Vtan, Vinf);
    Re = rho*Vo*c/mu;   Ma = Vo/Vson;
end
