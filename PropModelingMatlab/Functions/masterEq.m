function res = masterEq(alpha,beta,J,x,sigma,B)
% MASTEREQ  Equação mestra do BEM-Lowry
%   res = 0  →  solução do ângulo de ataque em cada (x,J)

    phi   = beta - alpha;  phi_r = phi * pi/180;
    [CL,CD] = aeroCoeffs(alpha);
    F   = prandtlF(B,x,phi_r);

    lamT = CL*cosd(phi) - CD*sind(phi);
    lamP = CL*sind(phi) + CD*cosd(phi);

    num = pi*x*(4*F*sin(phi_r).^2 - lamT*sigma);
    den =       4*F*sin(phi_r).*cos(phi_r) + lamP*sigma;

    res = num/den - J;
end
