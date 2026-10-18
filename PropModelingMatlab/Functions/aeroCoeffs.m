function [CL,CD] = aeroCoeffs(alpha)
% AEROCOEFFS  Polares analíticas simplificadas de Lowry
%   alpha em graus

    C1 = 0.09458 * (alpha + 4.8);
    CL = (1.5 + C1 - sqrt((1.5 - C1).^2 + 0.04)) / 2;
    CD = 0.015 + 0.03 * ((4*C1/3) - 1).^4;
end
