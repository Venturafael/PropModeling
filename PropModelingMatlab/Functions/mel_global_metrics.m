function [Tc, Pc, T, Q, CT, CP, ETA] = mel_global_metrics(P, z_tip, I1, I2, J1, J2)
% Coeficientes globais e métricas (Tc, Pc, T, Q, CT, CP, η).

    Tc = I1*z_tip - I2*z_tip^2;
    Pc = J1*z_tip + J2*z_tip^2;

    T   = 0.5*P.rho*P.V^2*pi*P.R^2 * Tc;
    Pae = 0.5*P.rho*P.V^3*pi*P.R^2 * Pc;
    Q   = Pae / (2*pi*P.n);

    CT  = T /(P.rho*P.n^2*P.D^4);
    CP  = (2*pi*Q)/(P.rho*P.n^2*P.D^5);
    ETA = P.J*CT/max(CP,eps);
end
