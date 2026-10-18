function [dCTdxi, dCPdxi, T_elem] = mel_distributions(P, z_tip, I1p, I2p, J1p, J2p, dxi)
% Distribuições radiais dCT/dξ, dCP/dξ e esforço elemental T(ξ).

    KT   = 0.5*P.rho*P.V^2*pi*P.R^2;
    dTC  = (I1p*z_tip - I2p*z_tip^2);
    T_elem = KT * dTC .* dxi;

    KCT  = (0.5*P.V^2*pi*P.R^2)/(P.n^2*P.D^4);
    KCP  = (0.5*P.V^3*pi*P.R^2)/(P.n^3*P.D^5);

    dCTdxi = KCT * dTC;
    dCPdxi = KCP * (J1p*z_tip + J2p*z_tip^2);
end
