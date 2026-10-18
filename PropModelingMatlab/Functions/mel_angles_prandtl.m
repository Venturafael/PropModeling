function ang = mel_angles_prandtl(z_tip, P)
% Prandtl F, ângulos e G(ξ) para o MEL.

    xi   = P.xi;  X = P.X;
    phi_t  = atan( P.lambda*(1 + 0.5*z_tip) );
    sphi_t = max(sin(phi_t), 1e-12);

    ft = (P.B/2).*(1 - xi)./sphi_t;
    F  = (2/pi).*acos( exp(-max(ft,0)) );
    F  = max(F, P.F_min);

    phi  = atan( tan(phi_t)./max(xi,1e-8) );
    sphi = max(sin(phi), P.sphi_floor);
    cphi = cos(phi);
    tphi = tan(phi);

    G = F .* X .* cphi .* sphi;

    ang = struct('xi',xi,'X',X,'phi',phi,'sphi',sphi,'cphi',cphi,'tphi',tphi,'F',F,'G',G);
end
