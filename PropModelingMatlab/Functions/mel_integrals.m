function [I1p,I2p,J1p,J2p,I1,I2,J1,J2,dxi] = mel_integrals(P, ang, epsv)
% Termos diferenciais e integrais do MEL (Adkins–Liebeck).

    xi = P.xi;  sphi = ang.sphi;  cphi = ang.cphi;  tphi = ang.tphi;  G = ang.G;

    I1p = 4 .* xi .* G .* (1 - epsv.*tphi);
    I2p = P.lambda .* (I1p ./ max(2*xi,1e-12)) .* (1 + epsv./max(tphi,1e-12)) .* sphi .* cphi;
    J1p = 4 .* xi .* G .* (1 + epsv./max(tphi,1e-12));
    J2p = 0.5 .* J1p .* (1 - epsv.*tphi) .* (cphi.^2);

    dxi = gradient(xi);
    I1  = sum(I1p .* dxi, 'omitnan');
    I2  = sum(I2p .* dxi, 'omitnan');
    J1  = sum(J1p .* dxi, 'omitnan');
    J2  = sum(J2p .* dxi, 'omitnan');
end
