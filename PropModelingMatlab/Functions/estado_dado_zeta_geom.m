function [S, COEFF] = estado_dado_zeta_geom(z_tip, pol, P, beta_deg_vec, c_vec_m)
% MEL com GEOMETRIA imposta: beta(ξ) [deg] e c(ξ) [m].

    beta_deg_vec = beta_deg_vec(:); c_vec_m = c_vec_m(:);
    N = numel(P.xi);
    assert(numel(beta_deg_vec)==N && numel(c_vec_m)==N, 'beta/c ≠ tamanho de xi');

    ang = mel_angles_prandtl(z_tip, P);   % φ, F, G, etc.

    a=zeros(N,1); ap=a; W=a; Re=a; Ma=a; Cd=a; Cl=a; alpha=a;
    for k=1:N
        [a(k),ap(k),W(k),Re(k),Ma(k),alpha(k),Cl(k),Cd(k)] = ...
            mel_local_section_geom(beta_deg_vec(k), c_vec_m(k), ...
                                   ang.phi(k), ang.sphi(k), ang.cphi(k), ang.tphi(k), ...
                                   ang.X(k), z_tip, P, pol);
    end

    epsv = max(Cd./max(Cl,P.eps_floor), P.eps_floor);
    [I1p,I2p,J1p,J2p,I1,I2,J1,J2,dxi] = mel_integrals(P, ang, epsv);
    [Tc,Pc,T,Q,CT,CP,ETA] = mel_global_metrics(P, z_tip, I1, I2, J1, J2);
    [dCTdxi,dCPdxi,T_elem] = mel_distributions(P, z_tip, I1p, I2p, J1p, J2p, dxi);

    X = ang.X; F = ang.F; G = ang.G; xi = ang.xi;

    S = struct('xi',xi,'r',P.R*xi,'X',X,'F',F,'zeta',z_tip+0*xi,'phi',ang.phi*180/pi, ...
        'alpha',alpha,'beta',beta_deg_vec,'a',a,'ap',ap,'W',W,'c',c_vec_m, ...
        'Cl',Cl,'Cd',Cd,'eps',epsv,'Re',Re,'Ma',Ma,'G',G,'I1p',I1p,'I2p',I2p, ...
        'J1p',J1p,'J2p',J2p,'I1',I1,'I2',I2,'J1',J1,'J2',J2,'TC',Tc,'PC',Pc, ...
        'T',T,'Q',Q,'CT',CT,'CP',CP,'ETA',ETA,'J',P.J,'zeta_tip',z_tip, ...
        'T_elem',T_elem,'dCTdxi',dCTdxi,'dCPdxi',dCPdxi);

    COEFF = struct('I1',I1,'I2',I2,'J1',J1,'J2',J2);
end
