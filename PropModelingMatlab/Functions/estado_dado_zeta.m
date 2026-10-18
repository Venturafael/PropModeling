function [S, COEFF] = estado_dado_zeta(z_tip, Cl_here, pol, P)
% Estado aerodinâmico pelo MEL (Adkins & Liebeck) para um ζ dado.
% Mantém compatibilidade com sua versão original.

    if isempty(Cl_here), error('Cl_here está vazio.'); end
    N = numel(P.xi);
    if isscalar(Cl_here), Cl = Cl_here*ones(N,1);
    else
        assert(numel(Cl_here)==N,'Cl_here deve ter N=%d elementos.',N);
        Cl = Cl_here(:);
    end

    % 1) Ângulos, Prandtl e G
    ang = mel_angles_prandtl(z_tip, P);   % struct: xi,X,phi,sphi,cphi,tphi,F,G

    % 2) Produto W*c (MEL)
    Wc = mel_compute_Wc(z_tip, Cl, P, ang);

    % 3) Loop por estação: resolve a, a', W, c, α, Cd, Re, Ma
    a   = zeros(N,1);  ap  = zeros(N,1);
    W   = zeros(N,1);  c   = zeros(N,1);
    Re  = zeros(N,1);  Ma  = zeros(N,1);
    Cd  = zeros(N,1);  alpha = zeros(N,1);

    for k = 1:N
        [a(k), ap(k), W(k), c(k), Re(k), Ma(k), alpha(k), Cd(k)] = ...
            mel_local_section( ...
                Cl(k), Wc(k), ang.sphi(k), ang.cphi(k), ang.tphi(k), ...
                ang.X(k), z_tip, P, pol);
    end

    % 4) Integrais do MEL e globais
    epsv = max(Cd./max(Cl,P.eps_floor), P.eps_floor);
    [I1p,I2p,J1p,J2p,I1,I2,J1,J2,dxi] = mel_integrals(P, ang, epsv);

    [Tc, Pc, T, Q, CT, CP, ETA] = mel_global_metrics(P, z_tip, I1, I2, J1, J2);

    % 5) Distribuições (por ξ)
    [dCTdxi, dCPdxi, T_elem] = mel_distributions(P, z_tip, I1p, I2p, J1p, J2p, dxi);

    beta = ang.phi*180/pi + alpha;

    % ===== Monta struct de saída =====
    S = struct('xi',P.xi,'r',P.R*P.xi,'X',ang.X,'F',ang.F,'zeta',z_tip+0*P.xi, ...
        'phi',ang.phi*180/pi,'alpha',alpha,'beta',beta,'a',a,'ap',ap, ...
        'W',W,'c',c,'Cl',Cl,'Cd',Cd,'eps',epsv,'Re',Re,'Ma',Ma,'G',ang.G, ...
        'I1p',I1p,'I2p',I2p,'J1p',J1p,'J2p',J2p, ...
        'I1',I1,'I2',I2,'J1',J1,'J2',J2, ...
        'TC',Tc,'PC',Pc,'T',T,'Q',Q,'CT',CT,'CP',CP,'ETA',ETA,'J',P.J, ...
        'zeta_tip',z_tip,'T_elem',T_elem,'dCTdxi',dCTdxi,'dCPdxi',dCPdxi);

    COEFF = struct('I1',I1,'I2',I2,'J1',J1,'J2',J2);
end
