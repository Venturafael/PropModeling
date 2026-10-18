function [a_out, ap_out, W_out, c_out, Re_out, Ma_out, alpha_out, Cd_out] = ...
    mel_local_section(Cl_target, Wc_k, sphi_k, cphi_k, tphi_k, Xk, z_tip, P, pol)
% Resolve a, a', W, c, α, Cd, Re, Ma em uma estação (iteração em a).

    a_k = 0.05;   % chute
    for it = 1:100
        W_k = P.V*(1 + a_k)/sphi_k;
        c_k = Wc_k / max(W_k, 1e-12);
        if P.c_min_tip > 0, c_k = max(c_k, P.c_min_tip); end

        Re_k = W_k*c_k / max(P.nu,1e-12);
        Ma_k = W_k / P.a_sound;

        alpha_k = P.alpha_bias_deg + alpha_for_Cl_target(Cl_target, Re_k, Ma_k, pol);
        [~, Cd_k] = lookup_clarky2(alpha_k, Re_k, Ma_k, pol, 'interp');
        Cd_k = P.Cd_factor * Cd_k;
        eps_k = max(Cd_k / max(Cl_target, P.eps_floor), P.eps_floor);

        a_new = 0.5*z_tip*(cphi_k^2)*(1 - eps_k*tphi_k);
        ap_k  = 0.5*z_tip*(cphi_k*sphi_k)*(1 + eps_k/max(tphi_k,1e-12))/max(Xk,1e-12);

        a_new = min(max(a_new, P.a_limits(1)),   P.a_limits(2));
        ap_k  = min(max(ap_k , P.ap_limits(1)),  P.ap_limits(2));

        a_upd = 0.7*a_k + 0.3*a_new;
        if abs(a_upd - a_k) < P.tol_loc
            a_k = a_upd;
            a_out = a_k; ap_out = ap_k;
            W_out = W_k; c_out  = c_k;
            Re_out = Re_k; Ma_out = Ma_k;
            alpha_out = alpha_k; Cd_out = Cd_k;
            return
        end
        a_k = a_upd;
    end

    % se não convergir em 100 iterações, retorna último estado
    a_out = a_k; ap_out = ap_k;
    W_out = W_k; c_out  = c_k;
    Re_out = Re_k; Ma_out = Ma_k;
    alpha_out = alpha_k; Cd_out = Cd_k;
end
