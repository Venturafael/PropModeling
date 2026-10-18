function Wc = mel_compute_Wc(z_tip, Cl_vec, P, ang)
% Produto W*c exigido pelo MEL para o Cl alvo.
    Wc = (4*pi*P.lambda .* ang.G .* P.V .* P.R .* z_tip) ./ (Cl_vec * P.B);
end
