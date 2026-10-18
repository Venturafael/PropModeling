function F = prandtlF(B,x,phi)
% PRANDTLF  Fator de ponta de Prandtl (Lowry)
%   F = prandtlF(B,x,phi)
%   B   – número de pás
%   x   – estação radial (r/R)
%   phi – ângulo de fluxo em radianos

    f = B/2 * (1-x) ./ (x .* sin(phi));
    F = max( 2/pi * acos(exp(-f)) , 0.05 );  % limite inferior p/ evitar zero
end
